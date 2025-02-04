#!/usr/bin/env python3
import argparse
import asyncio
import glob
import logging
import os
import shlex
import shutil
import signal
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Set, TextIO, Tuple

import coloredlogs
import ruamel.yaml
import ruamel.yaml.error

_log = logging.getLogger("AdapterID")


def quote(value):
    if isinstance(value, os.PathLike):
        value = os.fspath(value)

    return shlex.quote(str(value))


def determine_ar_executable(candidate: Optional[str]) -> Optional[str]:
    if candidate is None:
        known_executables = ("adapterremoval3", "AdapterRemoval")
        for executable in known_executables:
            executable = shutil.which(executable)
            if executable is not None:
                break
        else:
            _log.error("Could not find %s executable", "/".join(known_executables))
            executable = None
    else:
        executable = shutil.which(candidate)
        if executable is None:
            _log.error("Could not find %s executable", candidate)

    return executable


def filter_options(dd: Dict[Any, Any]) -> Iterator[Tuple[Any, Any]]:
    for key, value in dd.items():
        if key not in ("Options", "Prefixes", "Genomes"):
            yield key, value


def safe_load(stream: TextIO):
    yaml: Any = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)  # mypy: ignore

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ruamel.yaml.error.MantissaNoDotYAML1_1Warning)

        data: Any = yaml.load(stream)

    for _, samples in filter_options(data):
        for sample, libraries in filter_options(samples):
            for library, lanes in filter_options(libraries):
                for lane, tmpl in filter_options(lanes):
                    yield sample, library, lane, tmpl


@dataclass
class Command:
    dry_run: bool
    key: Tuple[str, ...]
    args: Tuple[str, ...]
    destination: str

    def __hash__(self):
        return hash(self.key)


def build_command(
    args: argparse.Namespace,
    sample: str,
    library: str,
    lane: str,
    tmpl: Any,
) -> Iterator[Tuple[bool, Optional[Command]]]:
    if not isinstance(tmpl, str):
        _log.warning("Skipping pre-processed data: %s > %s > %s", sample, library, lane)
        yield False, None
        return

    files_1 = sorted(glob.glob(tmpl.replace("{Pair}", "1")))
    files_2 = sorted(glob.glob(tmpl.replace("{Pair}", "2")))
    if not (files_1 or files_2):
        _log.error("No files found for %r", tmpl)
        return
    elif len(files_1) != len(files_2):
        _log.error("Mismatching R1/R2 files for %r", tmpl)
        return

    for nth, (file_1, file_2) in enumerate(zip(files_1, files_2), start=1):
        destination = args.output / f"{sample}_{library}_{lane}_{nth:03g}.txt"
        to_be_run = True

        if file_1 == file_2:
            _log.warning("Skipping SE data: %s > %s > %s", sample, library, lane)
            to_be_run = False
        elif destination.exists():
            _log.info("Skipping already done data: %s > %s > %s", sample, library, lane)
            to_be_run = False

        command: List[str] = [
            args.executable,
            "--file1",
            quote(file_1),
            "--file2",
            quote(file_2),
            "--adapter1",
            quote(args.adapter1),
            "--adapter2",
            quote(args.adapter2),
            "--identify-adapters",
            "--threads",
            str(args.threads),
        ]

        if args.executable == "adapterremoval3":
            command += ["--log-progress", "log"]

        yield to_be_run, Command(
            dry_run=args.dry_run,
            key=(sample, library, lane, str(nth)),
            destination=destination,
            args=tuple(command),
        )


async def worker(task: Command) -> bool:
    if task.dry_run:
        _log.info("Dry running %s", task.destination)
        _log.info("  Command = %s", " ".join(map(quote, task.args)))
        return True

    _log.info("Starting to create %s", task.destination)
    _log.debug("  Command = %s", " ".join(map(quote, task.args)))

    proc = await asyncio.create_subprocess_exec(
        *task.args,
        stdin=asyncio.subprocess.DEVNULL,
        stdout=asyncio.subprocess.PIPE,
        # stderr=asyncio.subprocess.PIPE,
        close_fds=True,
    )

    stdout, _ = await proc.communicate()
    if proc.returncode:
        _log.error("Error while creating %s:", task.destination)
        return False

    Path(task.destination).write_bytes(stdout)
    _log.info("Finished creating %s", task.destination)

    return True


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("yaml", nargs="+", type=Path)
    parser.add_argument("--instances", type=int, default=4)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output", default=Path("output"), type=Path)
    parser.add_argument("--adapter1", default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
    parser.add_argument("--adapter2", default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
    parser.add_argument("--executable")
    parser.add_argument("--dry-run", action="store_true")

    return parser.parse_args(argv)


async def main(argv: Sequence[str]):
    args = parse_args(argv)
    coloredlogs.install(
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    args.executable = determine_ar_executable(args.executable)
    if args.executable is None:
        return 1

    _log.info("Using AdapterRemoval executable %s", args.executable)

    tasks: List[Command] = []
    output_files: Set[str] = set()
    for filepath in args.yaml:
        _log.info("Reading YAML file %s", quote(filepath))

        for sample, library, lane, tmpl in sorted(safe_load(filepath), reverse=True):
            commands = build_command(
                args=args,
                sample=sample,
                library=library,
                lane=lane,
                tmpl=tmpl,
            )

            n_tasks_generated = 0
            for to_be_run, task in commands:
                n_tasks_generated += 1
                if not to_be_run:
                    continue

                assert task is not None
                if task.destination in output_files:
                    _log.error("File is clobbered: %s", task.destination)
                    return 1

                output_files.add(task.destination)
                tasks.append(task)

            if not n_tasks_generated:
                _log.error("Aborted")
                return 1

    if not args.dry_run:
        os.makedirs(args.output, exist_ok=True)

    n_tasks = len(tasks)
    n_done = 0
    n_errors = 0
    pending: Set[asyncio.Task[bool]] = set()
    interrupted = False

    loop = asyncio.get_event_loop()

    def _on_signal(sig):
        nonlocal interrupted

        interrupted = True
        for task in pending:
            task.cancel()

        loop.remove_signal_handler(sig)

    for sig in [signal.SIGINT, signal.SIGTERM]:
        loop.add_signal_handler(sig, _on_signal, sig)

    _log.info("Running %i tasks", len(tasks))
    while (tasks or pending) and not interrupted:
        while tasks and len(pending) < args.instances:
            pending.add(asyncio.create_task(worker(tasks.pop())))

        done, pending = await asyncio.wait(
            pending,
            return_when=asyncio.FIRST_COMPLETED,
        )

        n_done += len(done)
        for task in done:
            if task.cancelled():
                _log.warning("task was cancelled")
                n_errors += 1
            elif task.exception():
                _log.error("Error while running task: %s", task.exception())
                n_errors += 1
            elif not task.result():
                n_errors += 1

        _log.info(f"Status: %i of %i, %i errors", n_done, n_tasks, n_errors)

    if args.dry_run:
        _log.info("Dry run complete")

    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main(sys.argv[1:])))
