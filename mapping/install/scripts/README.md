# Scripts

* `tsv_files_to_yaml.py`
  Generates a sample table for PALEOMIX, using a TSV file containing the columns `Sample`, `File1`, and `File2`.
* `finalize_bam.py`  
  Script for filtering BAM files and/or collecting statistics. Writing of filtered BAMs can be disabled by omitting the corresponding output options, in order to do preliminary collection of statistics. Depending on filters used, name-order sorting may be required for the input BAM.

# Misc statistics

* `insert_size_quantiles.py`
* `insert_size_sanity_check.py`

# Adapter identification/verification

* `identify_adapters.py`
  Script for running `AdapterRemoval --identify-adapters` on a paleomix YAML file (BAM pipeline). Prints the shell commands used for running AdapterRemoval, which can be piped to parallel.
* `classify_adapters.py`
  Parses the output of `AdapterRemoval --identify-adapters` and tries to guess the adapter sequence by comparison to the recommended BGI and Illumina sequences. Identify is provided as the number of matching bases from the start of the sequence, allowing at most one mismatch.

## Plotting filtering metrics

* `metrics_collect.py`  
  Simple script demonstrating how structured JSON statistics can be converted to TSV for use with `R`, etc.
* `metrics_plot.r`  
  R-script for plotting insert sizes, match percentages, etc. from TSV files generated using `metrics_collect.py`

## Snakemake files

* `Snakefile.filter_bams`  
  Filters BAMs in `input/` and writing the resulting good/bad BAMs to `output` along with `samtools stats` and `samtools idxstats` results.
* `Snakefile.insert_sizes`
  Snakefile for inferring insert size distributions with `insert_size_quantiles.py` on a set of BAM files.

## Misc files

* `setup.cfg`
  Custom settings for [flake8](https://flake8.pycqa.org/) linting for use with [black](https://github.com/psf/black).
