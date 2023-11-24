# Mapping and filtering of wildebeests

The following describes the commands used to identify/verify adapters, map, and filter BAMs for one or more batches of samples.

----------------------------------------------------------------------------------------

## Software used

 * paleomix, 'pub/2022/africa' branch
 * python v3.9.9
 * adapterremoval v2.3.2
 * bwa v0.7.17
 * bwa-mem2 v2.2.1
 * samtools v1.11

### Python modules

 * coloredlogs
 * pysam
 * ruamel.yaml

See `install/scripts.requirements.txt` for exact module versions used.

----------------------------------------------------------------------------------------

## 1. Construction of `paleomix bam` YAML files

The YAML configuration files for the `paleomix bam` command are required not only for running the mapping pipeline, but are also used by the adapter identification step to locate PE FASTQ files. 

Initial YAML files were generated as follows:

```bash
    paleomix bam makefile > project.batch_1.yaml
    python3 scripts/install/tsv_files_to_yaml.py project.batch_1.tsv >> project.batch_1.yaml

    paleomix bam makefile > project.batch_2.yaml
    python3 scripts/install/tsv_files_to_yaml.py project.batch_2.tsv >> project.batch_2.yaml
```

The YAML files were tweaked to minimize read filtering and trimming of low quality bases. 

----------------------------------------------------------------------------------------

## 2. Identification of adapters

To ensure that the correct adapter sequences were trimmed from all samples, `AdapterRemoval --identify-adapters` was run on all PE reads and the output was compared with the recommended [BGI](https://en.mgi-tech.com/Download/download_file/id/71) or [Illumina](https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html) adapters sequences for trimming:

    Illumina forward: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    Illumina reverse: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

    BGI forward: AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
    BGI reverse: AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG

A makefile is provided for running the scripts for each batch:

```bash
    make adapters BATCH=batch_1
    make adapters BATCH=batch_2
```

Individual output files from `AdapterRemoval --identify-adapters` are written to `${BATCH}.adapters/` and the resulting `${BATCH}.adapters.tsv` files contains attemped automatic classifications of adapter sequences, reporting the best match with either BGI or Illumina sequences.

The following adapters were identified and the AdapterRemoval settings `--adapter1` and `--adapter2` were updated in the YAML files using those sequences:

 * Batch 1: Standard Illumina adapters.

----------------------------------------------------------------------------------------

# 3. Read mapping

Read mapping was performed using a development version of [PALEOMIX](https://github.com/mikkelschubert/paleomix). Minimal filtering is performed during this run (see above), with the final BAM file containing all input reads (some pairs of which may be merged into a single sequence), except for empty reads.

A makefile is provided for running the pipeline for each batch:

```bash
    make mapping BATCH=batch_1
```

The output BAMs, temporary files, and logs are written to `${BATCH}.raw_bams/`.

---------------------------------------------------------------------------------------

# 4. Insert size quantiles

Run the SnakeMake file to collect insert size quantiles:

```bash
    make insert_sizes BATCH=batch_1
    make insert_sizes BATCH=batch_2
```

The resulting tables are written to `${BATCH}.insert_sizes.tsv`.

---------------------------------------------------------------------------------------

# 5. Filtering and SAMTools statistics

Run the SnakeMake file to filter the BAMs and to collect SAMTools statistics:

```bash
    make filtering BATCH=batch_1
    make filtering BATCH=batch_2
```

This writes the filtered BAMs and statistics files to `${BATCH}/`.

---------------------------------------------------------------------------------------

# 6. Coverage statistics

Coverage of filtered and retained reads were gathered following filtering:


```bash
    make coverage BATCH=batch_1
    make coverage BATCH=batch_2
```

The resulting tables are written to `${BATCH}.coverage.tsv`.
