This pipeline with generate a series of QC statistics from your sequencing data. The pipeline requires that the sequencing data was mapped using the PALEOMIX pipeline. It is also required that the data was mapped to a distant and close reference.

The pipeline comes with a script (`workflow/scripts/make_config.py`) to generate the configuration file used for input. It is recommended to use the conda setup to ensure that all software is available.

# Requirements

The pipeline requires that the `snakemake` is installed

```bash
python3 -m pip install snakemake --user
```
Note that snakemake requires python 3.7+.

# Running the pipeline

It is possible to run the pipeline using either Mamba / Conda or by having all the required software in your path.

# setup Mamba / Conda

Download and install mamba from https://github.com/conda-forge/miniforge#mambaforge

To ensure that conda/mamba is not loaded at every log in and to set the channels it is recommended to add the following to your `~/.condarc`:
```yaml
auto_activate_base: false
## https://conda-forge.org/docs/user/tipsandtricks.html#using-multiple-channels
channel_priority: strict
channels:
- conda-forge
- defaults
```

# Setup configuration file

It is recommend to divide the pipeline into three steps:

1. Check fastq data *prior* to mapping
  ```bash
  cores=4
  PALEOMIX_YAML="/path/to/paleomix.yaml"

  python3 workflow/scripts/make_config.py \
    --paleomix_yaml ${PALEOMIX_YAML} \
    --output_yaml prior.yaml

  python3 -m snakemake --use-conda \
    --configfile prior.yaml \
    --snakefile /path/to/workflow/Snakefile \
    -p -c ${cores} RUN_multiqc_pre
  ```
  This will generate a multiQC report of all the sequencing data. It might be that you want to exclude lanes or entire samples.

2. Check mapping statistics after mapping to e.g. exclude samples without data
  ```bash
  cores=4
  PALEOMIX_YAML="/path/to/paleomix.yaml"
  BAM_DIR_PATH="/path/to/bam"

  python3 workflow/scripts/make_config.py \
    --paleomix_yaml ${PALEOMIX_YAML} \
    --bams ${BAM_DIR_PATH}/*bam \
    --close_ref_name GrantsGazelle \
    --distant_ref_name NAME_IN_PALEOMIX \
    --perfect_bam_name NAME_IN_PALEOMIX  \
    --output_yaml mapstat.yaml

  python3 -m snakemake --use-conda \
    --configfile mapstat.yaml \
    --snakefile /path/to/workflow/Snakefile \
    -p -c ${cores} RUN_json
  ```
  Based on mapping statistics you might want to excluded some samples.

3. Run the remaining part of the pipeline
  ```bash
  cores=4  ## bump this to a fun number
  PALEOMIX_YAML="/path/to/paleomix.yaml"
  BAM_DIR_PATH="/path/to/bam"

  python3 workflow/scripts/make_config.py \
    --paleomix_yaml ${PALEOMIX_YAML} \
    --bams ${BAM_PATH}/*bam --output_yaml main.yaml \
    --close_ref_name GrantsGazelle \
    --distant_ref_name NAME_IN_PALEOMIX \
    --perfect_bam_name NAME_IN_PALEOMIX  \
    --exclude SAMPLE_NAME1 SAMPLE_NAME2 ## excludes these two samples from all downstream analyses

  python3 -m snakemake --use-conda \
    --configfile main.yaml \
    --snakefile /path/to/workflow/Snakefile \
    -p -c ${cores} RUN_all
  ```

## Running the software

### with conda

```bash
python3 -m snakemake --use-conda --configfile conf.yaml RUN_all
```

### without conda

```bash
python3 -m snakemake --configfile conf.yaml RUN_all
```


#### Dependencies required when not using Mamba / Conda

Following software should be in your $PATH

- python3
  - numpy
  - pandas
- angsd
- fastqc
- multiqc
- samtools
- bcftools
- plink
- plink2
- winsfs
- r
  - ggplot2
  - data.table
  - reshape2
  - RColorBrewer
  - ape
  - bit64
  - mclust
