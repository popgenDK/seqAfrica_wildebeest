# Depth filtering based on genotype calls

Needs as input a bcf file with all sites, including non-variable, and ideally no filtering. It extracts the per site combined depth across all samples,
plots the depth distribution and generates bed files with sites to keep (`good.bed`) and remove (`bad.bed`) based on minimum and maximum depth thersholds.

The minimum and maximum depth thresholds need to be set in the config. To set them you can do a frist run to just extract and plot the depth distribution using the 
rule `check_distribution` as target rule:

```snakemake --snakefile filter_depth.snakefile --configfile config.yaml -j 20 check_distribtuion```

This will produce `depth_pos.txt` file with global depth for each positon, and ```depth_distribution.png``` with a histogram of per site depth distribution.
From here you can get the information to choose optimal depth thersholds, there is not a single way to set them, but they should keep the main central peak of depths
and exclude outliers (often you will see peaks in the very low depth and very high depth areas, that is what you want to get rid off).

After having set the fitlers, you can do a normal run with will apply it and produce the bed files:

```snakemake --snakefile filter_depth.snakefile --configfile config.yaml -j 20```
