.PHONEY = adapters

# Needed to support `source`, since default shell (sh) doesn't and for advanced globs
SHELL := /bin/bash -O extglob -c

# Prefix that ensures that commands run with activated modules
ACTIVATED := source install/activate.sh && nice -19 ionice -c 3

# Override by running `make BATCH=batch_1`
BATCH := batch_1
# Find YAML files for this batch. Requires extglob to match ${BATCH}.yaml and ${BATCH}.blah.yaml
YAMLFILES := ${BATCH}*(.*).yaml

help:
	@echo 'Usage: make <adapters/mapping/preliminary/filtering> BATCH=$${batch}'
	@echo "  e.g. make mapping BATCH=batch_1"

adapters:
	# Support the case where output is a symlink to a removed folder
	mkdir -p $(shell readlink -f ${BATCH}.adapters)
	${ACTIVATED} python3 install/scripts/identify_adapters.py $(YAMLFILES) --output ${BATCH}.adapters
	${ACTIVATED} python3 install/scripts/classify_adapters.py ${BATCH}.adapters | tee ${BATCH}.adapters.tsv

mapping:
	# Support the case where output is a symlink to a removed folder
	mkdir -p $(shell readlink -f ${BATCH}.raw_bams)
	mkdir -p ${BATCH}.raw_bams/logs
	${ACTIVATED} paleomix bam run \
		--destination ${BATCH}.raw_bams/ \
		--temp-root ${BATCH}.raw_bams/temp \
		--log-file ${BATCH}.raw_bams/logs/${HOSTNAME}_$(shell date +"%Y%m%d_%H%M%S").txt \
		$(YAMLFILES)

filtering:
	${ACTIVATED} snakemake --snakefile install/scripts/Snakefile.filter_bams --config batch=${BATCH} -j 32 -p

insert_sizes:
	${ACTIVATED} snakemake --snakefile install/scripts/Snakefile.insert_sizes --config batch=${BATCH} -j 32 -pn

coverage:
	${ACTIVATED} python3 install/scripts/coverage_statistics.py ${BATCH}/ > ${BATCH}.coverage.tsv

