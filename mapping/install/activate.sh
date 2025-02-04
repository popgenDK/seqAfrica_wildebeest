#!/bin/bash

# README:
#   To activate, run (if using bash)
#     source install/activate.sh
#   The `return` statements ensures that `make` aborts if any of these commands fail

# Location of (my) mapping pipeline
PIPELINE_PATH="/path/to/mapping"
# Location of custom modules
MODULES_PATH="/path/to/modulefiles/"

if [ ! -O "${PIPELINE_PATH}" ]; then
	echo "ERROR: You do not own the folder '${PIPELINE_PATH}'" >&2
	echo "       Please update your 'activate.sh' script to point to the location of *your* pipeline" >&2
	return 1
fi

module load python/3.9.9 || return

if [ ! -d "${MODULES_PATH}" ]; then
	echo "ERROR: Modules folder '${MODULES_PATH}' does not exist or you do not have access" >&2
	echo "       Please update your 'activate.sh' script to point to the location the requires modules" >&2
	return 1
fi

module use --prepend "${MODULES_PATH}" || return

# Used for adapter identification only
# module load adapterremoval/3.0.0-alpha1 || return

# Used by mapping pipeline
module load adapterremoval/2.3.2 || return
module load bwa/0.7.17 || return
module load bwa-mem2/2.2.1 || return
module load samtools/1.11 || return

# Used for post-filtering QC
module load mamba/1.3.1 || return
module load winsfs/0.7.0-ff3104f || return
