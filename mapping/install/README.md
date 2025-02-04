# Required software

## Installation using conda

```bash
# Create and activate conda environment with tools used for mapping and filtering
conda env create -f "conda.yaml" --prefix "./conda"
# Activates the conda environment; this script can be run from any location
source activate.sh

# Install paleomix in development mode; this means that changes made to the
# ./paleomix folder are reflected when running paleomix via conda
cd paleomix
python3 "./setup.py" develop
```

# Running tools

To run e.g. paleomix through the conda environment, remember to first activate the environment (if not already active):

```bash
source activate.sh
paleomix
python3 ./install/scripts/finalize_bam.py
```
