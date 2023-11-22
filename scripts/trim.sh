#!/bin/bash

set -e

# Load HCC modules
# . /util/opt/lmod/lmod/init/profile
# export -f module
# module use /util/opt/hcc-modules/Common/
# module load anaconda
source /opt/conda/etc/profile.d/conda.sh
conda activate /opt/ProkEvo_dir/prokevo

# trimmomatic "$@"
trimmomatic PE -threads 1 $1 $2 $3 $4 $5 $6 \
    HEADCROP:15 CROP:200 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:50

conda deactivate
