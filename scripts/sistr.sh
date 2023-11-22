#!/bin/bash

# Load HCC modules
# . /util/opt/lmod/lmod/init/profile
# export -f module
# module use /util/opt/hcc-modules/Common/
# module load anaconda
source /opt/conda/etc/profile.d/conda.sh
conda activate /opt/ProkEvo_dir/prokevo

# sistr "$@"
sistr --qc -vv --alleles-output $1 --novel-alleles $2 --cgmlst-profiles $3 \
    -f csv -o $4 $5

conda deactivate
