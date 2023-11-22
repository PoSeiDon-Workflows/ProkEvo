#!/bin/bash

# Load HCC modules
# . /util/opt/lmod/lmod/init/profile
# export -f module
# module use /util/opt/hcc-modules/Common/
# module load anaconda
source /opt/conda/etc/profile.d/conda.sh
conda activate /opt/ProkEvo_dir/prokevo

# fastqc "$@"
fastqc $1 $2 --extract

NAME=$(echo "$1" | sed -e 's/\.fastq/_fastq/')
mv ${NAME}c/summary.txt ${NAME}_summary.txt

NAME=$(echo "$2" | sed -e 's/\.fastq/_fastq/')
mv ${NAME}c/summary.txt ${NAME}_summary.txt

conda deactivate
