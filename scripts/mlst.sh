#!/bin/bash

module load anaconda
conda activate prokevo

# mlst "$@"
# Change scheme based on organism
mlst --legacy --scheme senterica --csv "$@"

conda deactivate
