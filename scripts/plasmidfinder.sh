#!/bin/bash

# Load HCC modules
# . /util/opt/lmod/lmod/init/profile
# export -f module
# module use /util/opt/hcc-modules/Common/
# module load anaconda
source /opt/conda/etc/profile.d/conda.sh
conda activate /opt/ProkEvo_dir/prokevo

# Download Plasmidfinder db
mkdir ./dbs
cd ./dbs
tar -xvf 9cdf35065947.tar.gz --strip-components 1 && rm *.tar.gz && python INSTALL.py
cd ..

export PLASMID_DB=../dbs/

mkdir $2
plasmidfinder.py -p $PLASMID_DB -i $1 -o $2
tar -czvf $2.tar.gz $2

conda deactivate
