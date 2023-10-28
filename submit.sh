#!/bin/bash

set -e

export PYTHONPATH=`pegasus-config --python`

TOPDIR=`pwd`

<<DO_NOT_COMMENT_OUT
You can use raw Illumina sequences stored locally instead of downloading them from NCBI.
To do this, you will need to add information about the input files in rc.yml.
You can add these files in a loop as shown below, or any way you prefer.
You should also comment out the dependency in root-dax.py to skip the download.
DO_NOT_COMMENT_OUT

<<COMM
while read line
do
echo ''${line}'_1.fastq file:///absolute_path_to_fastq_files/'${line}'_1.fastq site="local"' >> rc.yml
echo ''${line}'_2.fastq file:///absolute_path_to_fastq_files/'${line}'_2.fastq site="local"' >> rc.yml
done < sra_ids.txt
COMM

# Clean old directories and files
rm -rf data_tmp
rm -rf scratch
rm -rf outputs
rm -rf prokevo
rm -rf $USER
rm -rf root-pipeline.dax
rm -rf sites.xml
rm -rf rc.yml
cp rc.yml.org rc.yml

# Set working path to current directory
sed -i "s|ProkEvo_dir|$PWD|g" tc.yml
sed -i "s|ProkEvo_dir|$PWD|g" rc.yml
for i in scripts/*.sh
do
sed -i "s|ProkEvo_dir|$PWD|g" $i
done

export RUN_DIR=$TOPDIR/data_tmp
mkdir -p $RUN_DIR
./root-dax.py $RUN_DIR > root-pipeline.dax

# create the site catalog
# this section contains the information about the running site
cat > sites.yml <<EOF
---
pegasus: "5.0"
sites:
 -
  name: "local-hcc"
  arch: "x86_64"
  os.type: "linux"
  directories:
   -
    type: "localStorage"
    path: "${PWD}/outputs"
    sharedFileSystem: false
    fileServers:
     -
      operation: "all"
      url: "file://${PWD}/outputs"
   -
    type: "sharedScratch"
    path: "${PWD}/scratch"
    sharedFileSystem: false
    fileServers:
     -
      operation: "all"
      url: "file://${PWD}/scratch"
  profiles:
    env:
      PEGASUS_HOME: "/util/opt/pegasus-wms/5.0/"
    condor:
      grid_resource: "batch slurm"
      request_memory: "ifthenelse(isundefined(DAGNodeRetry) || DAGNodeRetry == 0,\
        \ 2000, 120000)"
    pegasus:
      auxillary.local: "true"
      queue: "batch"
      style: "glite"
EOF


# plan and submit the root workflow
pegasus-plan --conf pegasusrc --sites local-hcc --output-site local-hcc --dir ${PWD} --submit root-pipeline.dax # --cluster label

# to resume/restart fixed workflow
# pegasus-run <run_directory>
