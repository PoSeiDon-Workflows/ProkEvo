#!/usr/bin/env python

"""
NOTE:
- Comment out NCBI Download (sra_run) dependency line at the end if local Illumina sequences are used.
"""

import sys
import os
from pathlib import Path

from Pegasus.api import (
    Workflow,
    SubWorkflow,
    TransformationCatalog,
    ReplicaCatalog,
    SiteCatalog,
    File,
    Job,
    Transformation,
    Namespace,
    OS,
    Arch,
)


dax = Workflow("pipeline")
tc = TransformationCatalog()
rc = ReplicaCatalog()
sc = SiteCatalog()
base_dir = os.getcwd()

run_dir = sys.argv[1]

sra_run = []
trim_run = []
fastqc_run = []
spades_run = []
quast_run = []
filtering_run = []
prokka_run = []
plasmidfinder_run = []
forward_file = []
reverse_file = []
list_of_fastqc_files = []
list_of_gff_files = []
list_of_contig_files = []
list_of_filtererd_sra_ids = []


# Open input list and count files
input_file = open("sra_ids.txt")
rc.add_replica("local", "sra_ids.txt", Path("sra_ids.txt").resolve())
lines = input_file.readlines()
length = len(lines)

# Add file executable and job for sub-pipeline
c = File("sub-pipeline.yml")

# Add file for conda environment
conda_file = File("prokevo.yml")
rc.add_replica("local", conda_file.lfn, Path(conda_file.lfn).resolve())

# Add a job to analyze the output of split and generate a sub dax with correct number of parallelism based on output of previous job
generate = Job("ex_generate")
generate.add_args(run_dir)
generate.set_stdout(c)
# Rajiv generate.add_profiles(Namespace.ENV, key="PYTHONPATH", value=os.environ["PYTHONPATH"])
generate.add_profiles(Namespace.SELECTOR, key="execution.site", value="local")
dax.add_jobs(generate)

# Add a subdax job of type DAX that takes the runtime generated sub dax file in the previous step and runs the computation.
sub_dax = SubWorkflow(c)
sub_dax.add_planner_args(
    basename="sub-pipeline", sites=["local-hcc"], output_sites=["local-hcc"]
)
dax.add_jobs(sub_dax)


# Start analysis
# add job to create conda environment
conda_run = Job("ex_conda_run")
conda_run.add_args(conda_file)
conda_run.add_inputs(conda_file)
dax.add_jobs(conda_run)

for i in range(0, length):
    srr_id = lines[i].strip()

    forward_file.append(File(str(srr_id) + "_1.fastq"))
    reverse_file.append(File(str(srr_id) + "_2.fastq"))
    # Rajiv dax.add_inputs(forward_file[i])
    # Rajiv dax.add_inputs(reverse_file[i])

    # add job for downloading data from NCBI
    sra_run.append(Job("ex_sra_run"))
    sra_run[i].add_args(str(srr_id))
    sra_run[i].add_outputs(forward_file[i], stage_out=False)
    sra_run[i].add_outputs(reverse_file[i], stage_out=False)
    # add profile for download limit
    # Profile(PROPERTY_KEY[0], PROFILE KEY, PROPERTY_KEY[1])
    sra_run[i].add_profiles(Namespace.DAGMAN, "CATEGORY", "sradownload")
    # sra_run[i].add_profiles(Profile("pegasus", "label", str(srr_id)))
    dax.add_jobs(sra_run[i])

    # add job for Trimmomatic
    trim_run.append(Job("ex_trim_run"))
    trim_run[i].add_args(
        forward_file[i].lfn,
        reverse_file[i].lfn,
        str(srr_id) + "_pair_1_trimmed.fastq",
        str(srr_id) + "_unpair_1_trimmed.fastq",
        str(srr_id) + "_pair_2_trimmed.fastq",
        str(srr_id) + "_unpair_2_trimmed.fastq",
    )
    trim_run[i].add_inputs(forward_file[i])
    trim_run[i].add_inputs(reverse_file[i])
    trim_run[i].add_outputs(str(srr_id) + "_pair_1_trimmed.fastq", stage_out=False)
    trim_run[i].add_outputs(str(srr_id) + "_unpair_1_trimmed.fastq", stage_out=False)
    trim_run[i].add_outputs(str(srr_id) + "_pair_2_trimmed.fastq", stage_out=False)
    trim_run[i].add_outputs(str(srr_id) + "_unpair_2_trimmed.fastq", stage_out=False)
    # trim_run[i].add_profiles(Profile("pegasus", "label", str(srr_id)))
    dax.add_jobs(trim_run[i])

    # add job for FastQC
    fastqc_run.append(Job("ex_fastqc_run"))
    fastqc_run[i].add_args(
        str(srr_id) + "_pair_1_trimmed.fastq", str(srr_id) + "_pair_2_trimmed.fastq"
    )
    fastqc_run[i].add_inputs(str(srr_id) + "_pair_1_trimmed.fastq")
    fastqc_run[i].add_inputs(str(srr_id) + "_pair_2_trimmed.fastq")
    fastqc_run[i].add_outputs(
        str(srr_id) + "_pair_1_trimmed_fastqc/summary.txt",
        stage_out=False,
    )
    fastqc_run[i].add_outputs(
        str(srr_id) + "_pair_2_trimmed_fastqc/summary.txt",
        stage_out=False,
    )
    dax.add_jobs(fastqc_run[i])
    # add files
    f1 = File(str(srr_id) + "_pair_1_trimmed_fastqc/summary.txt")
    f2 = File(str(srr_id) + "_pair_2_trimmed_fastqc/summary.txt")
    list_of_fastqc_files.append(f1)
    list_of_fastqc_files.append(f2)

    # add job for Spades
    spades_run.append(Job("ex_spades_run"))
    spades_run[i].add_args(
        str(srr_id) + "_pair_1_trimmed.fastq",
        str(srr_id) + "_pair_2_trimmed.fastq",
        str(srr_id) + "_spades_output",
    )
    spades_run[i].add_inputs(str(srr_id) + "_pair_1_trimmed.fastq")
    spades_run[i].add_inputs(str(srr_id) + "_pair_2_trimmed.fastq")
    spades_run[i].add_outputs(str(srr_id) + "_spades_output/contigs.fasta")
    spades_run[i].add_profiles(Namespace.PEGASUS, "runtime", "3600")
    spades_run[i].add_profiles(Namespace.GLOBUS, "maxwalltime", "600")
    # spades_run[i].add_profiles(Profile("pegasus", "label", str(srr_id)))
    dax.add_jobs(spades_run[i])

    # add job for Quast
    quast_run.append(Job("ex_quast_run"))
    quast_run[i].add_args(
        str(srr_id) + "_quast_output", str(srr_id) + "_spades_output/contigs.fasta"
    )
    quast_run[i].add_inputs(str(srr_id) + "_spades_output/contigs.fasta")
    quast_run[i].add_outputs(str(srr_id) + "_quast_output/transposed_report.tsv")
    # quast_run[i].add_profiles(Profile("pegasus", "label", str(srr_id)))
    dax.add_jobs(quast_run[i])

    # add job for Filterig contigs
    filtering_run.append(Job("ex_filtering_run"))
    filtering_run[i].add_args(
        str(srr_id) + "_quast_output/transposed_report.tsv",
        str(srr_id) + "_spades_output/contigs.fasta",
        run_dir,
        str(srr_id),
    )
    filtering_run[i].add_inputs(str(srr_id) + "_quast_output/transposed_report.tsv")
    filtering_run[i].add_inputs(str(srr_id) + "_spades_output/contigs.fasta")
    # filtering_run[i].add_profiles(Profile("pegasus", "label", str(srr_id)))
    dax.add_jobs(filtering_run[i])

# add job for cat FastQC Fail
ex_cat = Transformation(
    namespace="dax",
    name="cat",
    version="4.0",
    site="local-hcc",
    pfn="/bin/cat",
    os_type=OS.LINUX,
    arch=Arch.X86_64,
    is_stageable=False,
)
tc.add_transformations(ex_cat)
output_fastqc_cat = File("fastqc_summary_all.txt")
cat = Job(ex_cat, namespace="dax")
cat.add_args(*list_of_fastqc_files)
for l in list_of_fastqc_files:
    cat.add_inputs(l)
cat.set_stdout(output_fastqc_cat, stage_out=True, register_replica=False)
dax.add_jobs(cat)

# add job for FastQC Fail
output_fastqc_fail = File("fastqc_summary_final.txt")
fastqc_fail_run = Job("ex_fastqc_fail_run")
fastqc_fail_run.add_args(output_fastqc_cat, output_fastqc_fail)
fastqc_fail_run.add_inputs(output_fastqc_cat)
fastqc_fail_run.add_outputs(output_fastqc_fail, stage_out=True)
dax.add_jobs(fastqc_fail_run)


for i in range(0, length):
    # Add control-flow dependencies
    # dax.add_dependency(trim_run[i], conda_run)
    # USE THE LINE ABOVE AND COMMENT OUT THE 2 LINES BELOW TO SKIP NCBI DOWNLOAD!!!
    dax.add_dependency(conda_run, children=[sra_run[i]])
    dax.add_dependency(sra_run[i], children=[trim_run[i]])
    dax.add_dependency(trim_run[i], children=[fastqc_run[i]])
    dax.add_dependency(fastqc_run[i], children=[cat])
    dax.add_dependency(trim_run[i], children=[spades_run[i]])
    dax.add_dependency(spades_run[i], children=[quast_run[i]])
    dax.add_dependency(quast_run[i], children=[filtering_run[i]])
    dax.add_dependency(filtering_run[i], children=[generate])
dax.add_dependency(cat, children=[fastqc_fail_run])
dax.add_dependency(generate, children=[sub_dax])


# Write the DAX to stdout
dax.add_transformation_catalog(tc)
dax.add_site_catalog(sc)
# dax.add_replica_catalog(rc)
dax.write(sys.stdout)
