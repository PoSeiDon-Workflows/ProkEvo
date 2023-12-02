#!/usr/bin/env python3

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
    Container,
)


dax = Workflow("pipeline")
tc = TransformationCatalog()
rc = ReplicaCatalog()
sc = SiteCatalog()
base_dir = os.getcwd()

sra_file = Path(sys.argv[1]).resolve()
run_dir = sra_file.parent.name
container = Container(
    "prokevo",
    Container.DOCKER,
    # image="docker://pegasus/prokevo:latest",
    image="file://" + str(Path("data/prokevo.tar").resolve()),
    image_site="local",
)
tc.add_containers(container)

trim_run = []
fastqc_run = []
spades_run = []
quast_run = []
filtering_run = []
prokka_run = []
plasmidfinder_run = []
sistr_run = []
forward_file = []
reverse_file = []
list_of_fastqc_files = []
list_of_gff_files = []
list_of_contig_files = []
list_of_sistr_files = []

output_filtering_contigs = []

# Open input list and count files
input_file = open(sra_file)
lines = input_file.readlines()
length = len(lines)

rc.add_regex_replica("local", ".*", (Path("data") / "[0]").resolve())

for i in range(0, length):
    srr_id = lines[i].strip()

    forward_file.append(File(f"{srr_id}_1.fastq"))
    reverse_file.append(File(f"{srr_id}_2.fastq"))

    # add job for Trimmomatic
    # rc.add_replica(
    #     "local", f"{srr_id}_1.fastq", (Path("data") / f"{srr_id}_1.fastq").resolve()
    # )
    # rc.add_replica(
    #     "local", f"{srr_id}_2.fastq", (Path("data") / f"{srr_id}_2.fastq").resolve()
    # )
    trim_run.append(Job("ex_trim_run"))
    trim_run[i].add_args(
        forward_file[i].lfn,
        reverse_file[i].lfn,
        f"{srr_id}_pair_1_trimmed.fastq",
        f"{srr_id}_unpair_1_trimmed.fastq",
        f"{srr_id}_pair_2_trimmed.fastq",
        f"{srr_id}_unpair_2_trimmed.fastq",
    )
    trim_run[i].add_inputs(forward_file[i])
    trim_run[i].add_inputs(reverse_file[i])
    trim_run[i].add_outputs(f"{srr_id}_pair_1_trimmed.fastq", stage_out=False)
    trim_run[i].add_outputs(f"{srr_id}_unpair_1_trimmed.fastq", stage_out=False)
    trim_run[i].add_outputs(f"{srr_id}_pair_2_trimmed.fastq", stage_out=False)
    trim_run[i].add_outputs(f"{srr_id}_unpair_2_trimmed.fastq", stage_out=False)
    # trim_run[i].add_profiles(Profile("pegasus", "label", srr_id))
    dax.add_jobs(trim_run[i])

    # add job for FastQC
    fastqc_run.append(Job("ex_fastqc_run"))
    fastqc_run[i].add_args(
        f"{srr_id}_pair_1_trimmed.fastq", f"{srr_id}_pair_2_trimmed.fastq"
    )
    fastqc_run[i].add_inputs(f"{srr_id}_pair_1_trimmed.fastq")
    fastqc_run[i].add_inputs(f"{srr_id}_pair_2_trimmed.fastq")
    fastqc_run[i].add_outputs(
        f"{srr_id}_pair_1_trimmed_fastq_summary.txt",
        stage_out=False,
    )
    fastqc_run[i].add_outputs(
        f"{srr_id}_pair_2_trimmed_fastq_summary.txt",
        stage_out=False,
    )
    dax.add_jobs(fastqc_run[i])
    # add files
    f1 = File(f"{srr_id}_pair_1_trimmed_fastq_summary.txt")
    f2 = File(f"{srr_id}_pair_2_trimmed_fastq_summary.txt")
    list_of_fastqc_files.append(f1)
    list_of_fastqc_files.append(f2)

    # add job for Spades
    spades_run.append(Job("ex_spades_run"))
    spades_run[i].add_args(
        f"{srr_id}_pair_1_trimmed.fastq",
        f"{srr_id}_pair_2_trimmed.fastq",
        f"{srr_id}_spades_output",
    )
    spades_run[i].add_inputs(f"{srr_id}_pair_1_trimmed.fastq")
    spades_run[i].add_inputs(f"{srr_id}_pair_2_trimmed.fastq")
    spades_run[i].add_outputs(f"{srr_id}_spades_output/contigs.fasta")
    # Rajiv spades_run[i].add_profiles(Namespace.PEGASUS, "runtime", "3600")
    # spades_run[i].add_profiles(Profile("pegasus", "label", srr_id))
    dax.add_jobs(spades_run[i])

    # add job for Quast
    quast_run.append(Job("ex_quast_run"))
    quast_run[i].add_args(
        f"{srr_id}_quast_output", f"{srr_id}_spades_output/contigs.fasta"
    )
    quast_run[i].add_inputs(f"{srr_id}_spades_output/contigs.fasta")
    quast_run[i].add_outputs(f"{srr_id}_quast_output/transposed_report.tsv")
    dax.add_jobs(quast_run[i])

    # add job for Filterig contigs
    filtering_run.append(Job("ex_filtering_run"))
    filtering_run[i].add_args(
        f"{srr_id}_quast_output/transposed_report.tsv",
        f"{srr_id}_spades_output/contigs.fasta",
        run_dir,
        srr_id,
    )
    filtering_run[i].add_inputs(f"{srr_id}_quast_output/transposed_report.tsv")
    filtering_run[i].add_inputs(f"{srr_id}_spades_output/contigs.fasta")
    filtering_run[i].add_outputs(f"{srr_id}_contigs.fasta")
    output_filtering_contigs.append(f"{srr_id}_contigs.fasta")
    list_of_contig_files.append(f"{srr_id}_contigs.fasta")
    dax.add_jobs(filtering_run[i])

    # ------------
    # Sub Workflow
    # ------------

    prokka_run.append(Job("ex_prokka_run"))
    prokka_run[i].add_args(
        srr_id, f"{srr_id}_prokka_output", output_filtering_contigs[i]
    )
    prokka_run[i].add_inputs(output_filtering_contigs[i])
    prokka_run[i].add_outputs(
        f"{srr_id}_prokka_output/{srr_id}.gff",
        stage_out=True,
    )
    prokka_run[i].add_outputs(f"{srr_id}_prokka_output.tar.gz", stage_out=True)
    # prokka_run[i].add_profiles(Namespace.PEGASUS, "label", srr_id)
    # Rajiv prokka_run[i].add_profiles(Namespace.PEGASUS, "runtime", "14400")
    dax.add_jobs(prokka_run[i])
    # add files
    f = File(f"{srr_id}_prokka_output/{srr_id}.gff")
    list_of_gff_files.append(f)

    # add job for plasmidfinder
    plasmidfinder_run.append(Job("ex_plasmidfinder_run"))
    plasmidfinder_run[i].add_args(
        output_filtering_contigs[i], f"{srr_id}_plasmidfinder_output"
    )
    plasmidfinder_run[i].add_inputs("9cdf35065947.tar.gz")
    # rc.add_replica(
    #     "local", "9cdf35065947.tar.gz", (Path("data") / "9cdf35065947.tar.gz").resolve()
    # )
    plasmidfinder_run[i].add_inputs(output_filtering_contigs[i])
    plasmidfinder_run[i].add_outputs(
        f"{srr_id}_plasmidfinder_output.tar.gz", stage_out=True
    )
    dax.add_jobs(plasmidfinder_run[i])

    # add job for sistr
    sistr_run.append(Job("ex_sistr_run"))
    sistr_run[i].add_args(
        f"{srr_id}_allele_results.json",
        f"{srr_id}_novel_alleles.fasta",
        f"{srr_id}_cgmlst_profiles.csv",
        f"{srr_id}_sistr_output.csv",
        output_filtering_contigs[i],
    )
    sistr_run[i].add_inputs(output_filtering_contigs[i])
    sistr_run[i].add_outputs(f"{srr_id}_sistr_output.csv", stage_out=False)
    dax.add_jobs(sistr_run[i])
    list_of_sistr_files.append(f"{srr_id}_sistr_output.csv")


# add job for cat FastQC Fail
ex_cat = Transformation(
    namespace="dax",
    name="cat",
    version="4.0",
    site="condorpool",
    pfn="/bin/cat",
    os_type=OS.LINUX,
    arch=Arch.AARCH64,
    is_stageable=False,
    container=container,
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

# ------------
# Sub Workflow
# ------------

# add job for mlst
mlst_run = Job("ex_mlst_run")
mlst_run.add_args(*list_of_contig_files)
o = File("mlst_output.csv")
mlst_run.set_stdout(o, stage_out=True)
# Rajiv mlst_run.add_profiles(Namespace.PEGASUS, "runtime", "108000")
# mlst_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id))
dax.add_jobs(mlst_run)


# add job for abricate vfdb
abricate_vfdb_run = Job("ex_abricate_run")
abricate_vfdb_run.add_args("vfdb", *list_of_contig_files)

o = File("sabricate_vfdb_output.csv")
abricate_vfdb_run.set_stdout(o, stage_out=True)
dax.add_jobs(abricate_vfdb_run)

# add job for abricate argannot
abricate_argannot_run = Job("ex_abricate_run")
abricate_argannot_run.add_args("argannot", *list_of_contig_files)
o = File("sabricate_argannot_output.csv")
abricate_argannot_run.set_stdout(o, stage_out=True)
dax.add_jobs(abricate_argannot_run)

# add job for abricate card
abricate_card_run = Job("ex_abricate_run")
abricate_card_run.add_args("card", *list_of_contig_files)
o = File("sabricate_card_output.csv")
abricate_card_run.set_stdout(o, stage_out=True)
dax.add_jobs(abricate_card_run)

# add job for abricate ncbi
abricate_ncbi_run = Job("ex_abricate_run")
abricate_ncbi_run.add_args("ncbi", *list_of_contig_files)
o = File("sabricate_ncbi_output.csv")
abricate_ncbi_run.set_stdout(o, stage_out=True)
dax.add_jobs(abricate_ncbi_run)

# add job for abricate plasmidfinder
abricate_plasmidfinder_run = Job("ex_abricate_run")
abricate_plasmidfinder_run.add_args("plasmidfinder", *list_of_contig_files)
o = File("sabricate_plasmidfinder_output.csv")
abricate_plasmidfinder_run.set_stdout(o, stage_out=True)
dax.add_jobs(abricate_plasmidfinder_run)

# add job for abricate resfinder
abricate_resfinder_run = Job("ex_abricate_run")
abricate_resfinder_run.add_args("resfinder", *list_of_contig_files)
o = File("sabricate_resfinder_output.csv")
abricate_resfinder_run.set_stdout(o, stage_out=True)
dax.add_jobs(abricate_resfinder_run)

for l in list_of_contig_files:
    mlst_run.add_inputs(l)
    abricate_vfdb_run.add_inputs(l)
    abricate_argannot_run.add_inputs(l)
    abricate_card_run.add_inputs(l)
    abricate_ncbi_run.add_inputs(l)
    abricate_plasmidfinder_run.add_inputs(l)
    abricate_resfinder_run.add_inputs(l)

# add job for Roary
roary_run = Job("ex_roary_run")
roary_run.add_args("roary_output", *list_of_gff_files)
for l in list_of_gff_files:
    roary_run.add_inputs(l)
roary_run.add_outputs("roary_output/core_gene_alignment.aln", stage_out=True)
roary_run.add_outputs("roary_output.tar.gz", stage_out=True)
# Rajiv roary_run.add_profiles(Namespace.PEGASUS, "runtime", "604800")
# Rajiv roary_run.add_profiles(Namespace.CONDOR, "request_memory", "970000")
# Rajiv roary_run.add_profiles(Namespace.CONDOR, "memory", "970000")
roary_run.add_profiles(Namespace.PEGASUS, "memory", "2048")
dax.add_jobs(roary_run)

# add job for baps_run
# R script, wrapper
fastbaps_run = Job("ex_fastbaps_run")
fastbaps_output = File("fastbaps_baps.csv")
fastbaps_run.add_args("roary_output/core_gene_alignment.aln", fastbaps_output)
fastbaps_run.add_inputs("roary_output/core_gene_alignment.aln")
fastbaps_run.add_outputs(fastbaps_output, stage_out=True)
# Rajiv fastbaps_run.add_profiles(Namespace.CONDOR, "request_memory", "30000")
# fastbaps_run.add_profiles(Namespace.GLOBUS, "maxmemory", "30000")
fastbaps_run.add_profiles(Namespace.PEGASUS, "memory", "2048")
# Remove as R code segfaults
# dax.add_jobs(fastbaps_run)

output_sistr_cat = File("sistr_all.csv")
cat = Job(ex_cat, namespace="dax")
cat.add_args(*list_of_sistr_files)
for l in list_of_sistr_files:
    cat.add_inputs(l)
cat.set_stdout(output_sistr_cat, stage_out=True, register_replica=False)
dax.add_jobs(cat)

# add job for sistr output filtering
output_sistr_merge_cat = File("sistr_all_merge.csv")
merge_sistr_run = Job("ex_merge_sistr")
merge_sistr_run.add_args(output_sistr_cat, output_sistr_merge_cat)
merge_sistr_run.add_inputs(output_sistr_cat)
merge_sistr_run.add_outputs(output_sistr_merge_cat, stage_out=True)
dax.add_jobs(merge_sistr_run)


# Write the DAX to stdout
dax.add_transformation_catalog(tc)
# dax.add_site_catalog(sc)
# dax.add_replica_catalog(rc)
rc.write("replicas.yml")
dax.write(sys.stdout)
