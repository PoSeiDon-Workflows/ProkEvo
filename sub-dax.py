#!/usr/bin/python

"""
USAGE: # ./sub-dax.py $RUN_DIR > sub-pipeline.dax
"""

"""
NOTE:
- Comment out SISTR dependency line at the end if non-Salmonella organism is used.
"""

import sys
import os

from Pegasus.api import (
    Workflow,
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

prokka_run = []
plasmidfinder_run = []
sistr_run = []
list_of_gff_files = []
list_of_contig_files = []
list_of_filtererd_sra_ids = []
list_of_sistr_files = []


# Get list of contigs after filtering
for file_name in os.listdir(run_dir):
    rc.add_replica("local-hcc", file_name, f"file://{run_dir}/{file_name}")
    list_of_contig_files.append(File(file_name))


i = 0
for output_filtering_contigs in list_of_contig_files:
    srr_id = output_filtering_contigs.lfn.split("_")[0]

    # add job for Prokka
    prokka_run.append(Job("ex_prokka_run"))
    prokka_run[i].add_args(
        srr_id, str(srr_id) + "_prokka_output", output_filtering_contigs
    )
    prokka_run[i].add_inputs(output_filtering_contigs)
    prokka_run[i].add_outputs(
        str(srr_id) + "_prokka_output/" + str(srr_id) + ".gff",
        stage_out=True,
    )
    prokka_run[i].add_outputs(str(srr_id) + "_prokka_output.tar.gz", stage_out=True)
    # prokka_run[i].add_profiles(Namespace.PEGASUS, "label", str(srr_id))
    prokka_run[i].add_profiles(Namespace.PEGASUS, "runtime", "14400")
    prokka_run[i].add_profiles(Namespace.GLOBUS, "maxwalltime", "240")
    dax.add_jobs(prokka_run[i])
    # add files
    f = File(str(srr_id) + "_prokka_output/" + str(srr_id) + ".gff")
    list_of_gff_files.append(f)

    # add job for plasmidfinder
    plasmidfinder_run.append(Job("ex_plasmidfinder_run"))
    plasmidfinder_run[i].add_args(
        output_filtering_contigs, str(srr_id) + "_plasmidfinder_output"
    )
    plasmidfinder_run[i].add_inputs(output_filtering_contigs)
    plasmidfinder_run[i].add_outputs(
        str(srr_id) + "_plasmidfinder_output.tar.gz", stage_out=True
    )
    # plasmidfinder_run[i].add_profiles(Namespace.PEGASUS, "label", str(srr_id))
    dax.add_jobs(plasmidfinder_run[i])

    # add job for sistr
    sistr_run.append(Job("ex_sistr_run"))
    sistr_run[i].add_args(
        str(srr_id) + "_allele_results.json",
        str(srr_id) + "_novel_alleles.fasta",
        str(srr_id) + "_cgmlst_profiles.csv",
        str(srr_id) + "_sistr_output.csv",
        output_filtering_contigs,
    )
    sistr_run[i].add_inputs(output_filtering_contigs)
    sistr_run[i].add_outputs(str(srr_id) + "_sistr_output.csv", stage_out=False)
    dax.add_jobs(sistr_run[i])
    list_of_sistr_files.append(str(srr_id) + "_sistr_output.csv")

    i = i + 1


# add job for mlst
mlst_run = Job("ex_mlst_run")
mlst_run.add_args(*list_of_contig_files)
for l in list_of_contig_files:
    mlst_run.add_inputs(l)
o = File("mlst_output.csv")
mlst_run.set_stdout(o, stage_out=True)
mlst_run.add_profiles(Namespace.PEGASUS, "runtime", "108000")
mlst_run.add_profiles(Namespace.GLOBUS, "maxwalltime", "1800")
# mlst_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id))
dax.add_jobs(mlst_run)

# add job for abricate vfdb
abricate_vfdb_run = Job("ex_abricate_run")
abricate_vfdb_run.add_args("vfdb", *list_of_contig_files)
for l in list_of_contig_files:
    abricate_vfdb_run.add_inputs(l)
o = File("sabricate_vfdb_output.csv")
abricate_vfdb_run.set_stdout(o, stage_out=True)
# abricate_vfdb_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(abricate_vfdb_run)

# add job for abricate argannot
abricate_argannot_run = Job("ex_abricate_run")
abricate_argannot_run.add_args("argannot", *list_of_contig_files)
for l in list_of_contig_files:
    abricate_argannot_run.add_inputs(l)
o = File("sabricate_argannot_output.csv")
abricate_argannot_run.set_stdout(o, stage_out=True)
# abricate_argannot_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(abricate_argannot_run)

# add job for abricate card
abricate_card_run = Job("ex_abricate_run")
abricate_card_run.add_args("card", *list_of_contig_files)
for l in list_of_contig_files:
    abricate_card_run.add_inputs(l)
o = File("sabricate_card_output.csv")
abricate_card_run.set_stdout(o, stage_out=True)
# abricate_card_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(abricate_card_run)

# add job for abricate ncbi
abricate_ncbi_run = Job("ex_abricate_run")
abricate_ncbi_run.add_args("ncbi", *list_of_contig_files)
for l in list_of_contig_files:
    abricate_ncbi_run.add_inputs(l)
o = File("sabricate_ncbi_output.csv")
abricate_ncbi_run.set_stdout(o, stage_out=True)
# abricate_ncbi_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(abricate_ncbi_run)

# add job for abricate plasmidfinder
abricate_plasmidfinder_run = Job("ex_abricate_run")
abricate_plasmidfinder_run.add_args("plasmidfinder", *list_of_contig_files)
for l in list_of_contig_files:
    abricate_plasmidfinder_run.add_inputs(l)
o = File("sabricate_plasmidfinder_output.csv")
abricate_plasmidfinder_run.set_stdout(o, stage_out=True)
# abricate_plasmidfinder_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(abricate_plasmidfinder_run)

# add job for abricate resfinder
abricate_resfinder_run = Job("ex_abricate_run")
abricate_resfinder_run.add_args("resfinder", *list_of_contig_files)
for l in list_of_contig_files:
    abricate_resfinder_run.add_inputs(l)
o = File("sabricate_resfinder_output.csv")
abricate_resfinder_run.set_stdout(o, stage_out=True)
# abricate_resfinder_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(abricate_resfinder_run)

# add job for Roary
roary_run = Job("ex_roary_run")
roary_run.add_args("roary_output", *list_of_gff_files)
for l in list_of_gff_files:
    roary_run.add_inputs(l)
roary_run.add_outputs("roary_output/core_gene_alignment.aln", stage_out=True)
roary_run.add_outputs("roary_output.tar.gz", stage_out=True)
roary_run.add_profiles(Namespace.PEGASUS, "runtime", "604800")
roary_run.add_profiles(Namespace.GLOBUS, "maxwalltime", "10080")
roary_run.add_profiles(Namespace.CONDOR, "request_memory", "970000")
roary_run.add_profiles(Namespace.CONDOR, "memory", "970000")
# roary_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(roary_run)

# add job for baps_run
# R script, wrapper
fastbaps_run = Job("ex_fastbaps_run")
fastbaps_output = File("fastbaps_baps.csv")
fastbaps_run.add_args("roary_output/core_gene_alignment.aln", fastbaps_output)
fastbaps_run.add_inputs("roary_output/core_gene_alignment.aln")
fastbaps_run.add_outputs(fastbaps_output, stage_out=True)
fastbaps_run.add_profiles(Namespace.CONDOR, "request_memory", "30000")
fastbaps_run.add_profiles(Namespace.GLOBUS, "maxmemory", "30000")
fastbaps_run.add_profiles(Namespace.PEGASUS, "memory", "30000")
# fastbaps_run.add_profiles(Namespace.PEGASUS, "label", str(srr_id)))
dax.add_jobs(fastbaps_run)

# ls
ls_run = Job("ex_ls")
ls_run.add_args(run_dir)
dax.add_jobs(ls_run)

# add job for cat sistr files
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


length = len(list_of_contig_files)
for i in range(0, length):
    # Add control-flow dependencies
    dax.add_dependency(plasmidfinder_run[i], children=[ls_run])
    dax.add_dependency(prokka_run[i], children=[roary_run])
    # COMMENT OUT THE LINE BELOW TO SKIP SISTR IF NON SALMONELLA ORGANISM IS USED!!!
    dax.add_dependency(sistr_run[i], children=[cat])
dax.add_dependency(mlst_run, children=[ls_run])
dax.add_dependency(abricate_argannot_run, children=[ls_run])
dax.add_dependency(abricate_card_run, children=[ls_run])
dax.add_dependency(abricate_ncbi_run, children=[ls_run])
dax.add_dependency(abricate_plasmidfinder_run, children=[ls_run])
dax.add_dependency(abricate_resfinder_run, children=[ls_run])
dax.add_dependency(abricate_vfdb_run, children=[ls_run])
dax.add_dependency(roary_run, children=[fastbaps_run])
dax.add_dependency(cat, children=[merge_sistr_run])


# Write the DAX to stdout
dax.add_transformation_catalog(tc)
dax.add_site_catalog(sc)
dax.add_replica_catalog(rc)
dax.write(sys.stdout)
