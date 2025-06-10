######################################
# wrapper for rule: cnvkit_prepare_reference
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnv reference \n##\n")
f.close()

shell.executable("/bin/bash")

# –no-edge if WGS


if snakemake.params.analysis_mode=="somatic_paired":
        command = "cnvkit.py reference " + " ".join(snakemake.input.normal_coverage_inputs) + \
                  " --fasta " + snakemake.input.reference + \
                  " -o " + snakemake.output.reference_cnn + \
                  " >> " + log_filename + " 2>&1"


if snakemake.params.analysis_mode=="reference_cohort":
        command = "cnvkit.py reference " + " ".join(snakemake.input.normal_coverage_inputs) + \
                  " --fasta " + snakemake.input.reference + \
                  " -o " + snakemake.output.reference_cnn + \
                  " >> " + log_filename + " 2>&1"


if snakemake.params.analysis_mode=="baseline":
        command = "cnvkit.py reference " + " ".join(snakemake.input.normal_coverage_inputs) + \
                  " --fasta " + snakemake.input.reference + \
                  " -o " + snakemake.output.reference_cnn + \
                  " >> " + log_filename + " 2>&1"

# else: 
#      command = "cnvkit.py reference " + \
#                   " --male-reference " + \
#                   " --fasta " + snakemake.input.reference + \
#                   " -o " + snakemake.output.reference_cnn + \
#                   " -t " + snakemake.input.target + \
#                   " -a " + snakemake.input.antitarget + \
#                   " >> " + log_filename + " 2>&1"
#                 #" --no-edge " + \ # –no-edge if WGS
# https://cnvkit.readthedocs.io/en/stable/pipeline.html#reference
# or cnvkit.py reference -o FlatReference.cnn -f ucsc.hg19.fa -t targets.bed -a antitargets.bed


f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)


# if snakemake.input.normal_coverage_inputs:
#     if snakemake.params.scope == "wgs" or snakemake.params.scope == "WGS":
#         command = "cnvkit.py reference " + " ".join(snakemake.input.normal_coverage_inputs) + \
#                   " --fasta " + snakemake.input.reference + \
#                   " -o " + snakemake.output.reference_cnn + \
#                   " --no-edge " + \
#                   " >> " + log_filename + " 2>&1"
#     else:
#         command = "cnvkit.py reference " + " ".join(snakemake.input.normal_coverage_inputs) + \
#                   " --fasta " + snakemake.input.reference + \
#                   " -o " + snakemake.output.reference_cnn + \
#                   " >> " + log_filename + " 2>&1"

# else:
#     if snakemake.params.scope == "wgs" or snakemake.params.scope == "WGS":
#         command = "cnvkit.py reference " + \
#                   " --male-reference " + \
#                   " --fasta " + snakemake.input.reference + \
#                   " -o " + snakemake.output.reference_cnn + \
#                   " -t " + snakemake.input.target + \
#                   " -a " + snakemake.input.antitarget + \
#                   " --no-edge " + \
#                   " >> " + log_filename + " 2>&1"
#     else:
#         command = "cnvkit.py reference " + \
#                   " --male-reference " + \
#                   " --fasta " + snakemake.input.reference + \
#                   " -o " + snakemake.output.reference_cnn + \
#                   " -t " + snakemake.input.target + \
#                   " -a " + snakemake.input.antitarget + \
#                   " >> " + log_filename + " 2>&1"