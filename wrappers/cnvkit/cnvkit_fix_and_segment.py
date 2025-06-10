######################################
# wrapper for rule: cnvkit_fix_and_segment
######################################
from snakemake.shell import shell
import os

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnv cnvkit fix and segment \n##\n")
f.close()

shell.executable("/bin/bash")

# fix cnvkit.py fix Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn Reference.cnn -o Sample.cnr
#               " --no-edge " + \ for WGS

command = "cnvkit.py fix" + \
            " " + snakemake.input.targetcoverage + \
            " " + snakemake.input.antitargetcoverage + \
            " " + snakemake.input.cnv_reference + \
            " -o " + snakemake.output.fix + \
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

# cnvkit.py segment Sample.cnr -o Sample.cns
#  " -t 1e-6 " + \ for WGS
command = "cnvkit.py segment" + \
            " " + snakemake.output.fix + \
            " -o " + snakemake.output.segments + \
            " -p " + str(snakemake.threads) + \
            " -m " + str(snakemake.params.method) + \
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

## Annotate with different GTF
# cnv_annotate.py [-o OUTPUT] annotate cnv_file

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = str(dir_path+"/scripts/")

if snakemake.params.annotate and os.path.exists(snakemake.params.gtf):
    shell("cp {snakemake.output.segments} {snakemake.output.segments}.noannot")
    shell("cp {snakemake.output.fix} {snakemake.output.fix}.noannot")
    #fix
    command = dir_path+"cnv_annotate.py" + \
              " -o " + snakemake.output.fix  + \
              " " + snakemake.params.gtf + \
              " " + snakemake.output.fix+ ".noannot " + \
              " >> " + log_filename + " 2>&1"

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    #segments
    command = dir_path+"cnv_annotate.py" + \
              " -o " + snakemake.output.segments + \
              " " + snakemake.params.gtf + \
              " " + snakemake.output.segments+ ".noannot " + \
              " >> " + log_filename + " 2>&1"

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
