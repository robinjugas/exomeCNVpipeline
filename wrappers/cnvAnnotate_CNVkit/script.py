######################################
# wrapper for rule: cnvAnnotateCNVkit
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## cnvAnnotateCNVkit \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript " + os.path.abspath(os.path.dirname(__file__))+"/cnvAnnotateCNVkit.R " +\
    snakemake.input.tsv + " " + \
    snakemake.input.classifyCNV_txt + " " + \
    snakemake.input.gtf + " " + \
    snakemake.output.tsv + " " + \
    snakemake.output.xlsx + " " + \
    " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()


shell(command)
