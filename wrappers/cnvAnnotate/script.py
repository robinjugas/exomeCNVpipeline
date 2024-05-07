######################################
# wrapper for rule: cnvExonAnnotate
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## cnvExonAnnotate \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript " + os.path.abspath(os.path.dirname(__file__))+"/cnvExonAnnotate.R " +\
    snakemake.input.tsv + " " + \
    snakemake.input.txt + " " + \
    snakemake.output.tsv + " " + \
    " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()


shell(command)
