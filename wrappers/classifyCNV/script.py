######################################
# wrapper for rule: cnvExonMerge
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## cnvExonMerge \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript " + os.path.abspath(os.path.dirname(__file__))+"/cnvExonMerge.R " +\
    snakemake.output.tsv + " " + \
    snakemake.input.gtf + " " + \
    str(snakemake.params.min_sv_length) + " " + \
    str(snakemake.params.max_sv_length) + " " + \
    str(snakemake.params.mergeCloseEventsDistanceThreshold) + " " + \
    " ".join(snakemake.input.tsv) +\
    " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()


shell(command)
