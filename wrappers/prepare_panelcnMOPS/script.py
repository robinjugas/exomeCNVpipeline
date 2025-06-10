#############################################################
# wrapper for rule: prepare_panelcnMOPS
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: prepare_panelcnMOPS \n##\n")
f.close()

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/prepare_panelCNMOPS_wrapper.R "+\
            snakemake.input.bed + " " +\
            snakemake.output.Rdata + " " +\
            " ".join(snakemake.input.bams) +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)
