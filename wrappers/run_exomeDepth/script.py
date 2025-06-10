#############################################################
# wrapper for rule: exomeDepth
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: exomeDepth \n##\n")
f.close()


command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/ExomeDepth_wrapper.R "+\
            snakemake.input.bam + " " +\
            snakemake.input.cohort + " " +\
            snakemake.params.sampleName + " " +\
            snakemake.input.ref + " " +\
            snakemake.output.tsv +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)
