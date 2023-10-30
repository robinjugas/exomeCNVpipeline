######################################
# wrapper for rule: vardict
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: vardict \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'a+')
f.write("## VERSION: vardict-java "+version+"\n")
f.close()

#AttributeError: 'Wildcards' object has no attribute 'full_name'
# replace 'full_name' with fullname

TEST_STRAND_BIAS = os.path.abspath(os.path.dirname(__file__)) + "/teststrandbias.R"
VAR2VCF = os.path.abspath(os.path.dirname(__file__)) + "/var2vcf_valid.pl"

command = "vardict-java -G " + snakemake.input.ref + \
          " -th " + str(snakemake.threads) + \
          " -N " + snakemake.wildcards.sample_name + \
          " -b " + snakemake.input.bam + \
          " -c 1 -S 2 -E 3 -g 4 " + snakemake.input.regions + \
          " 2>> " + log_filename + \
          " | " + TEST_STRAND_BIAS + \
          " 2>> " + log_filename + \
          " | " + VAR2VCF + \
          " -m 7 -c 1 -N " + snakemake.wildcards.sample_name + \
          " -f " + str(snakemake.params.AF_threshold) + " > " + snakemake.output.vcf + \
          " 2>> " + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "sed -i '/TYPE=[DI][UNE][PVL]/d' " + snakemake.output.vcf
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)


