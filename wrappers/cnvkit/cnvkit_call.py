######################################
# wrapper for rule: cnvkit_call
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnvkit call \n##\n")
f.close()

shell.executable("/bin/bash")


# cnvkit.py call Sample.cns -y -v Sample.vcf -o Sample.call.cns
if snakemake.params.analysis_mode == "somatic_paired" :
    command = "cnvkit.py call " + \
              snakemake.input.segment + \
              " -m " + str(snakemake.params.method) + \
              " -v " + snakemake.input.vcf + \
              " -o " + snakemake.output.calls + \
              " >> " + log_filename + " 2>&1"
else:
    command = "cnvkit.py call " + \
                snakemake.input.segment + \
                " -m " + str(snakemake.params.method) + \
                " -o " + snakemake.output.calls + \
                " >> " + log_filename + " 2>&1"
 # " --male-reference " + \



f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
