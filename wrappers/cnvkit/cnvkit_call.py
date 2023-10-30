######################################
# wrapper for rule: cnvkit_call
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnvkit call \n##\n")
f.close()

shell.executable("/bin/bash")

# call -y -m clonal {input.segment} -v {input.vcf} -o {output.calls} --purity {params.TC} {params.extra}) &> {log}
# " --purity " + snakemake.params.TC + \


if snakemake.params.scope == "wgs" or snakemake.params.scope == "WGS":
    command = "cnvkit.py call " + \
              snakemake.input.segment + \
              " --male-reference " + \
              " -m clonal " + \
              " -v " + snakemake.input.vcf + \
              " -o " + snakemake.output.calls + \
              " >> " + log_filename + " 2>&1"
else:
    command = "cnvkit.py call " + \
              snakemake.input.segment + \
              " --male-reference " + \
              " -m clonal " + \
              " -v " + snakemake.input.vcf + \
              " -o " + snakemake.output.calls + \
              " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
