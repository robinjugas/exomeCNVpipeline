######################################
# wrapper for rule: cnvAnnotateTargetRegions
######################################

import os
from snakemake.shell import shell

input_vcf = snakemake.input.vcf
output_tsv = snakemake.output.tsv
output_bed = snakemake.output.bed

script = os.path.join(os.path.dirname(__file__), "vcf_to_tsv.py")

shell(f"{script} {input_vcf} {output_tsv} {output_bed}")
