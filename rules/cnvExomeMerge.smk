
rule cnvMerge:
    input:
        tsv=expand("variant_calls/{{sample}}/{cnv_caller}/{cnv_caller}_CNV_{{sample}}.tsv",cnv_caller=used_SV_callers),
    output:
        bed="CNV_merged/{sample}.callers_merged.bed",
        tsv="CNV_merged/{sample}.callers_merged.tsv"
    log:
        "logs/{sample}/cnvMerge/cnvMerge.log"
    threads: 6
    conda:
        "../wrappers/cnvExonMerge/env.yaml" #stejne jako u cnvExonMerge
    script:
        "../wrappers/cnvMerge/script.py"



snake_dir = workflow.basedir
rule classifyCNV:
    input:
        bed="CNV_merged/{sample}.callers_merged.bed",
    output:
        txt="CNV_merged/{sample}.classified.txt",
        dir=directory("CNV_merged/{sample}/")
    log:
        "logs/{sample}/classifyCNV/classifyCNV.log"
    threads: 6
    conda:
        "../wrappers/classifyCNV/env.yaml"
    shell:
        """
        python3 {snake_dir}/wrappers/classifyCNV/ClassifyCNV/ClassifyCNV.py --infile {input.bed} --outdir {output.dir} --GenomeBuild hg19 --precise
        cp {output.dir}/Scoresheet.txt {output.txt}
        """

rule cnvAnnotate:
    input:
        tsv="CNV_merged/{sample}.callers_merged.tsv",
        txt="CNV_merged/{sample}.classified.txt"
    output:
        tsv="CNV_merged/{sample}_final.tsv"
    log:
        "logs/{sample}/cnvAnnotate/cnvAnnotate.log"
    threads: 6
    conda:
        "../wrappers/cnvAnnotate/env.yaml"
    script:
        "../wrappers/cnvAnnotate/script.py"