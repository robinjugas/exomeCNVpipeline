rule cnvExonMerge:
    input:
        tsv=expand("variant_calls/{{sample}}/{cnv_caller}/{cnv_caller}_CNV_{{sample}}.tsv",cnv_caller=used_SV_callers),
        gtf="/home/rj/4TB/CEITEC/GTFs_GRCh37/gencode.v41lift37.basic.annotation.gtf"  #upravit
    output:
        bed="CNV_exon_merged/{sample}.callers_merged.bed",
        tsv="CNV_exon_merged/{sample}.callers_merged.tsv"
    log:
        "logs/{sample}/cnvExonMerge/cnvExonMerge.log"
    threads: 6
    conda:
        "../wrappers/cnvExonMerge/env.yaml"
    script:
        "../wrappers/cnvExonMerge/script.py"


snake_dir = workflow.basedir
rule classifyCNV:
    input:
        bed="CNV_exon_merged/{sample}.callers_merged.bed",
    output:
        txt="CNV_exon_merged/{sample}.classified.txt",
        dir=directory("CNV_exon_merged/{sample}/")
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

rule cnvAnnotateExons:
    input:
        tsv="CNV_exon_merged/{sample}.callers_merged.tsv",
        txt="CNV_exon_merged/{sample}.classified.txt"
    output:
        tsv="CNV_exon_merged/{sample}_final.tsv"
    log:
        "logs/{sample}/cnvAnnotateExons/cnvAnnotateExons.log"
    threads: 6
    conda:
        "../wrappers/cnvAnnotateExons/env.yaml"
    script:
        "../wrappers/cnvAnnotateExons/script.py"