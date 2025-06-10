rule cnvTargetRegionsMerge:
    input:
        tsv=expand("variant_calls/{{sample}}/{cnv_caller}/{cnv_caller}_CNV_{{sample}}.tsv",cnv_caller=used_SV_callers),
        gtf=config["gtf_file"],
        bed=config["bed_file"]
    output:
        bed="CNV_TargetRegions/{sample}.callers_merged.bed",
        tsv="CNV_TargetRegions/{sample}.callers_merged.tsv"
    log:
        "logs/{sample}/{sample}_cnvTargetRegionsMerge.log"
    threads: 6
    conda:
        "../wrappers/cnvTargetRegionsMerge/env.yaml"
    script:
        "../wrappers/cnvTargetRegionsMerge/script.py"


snake_dir = workflow.basedir
rule classifyCNV_exons:
    input:
        bed="CNV_TargetRegions/{sample}.callers_merged.bed",
    output:
        txt="CNV_TargetRegions/{sample}.classified.txt",
        dir=directory("CNV_TargetRegions/{sample}/classifyCNV/")
    log:
        "logs/{sample}/{sample}_classifyCNV_TargetRegions.log"
    threads: 6
    conda:
        "../wrappers/classifyCNV/env.yaml"
    shell:
        """
        python3 {snake_dir}/wrappers/classifyCNV/ClassifyCNV/ClassifyCNV.py --infile {input.bed} --outdir {output.dir} --GenomeBuild hg19 --precise
        cp {output.dir}/Scoresheet.txt {output.txt}
        """

rule cnvAnnotateTargetRegions:
    input:
        tsv="CNV_TargetRegions/{sample}.callers_merged.tsv",
        classifyCNV_txt="CNV_TargetRegions/{sample}.classified.txt"
    output:
        tsv="CNV_TargetRegions/{sample}_called_TargetRegions.tsv"
    log:
        "logs/{sample}/{sample}_cnvAnnotateTargetRegions.log"
    threads: 6
    conda:
        "../wrappers/cnvTargetRegionsMerge/env.yaml"
    script:
        "../wrappers/cnvAnnotateTargetRegions/script.py"