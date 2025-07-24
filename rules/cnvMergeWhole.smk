
rule cnvWholeMerge:
    input:
        tsv=expand("variant_calls/{{sample}}/{cnv_caller}/{cnv_caller}_CNV_{{sample}}.tsv",cnv_caller=used_SV_callers),
        bed=config["bed_file"]
    output:
        bed="CNV_Whole/{sample}.callers_merged.bed",
        tsv="CNV_Whole/{sample}.callers_merged.tsv"
    log:
        "logs/{sample}/{sample}_cnvWholeMerge.log"
    threads: 3
    conda:
        "../wrappers/cnvWholeMerge/env.yaml"
    script:
        "../wrappers/cnvWholeMerge/script.py"



snake_dir = workflow.basedir
rule classifyCNV_whole:
    input:
        bed="CNV_Whole/{sample}.callers_merged.bed",
    output:
        txt="CNV_Whole/{sample}.classified.txt",
        dir=directory("CNV_Whole/{sample}/classifyCNV/")
    log:
        "logs/{sample}/{sample}_classifyCNV_TargetRegions.log"
    threads: 3
    conda:
        "../wrappers/classifyCNV/env.yaml"
    shell:
        """
        python3 {snake_dir}/wrappers/classifyCNV/ClassifyCNV/ClassifyCNV.py --infile {input.bed} --outdir {output.dir} --GenomeBuild hg19 --precise
        cp {output.dir}/Scoresheet.txt {output.txt}
        """

rule cnvAnnotateWhole:
    input:
        tsv="CNV_Whole/{sample}.callers_merged.tsv",
        classifyCNV_txt="CNV_Whole/{sample}.classified.txt"
    output:
        tsv="mergedCNVs_final/{sample}_called_CNVs.tsv"
    log:
        "logs/{sample}/{sample}_cnvAnnotateWhole.log"
    threads: 3
    conda:
        "../wrappers/cnvWholeMerge/env.yaml"
    script:
        "../wrappers/cnvAnnotateWhole/script.py"