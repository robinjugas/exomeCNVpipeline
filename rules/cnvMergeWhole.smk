
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
        "../wrappers/cnvMergeWhole/env.yaml"
    script:
        "../wrappers/cnvMergeWhole/script.py"



snake_dir = workflow.basedir
rule classifyCNV_whole:
    input:
        bed="CNV_Whole/{sample}.callers_merged.bed",
    output:
        txt="CNV_Whole/{sample}.classified.txt",
        dir=directory("CNV_Whole/{sample}/classifyCNV/")
    params:
        GenomeBuild=str(config["GenomeBuild"]),
    log:
        "logs/{sample}/{sample}_classifyCNV_TargetRegions.log"
    threads: 3
    conda:
        "../wrappers/classifyCNV/env.yaml"
    shell:
        """
        python3 {snake_dir}/wrappers/classifyCNV/ClassifyCNV/ClassifyCNV.py --infile {input.bed} --outdir {output.dir} --GenomeBuild {params.GenomeBuild} --precise
        cp {output.dir}/Scoresheet.txt {output.txt}
        """

rule cnvAnnotateWhole:
    input:
        tsv="CNV_Whole/{sample}.callers_merged.tsv",
        classifyCNV_txt="CNV_Whole/{sample}.classified.txt",
        gtf="processed_GFT.tsv",
    output:
        tsv="mergedCNVs_final/{sample}_final_CNVs_annotated.tsv",
        xlsx="mergedCNVs_final/{sample}_final_CNVs_annotated.xlsx",
    log:
        "logs/{sample}/{sample}_cnvAnnotateWhole.log"
    threads: 3
    conda:
        "../wrappers/cnvAnnotateWhole/env.yaml"
    script:
        "../wrappers/cnvAnnotateWhole/script.py"