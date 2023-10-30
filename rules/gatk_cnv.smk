
rule gatk_cnv_collect_allelic_counts:
    input:
        bam=get_bam_input,
        interval=expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
        ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output:
        "variant_calls/{sample}/gatk_cnv/{tumor_normal}_clean.allelicCounts.tsv",
    params:
        extra="",
    log:
        "logs/{sample}/gatk_cnv/collect_allelic_counts_{tumor_normal}.log",
    threads: 8
    resources: mem=10
    conda:
        "../wrappers/gatk/env.yaml"
    shell:
        "(gatk --java-options '-Xmx8g' CollectAllelicCounts "
        "-L {input.interval} "
        "-I {input.bam} "
        "-R {input.ref} "
        "-O {output}"
        "{params.extra}) &> {log}"

rule gatk_cnv_collect_read_counts:
    input:
        bam=get_bam_input,
        interval=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    output:
        "variant_calls/{sample}/gatk_cnv/{tumor_normal}_read_counts.hdf5",
    params:
        mergingRule="OVERLAPPING_ONLY",
        extra="",
    log:
        "logs/{sample}/gatk_cnv/collect_read_counts_{tumor_normal}.log",
    threads: 8
    resources: mem=10
    conda:
        "../wrappers/gatk/env.yaml"
    shell:
        "(gatk --java-options '-Xmx8g' CollectReadCounts "
        "-I {input.bam} "
        "-L {input.interval} "
        "--interval-merging-rule {params.mergingRule} "
        "{params.extra} "
        "-O {output}) &> {log}"

def normal_read_counts_input(wildcards):
    if config["tumor_normal_paired"] == True:
        return expand("variant_calls/{sample}/gatk_cnv/normal_read_counts.hdf5",sample=sample_tab.loc[
                sample_tab.tumor_normal == "normal", "donor"].tolist())
    else:
        return expand("variant_calls/{sample}/gatk_cnv/tumor_read_counts.hdf5",sample=sample_tab.sample.tolist())


rule gatk_create_panel_of_normals:
    input:
        germinal_read_counts = normal_read_counts_input,
    output:
        hdf5PoN="variant_calls/all_samples/gatk_cnv/panel_of_normals.hdf5"
    log:
        "logs/all_samples/gatk_cnv/create_panel_of_normals.log",
    threads: 8
    resources: mem=10
    conda:
        "../wrappers/gatk/env.yaml"
    script:
        "../wrappers/gatk/create_panel_of_normals.py"


rule gatk_cnv_denoise_read_counts:
    input:
        hdf5PoN="variant_calls/all_samples/gatk_cnv/panel_of_normals.hdf5",
        hdf5Tumor="variant_calls/{sample}/gatk_cnv/tumor_read_counts.hdf5",
    output:
        denoisedCopyRatio="variant_calls/{sample}/gatk_cnv/clean.denoisedCR.tsv",
        stdCopyRatio="variant_calls/{sample}/gatk_cnv/clean.standardizedCR.tsv",
    params:
        extra="",
    log:
        "logs/{sample}/gatk_cnv/denoiseCR.log",
    threads: 8
    resources: mem=10
    conda:
        "../wrappers/gatk/env.yaml"
    shell:
        "(gatk --java-options '-Xmx8g' DenoiseReadCounts -I {input.hdf5Tumor} "
        "--count-panel-of-normals {input.hdf5PoN} "
        "--standardized-copy-ratios {output.stdCopyRatio} "
        "--denoised-copy-ratios {output.denoisedCopyRatio} "
        "{params.extra}) &> {log}"


rule gatk_cnv_model_segments:
    input:
        denoisedCopyRatio="variant_calls/{sample}/gatk_cnv/clean.denoisedCR.tsv",
        allelicCounts="variant_calls/{sample}/gatk_cnv/tumor_clean.allelicCounts.tsv",
    output:
        "variant_calls/{sample}/gatk_cnv/clean.modelFinal.seg",
        temp("variant_calls/{sample}/gatk_cnv/clean.cr.seg"),
        temp("variant_calls/{sample}/gatk_cnv/clean.af.igv.seg"),
        temp("variant_calls/{sample}/gatk_cnv/clean.cr.igv.seg"),
        temp("variant_calls/{sample}/gatk_cnv/clean.hets.tsv"),
        temp("variant_calls/{sample}/gatk_cnv/clean.modelBegin.cr.param"),
        temp("variant_calls/{sample}/gatk_cnv/clean.modelBegin.af.param"),
        temp("variant_calls/{sample}/gatk_cnv/clean.modelBegin.seg"),
        temp("variant_calls/{sample}/gatk_cnv/clean.modelFinal.af.param"),
        temp("variant_calls/{sample}/gatk_cnv/clean.modelFinal.cr.param"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        outprefix="clean",
        extra="",
    log:
        "logs/{sample}/gatk_cnv/modelFinal.log",
    threads: 8
    resources: mem=10
    conda:
        "../wrappers/gatk/env.yaml"
    shell:
        "(gatk --java-options '-Xmx8g' ModelSegments "
        "--denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} "
        "--output {params.outdir} "
        "--output-prefix {params.outprefix}"
        "{params.extra}) &> {log}"

rule gatk_cnv_call_copy_ratio_segments:
    input:
        "variant_calls/{sample}/gatk_cnv/clean.cr.seg",
    output:
        segments="variant_calls/{sample}/gatk_cnv/clean.calledCNVs.seg",
        igv_segments="variant_calls/{sample}/gatk_cnv/clean.calledCNVs.igv.seg",
    params:
        extra="",
    log:
        "logs/{sample}/gatk_cnv/calledCNVs.seg.log",
    threads: 8
    resources: mem=10
    conda:
        "../wrappers/gatk/env.yaml"
    shell:
        "(gatk --java-options '-Xmx8g' CallCopyRatioSegments "
        "--input {input} "
        "--output {output.segments} "
        "{params.extra}) &> {log}"

rule gatk_cnv_vcf:
    input:
        segment="variant_calls/{sample}/gatk_cnv/clean.modelFinal.seg",
    output:
        vcf="variant_calls/{sample}/gatk_cnv/result_SV.vcf",
    params:
        sample_id="{sample}",
        hom_del_limit=0.5, #dat do workflow.json
        het_del_limit=1.5, #dat do workflow.json
        dup_limit=2.5, #dat do workflow.json
        TC=0.5,#lambda wildcards: get_sample(samples, wildcards)["TC"],
    log:
        "logs/{sample}/gatk_cnv/convert_to_vcf.log",
    threads: 8
    resources: mem=10
    # conda:
    #     "../wrappers/gatk/env_python.yaml"
    script:
        "../wrappers/gatk/gatk_cnv_vcf.py"