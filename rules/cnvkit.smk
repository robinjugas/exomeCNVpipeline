
def get_bam_input(wildcards):
    if config["tumor_normal_paired"] == True:
        return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample), "sample"])[0],
    else:
        return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample)[0]


rule cnvkit_prepare_region_beds:
    input:
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output:
        reference_bed="variant_calls/cnvkit_prepare_reference/reference_bed.bed",
        target="variant_calls/all_samples/cnvkit/target.bed",
        antitarget="variant_calls/all_samples/cnvkit/antitarget.bed"
    params:
        normal_bams="mapped/*.bam",
        target=expand("{lib_ROI}.target.bed",lib_ROI=config["lib_ROI"]),
        antitarget=expand("{lib_ROI}.antitarget.bed",lib_ROI=config["lib_ROI"])
    log:
        "logs/cnvkit_prepare_reference/prepare_reference_1.log",
    threads: workflow.cores
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py access {input.reference} -o {output.reference_bed} &> {log}
        cnvkit.py autobin {params.normal_bams} -t {input.regions} -g {output.reference_bed} &>> {log}    
        echo "$PWD" 
        mv $PWD/{params.target} {output.target}
        mv $PWD/{params.antitarget} {output.antitarget}
        """

rule cnvkit_get_coverage:
    input:
        bam=get_bam_input,
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        target="variant_calls/all_samples/cnvkit/target.bed",
        antitarget="variant_calls/all_samples/cnvkit/antitarget.bed"
    output:
        targetcoverage="variant_calls/{sample}/cnvkit/{tumor_normal}.targetcoverage.cnn",
        antitargetcoverage="variant_calls/{sample}/cnvkit/{tumor_normal}.antitargetcoverage.cnn",
    log:
        "logs/{sample}/cnvkit/{tumor_normal}_get_coverage.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py coverage -p {threads} {input.bam} {input.target} -o {output.targetcoverage} &> {log}
        cnvkit.py coverage -p {threads} {input.bam} {input.antitarget} -o {output.antitargetcoverage} &>> {log}
        """

def normal_coverage_inputs(wildcards):
    if config["tumor_normal_paired"] == True: #normal reference
        return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/normal.{tag}targetcoverage.cnn",sample=sample_tab.loc[
            sample_tab.tumor_normal == "normal", "donor"].tolist(),tag = ["","anti"]))}
    else:
        return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/tumor.{tag}targetcoverage.cnn",sample=sample_tab.sample.tolist(),tag = ["","anti"]))}



rule cnvkit_prepare_reference:
    input:
        unpack(normal_coverage_inputs),
        reference=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output:
        reference_cnn="variant_calls/all_samples/cnvkit/normal_reference.cnn",
    params:
        scope=config["lib_ROI"]
    log:
        "logs/all_samples/cnvkit_prepare_reference.log"
    threads: workflow.cores
    conda:
        "../wrappers/cnvkit/env.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_reference.py"


rule cnvkit_fix_and_segment:
    input:
        targetcoverage="variant_calls/{sample}/cnvkit/tumor.targetcoverage.cnn",
        antitargetcoverage="variant_calls/{sample}/cnvkit/tumor.antitargetcoverage.cnn",
        cnv_reference="variant_calls/all_samples/cnvkit/normal_reference.cnn",
    output:
        fix="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
        segments="variant_calls/{sample}/cnvkit/segmented_cov.cns",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        method="hybrid",
        extra="",
        scope=config["lib_ROI"],
        annotate=True, # annotate with GTF instead of BED
        gtf=expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir=reference_directory,ref_name=config["reference"])[0],
    log:
        "logs/{sample}/cnvkit/cnvkit_fix_and_segment.log"
    threads: 10
    resources: mem=8
    conda:
        "../wrappers/cnvkit/env.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_fix_and_segment.py"


rule vardict:
    input:  bam = "mapped/{bam_name}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
            regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: vcf="variant_calls/{sample}/cnvkit/vardict_SNV_{bam_name}.vcf",
    log: "logs/{sample}/cnvkit/{bam_name}_vardict.log"
    threads: 5
    resources: mem=8
    params:
        AF_threshold=0.05
    conda: "../wrappers/vardict/env.yaml"
    script: "../wrappers/vardict/script.py"


def vardict_SNV_vcf_input(wildcards):
    if config["tumor_normal_paired"] == True:
        return expand("variant_calls/{sample}/cnvkit/vardict_SNV_{input_bam}.vcf",sample = wildcards.sample,input_bam=sample_tab.loc[(sample_tab["donor"] == wildcards.sample) & (sample_tab["tumor_normal"] == "normal"), "sample"])[0]
    else:
        return expand("variant_calls/{sample}/cnvkit/vardict_SNV_{sample}.vcf",sample = wildcards.sample)[0]


rule cnvkit_call: #method parameter ; threshold or clonal https://cnvkit.readthedocs.io/en/stable/germline.html
    input:
        vcf = vardict_SNV_vcf_input,
        segment="variant_calls/{sample}/cnvkit/segmented_cov.cns",
    output:
        calls="variant_calls/{sample}/cnvkit/CNV_calls.cns",
    params:
        TC=0.5, #lambda wildcards: sample_tab.loc[wildcards.sample, 'donor'], #tumor content?
        scope=config["lib_ROI"],
        method="threshold"
    log:
        "logs/{sample}/cnvkit/cnvkit_call.log"
    threads: 10
    resources: mem=10
    conda: "../wrappers/cnvkit/env.yaml"
    script: "../wrappers/cnvkit/cnvkit_call.py"
    # shell:
    #     "(cnvkit.py call -y -m clonal {input.segment} -v {input.vcf} -o {output.calls} --purity {params.TC} {params.extra}) &> {log}"


# we can try genemetrics both with and without the segment files, take the intersection of those as a list of “trusted” genes, and visualize each of them with scatter:
rule cnvkit_genemetrics:
    input:
        cns="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
    output:
        trusted_genes="variant_calls/{sample}/cnvkit/trusted-genes.txt",
        ratio_genes="variant_calls/{sample}/cnvkit/ratio-genes.txt",
        segment_genes="variant_calls/{sample}/cnvkit/segment-genes.txt"
    params:
        threshold=0.2,
        min_probes=3,
        working_directory="variant_calls/{sample}/cnvkit/"
    log:
        "logs/{sample}/cnvkit/cnvkit_genemetrics.log"
    threads: 1
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        (cnvkit.py genemetrics -y {input.cnr} -s {input.cns} -t {params.threshold} -m {params.min_probes} | tail -n+2 | cut -f1 | sort > {output.segment_genes}) &> {log}
        (cnvkit.py genemetrics -y {input.cnr} | tail -n+2 | cut -f1 | sort > {output.ratio_genes}) &> {log}
        comm -12 {output.ratio_genes} {output.segment_genes} > {output.trusted_genes}
        """

        # for gene in `cat {output.trusted_genes}`
        # do
        #     cnvkit.py scatter -s {wildcards.sample}.cn{s,r} -g $gene -o Sample-$gene-scatter.pdf
        # done

# cnv annotate


########################################################################################################################
# PLOTS
def cnvkit_heatmap_input(wildcards):
    return expand("variant_calls/{sample}/cnvkit/CNV_calls.cns", sample=sample_tab.sample)


rule cnvkit_heatmap:
    input:
        cns=cnvkit_heatmap_input
    output:
        plot="variant_calls/all_samples/cnvkit_heatmap.png",
    params:
        extra="-d",
    log:
        "logs/all_samples/cnvkit_heatmap.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py heatmap {input.cns} -o {output.plot} {params.extra}) &> {log}"

rule cnvkit_diagram:
    input:
        cns="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
    output:
        pdf="variant_calls/{sample}/cnvkit/cnvkit_diagram.pdf",
    params:
        extra="",
    log:
        "logs/{sample}/cnvkit/cnvkit_diagram.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} {params.extra}) &> {log}"

rule cnvkit_scatter:
    input:
        vcf = vardict_SNV_vcf_input,
        cns="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
    output:
        plot="variant_calls/{sample}/cnvkit/cnvkit_scatter.png",
    params:
        extra="",
    log:
        "logs/{sample}/cnvkit/cnvkit_scatter.log"
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py scatter {input.cnr} -s {input.cns} -v {input.vcf} -o {output.plot} {params.extra}) &> {log}"


rule cnvkit_convert_to_vcf:
    input:
        segment="variant_calls/{sample}/cnvkit/CNV_calls.cns",
    output:
        vcf="variant_calls/{sample}/cnvkit/result_SV.vcf",
    params:
        sample="{sample}",
        hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
    log:
        "logs/{sample}/cnvkit/convert_to_vcf.log",
    threads: 10
    resources: mem=10
    conda:
        "../wrappers/cnvkit/env_python.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_vcf.py"


rule cnvkit_rename:
    input:
        segment="variant_calls/{sample}/cnvkit/CNV_calls.cns"
    output:
        tsv="variant_calls/{sample}/cnvkit/cnvkit_CNV_{sample}.tsv"
    threads: 1
    shell:
        "cp {input.segment} {output.tsv}"

