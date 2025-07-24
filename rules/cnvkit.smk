###########################################################################################################################################

# Prepare a BED file of baited regions for use with CNVkit.
rule cnvkit_prepare_region_beds:
    input:
        reference=config["genome"],
        regions=config["bed_file"],
        normal_bams = cohort_bam_inputs,
    output:
        reference_bed="variant_calls/all_samples/reference_bed.bed",
        target="variant_calls/all_samples/target.bed",
        antitarget="variant_calls/all_samples/antitarget.bed"
    params:
        #normal_bams="mapped/*.bam",
        target=re.sub(r"\.bed$", ".target.bed", os.path.basename(config["bed_file"])),
        antitarget=re.sub(r"\.bed$", ".antitarget.bed", os.path.basename(config["bed_file"])),
    log:
        "logs/cnvkit_prepare_reference/prepare_reference_1.log", 
    threads: 20
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py access {input.reference} -o {output.reference_bed}.temp &> {log}
        awk 'NR==FNR {{chr[$1]; next}} $1 in chr' {input.regions} {output.reference_bed}.temp > {output.reference_bed}
        cnvkit.py autobin {input.normal_bams} -m hybrid -t {input.regions} -g {output.reference_bed} &>> {log}    
        echo "$PWD" 
        mv $PWD/{params.target} {output.target}
        mv $PWD/{params.antitarget} {output.antitarget}
        """
##        cnvkit.py antitarget {input.regions} -g {output.reference_bed} -o {params.antitarget}

# test
# awk 'NR==FNR {{chr[$1]; next}} $1 in chr' OvarianCancer_GRCh38.bed reference_bed.bed > reference_bed.bed.new


###########################################################################################################################################

# For each sample...
rule cnvkit_get_coverage:
    input:
        bam="mapped/{sample}.bam",
        reference=config["genome"],
        target="variant_calls/all_samples/target.bed",
        antitarget="variant_calls/all_samples/antitarget.bed"
    output:
        targetcoverage="variant_calls/{sample}/cnvkit/{status}.targetcoverage.cnn",
        antitargetcoverage="variant_calls/{sample}/cnvkit/{status}.antitargetcoverage.cnn",
    log:
        "logs/{sample}/{status}_get_coverage.log"  #"logs/{sample}/{sample}_cnMOPS.log"
    threads: 6
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py coverage -p {threads} {input.bam} {input.target} -o {output.targetcoverage} &> {log}
        cnvkit.py coverage -p {threads} {input.bam} {input.antitarget} -o {output.antitargetcoverage} &>> {log}
        """

###########################################################################################################################################

# To analyze a cohort sequenced on a single platform, we recommend combining all normal samples into a pooled reference, even if matched tumor-normal pairs were sequenced
#  – our benchmarking showed that a pooled reference performed slightly better than constructing a separate reference for each matched tumor-normal pair. 
# Furthermore, even matched normals from a cohort sequenced together can exhibit distinctly different copy number biases (see Plagnol et al. 2012 and Backenroth et al. 2014);
#  reusing a pooled reference across the cohort provides some consistency to help diagnose such issues.
def normal_coverage_inputs(wildcards):
    if config["analysis_mode"] == "somatic_paired":
         return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/normal.{tag}targetcoverage.cnn",sample=sample_tab.loc[
            sample_tab.status == "normal", "sample_name"].tolist(),tag = ["","anti"]))}
        # return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/normal.{tag}targetcoverage.cnn",sample=sample_tab.loc[
        #     sample_tab.tumor_normal == "normal", "donor"].tolist(),tag = ["","anti"]))}
    elif config["analysis_mode"] == "reference_cohort":
        return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/control.{tag}targetcoverage.cnn",sample=sample_tab.loc[
            sample_tab.status == "control", "sample_name"].tolist(),tag = ["","anti"]))}
    elif config["analysis_mode"] == "baseline":
        return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/normal.{tag}targetcoverage.cnn",sample=sample_tab.loc[
            sample_tab.status == "normal", "sample_name"].tolist(),tag = ["","anti"]))}
    else:
        return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/normal.{tag}targetcoverage.cnn",sample=sample_tab.loc[
            sample_tab.status == "normal", "sample_name"].tolist(),tag = ["","anti"]))}
    # if config["tumor_normal_paired"] == True: #normal reference
    #     return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/normal.{tag}targetcoverage.cnn",sample=sample_tab.loc[
    #         sample_tab.tumor_normal == "normal", "donor"].tolist(),tag = ["","anti"]))}
    # else:
    #     return {'normal_coverage_inputs': set(expand("variant_calls/{sample}/cnvkit/tumor.{tag}targetcoverage.cnn",sample=sample_tab.sample.tolist(),tag = ["","anti"]))}



# With all normal/reference samples...
rule cnvkit_prepare_reference:
    input:
        unpack(normal_coverage_inputs),
        reference=config["genome"],
    output:
        reference_cnn="variant_calls/all_samples/normal_reference.cnn",
    params:
        analysis_mode=config["analysis_mode"],
        scope=config["bed_file"]
    log:
        "logs/all_samples/cnvkit_prepare_reference.log"
    threads: 20
    conda:
        "../wrappers/cnvkit/env.yaml"
    script:
        "../wrappers/cnvkit/cnvkit_reference.py"

###########################################################################################################################################

# For each tumor/ sample...
if config["analysis_mode"] == "somatic_paired":
    rule cnvkit_fix_and_segment_somatic_paired:
        input:
            targetcoverage="variant_calls/{sample}/cnvkit/tumor.targetcoverage.cnn",
            antitargetcoverage="variant_calls/{sample}/cnvkit/tumor.antitargetcoverage.cnn",
            cnv_reference="variant_calls/all_samples/normal_reference.cnn",
        output:
            fix="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
            segments="variant_calls/{sample}/cnvkit/segmented_cov.cns",
        params:
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            method="cbs",
            annotate=True, # annotate with GTF instead of BED
            gtf=config["gtf_file"],
        log:
            "logs/{sample}/cnvkit_fix_and_segment.log"
        threads: 6
        conda:
            "../wrappers/cnvkit/env.yaml"
        script:
            "../wrappers/cnvkit/cnvkit_fix_and_segment.py"

if config["analysis_mode"] == "reference_cohort":
    rule cnvkit_fix_and_segment_reference_cohort:
        input:
            targetcoverage="variant_calls/{sample}/cnvkit/normal.targetcoverage.cnn",
            antitargetcoverage="variant_calls/{sample}/cnvkit/normal.antitargetcoverage.cnn",
            cnv_reference="variant_calls/all_samples/normal_reference.cnn",
        output:
            fix="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
            segments="variant_calls/{sample}/cnvkit/segmented_cov.cns",
        params:
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            method="cbs",
            annotate=True, # annotate with GTF instead of BED
            gtf=config["gtf_file"],
        log:
            "logs/{sample}/cnvkit_fix_and_segment.log"
        threads: 6
        conda:
            "../wrappers/cnvkit/env.yaml"
        script:
            "../wrappers/cnvkit/cnvkit_fix_and_segment.py"

if config["analysis_mode"] == "baseline":
    rule cnvkit_fix_and_segment_baseline:
        input:
            targetcoverage="variant_calls/{sample}/cnvkit/normal.targetcoverage.cnn",
            antitargetcoverage="variant_calls/{sample}/cnvkit/normal.antitargetcoverage.cnn",
            cnv_reference="variant_calls/all_samples/normal_reference.cnn",
        output:
            fix="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
            segments="variant_calls/{sample}/cnvkit/segmented_cov.cns",
        params:
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            method="cbs",
            annotate=True, # annotate with GTF instead of BED
            gtf=config["gtf_file"],
        log:
            "logs/{sample}/cnvkit_fix_and_segment.log"
        threads: 6
        conda:
            "../wrappers/cnvkit/env.yaml"
        script:
            "../wrappers/cnvkit/cnvkit_fix_and_segment.py"

###########################################################################################################################################
## for tumor-normal paired only, change the cnvkit-call rule script

def vardict_SNV_vcf_input(wildcards):
    return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="tumor", "sample_name"])[1]

rule vardict:
    input:  bam = vardict_SNV_vcf_input,
            ref=config["genome"],
            refdict=config["refdict"],
            regions=config["bed_file"],
    output: vcf="variant_calls/{sample}/cnvkit/{sample}_tumor_vardict_SNP.vcf",
    log: "logs/{sample}/vardict.log"
    threads: 6
    resources: mem=8
    params:
        AF_threshold=0.05
    conda: "../wrappers/vardict/env.yaml"
    script: "../wrappers/vardict/script.py"



###########################################################################################################################################

 #method parameter ; threshold or clonal https://cnvkit.readthedocs.io/en/stable/germline.html
if config["analysis_mode"] == "somatic_paired":
    rule cnvkit_call_tumor:
        input:
            vcf = "variant_calls/{sample}/cnvkit/{sample}_tumor_vardict_SNP.vcf",
            segment="variant_calls/{sample}/cnvkit/segmented_cov.cns",
        output:
            calls="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        params:
            analysis_mode=config["analysis_mode"],
            method="threshold"
        log:
            "logs/{sample}/cnvkit_call.log"
        threads: 6
        conda: "../wrappers/cnvkit/env.yaml"
        script: "../wrappers/cnvkit/cnvkit_call.py"
        # shell:
        #     "(cnvkit.py call -y -m clonal {input.segment} -v {input.vcf} -o {output.calls} --purity {params.TC} {params.extra}) &> {log}"

else:
    rule cnvkit_call_else:
        input:
            segment="variant_calls/{sample}/cnvkit/segmented_cov.cns",
        output:
            calls="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        params:
            analysis_mode=config["analysis_mode"],
            method="threshold"
        log:
            "logs/{sample}/cnvkit_call.log"
        threads: 6
        conda: "../wrappers/cnvkit/env.yaml"
        script: "../wrappers/cnvkit/cnvkit_call.py"
        # shell:
        #     "(cnvkit.py call -y -m clonal {input.segment} -v {input.vcf} -o {output.calls} --purity {params.TC} {params.extra}) &> {log}"

########################################################################################################################
rule cnvkit_export:
    input:
        segment="variant_calls/{sample}/cnvkit/CNV_calls.cns",
    output:
        vcf="variant_calls/{sample}/cnvkit/CNV_calls.vcf",
        bed="variant_calls/{sample}/cnvkit/CNV_calls.bed",
    params:
        sample="{sample}",
        hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
    log:
        "logs/{sample}/convert_to_vcf.log",
    threads: 6
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        cnvkit.py export bed {input.segment} --show all -o {output.bed}
        cnvkit.py export vcf {input.segment} -i {params.sample} -o {output.vcf}
        """
    # script:
    #     "../wrappers/cnvkit/cnvkit_vcf.py"


rule cnvkit_2tsv:
    input:
        bed="variant_calls/{sample}/cnvkit/CNV_calls.bed",
    output:
        tsv="variant_calls/{sample}/cnvkit/cnvkit_CNV_{sample}.tsv"
    threads: 3
    shell:
        "cp {input.bed} {output.tsv}"


########################################################################################################################
# Identify targeted genes with copy number gain or loss above or below a threshold.
rule cnvkit_genemetrics:
    input:
        cns="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
    output:
        trusted_genes="variant_calls/{sample}/cnvkit/trusted-genes.txt",
        ratio_genes="variant_calls/{sample}/cnvkit/ratio-genes.txt",
        segment_genes="variant_calls/{sample}/cnvkit/segment-genes.txt"
    params:
        sample="{sample}",
        threshold=0.2,
        min_probes=3,
        working_directory="variant_calls/{sample}/cnvkit/"
    log:
        "logs/{sample}/cnvkit_genemetrics.log"
    threads: 3
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        """
        (cnvkit.py genemetrics -y {input.cnr} -s {input.cns} -t {params.threshold} -m {params.min_probes} | tail -n+2 | cut -f1 | sort > {output.segment_genes}) &> {log}
        (cnvkit.py genemetrics -y {input.cnr} | tail -n+2 | cut -f1 | sort > {output.ratio_genes}) &> {log}
        comm -12 {output.ratio_genes} {output.segment_genes} > {output.trusted_genes}
        for gene in `cat trusted-genes.txt`
        do
            cnvkit.py scatter -s {params.sample}.cn{{s,r}} -g $gene -o Sample-$gene-scatter.pdf
        done
        """



########################################################################################################################
# PLOTS

rule cnvkit_scatter:
    input:
        cns="variant_calls/{sample}/cnvkit/CNV_calls.cns",
        cnr="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
    output:
        plot="variant_calls/{sample}/cnvkit/cnvkit_scatter.png",
    params:
        extra=" --fig-size 18 6 --trend",
        sample="{sample}"
    log:
        "logs/{sample}/cnvkit_scatter.log"
    threads: 6
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py scatter -s {input.cns} {input.cnr} -o {output.plot} --title {params.sample} {params.extra}) &> {log}"




def cnvkit_heatmap_input(wildcards):
    return expand("variant_calls/{sample}/cnvkit/CNV_calls.cns", sample=SAMPLES_LIST)


rule cnvkit_heatmap:
    input:
        cns=cnvkit_heatmap_input
    output:
        plot="variant_calls/all_samples/cnvkit_heatmap.png",
    params:
        extra="-d",
    log:
        "logs/all_samples/cnvkit_heatmap.log"
    threads: 6
    conda:
        "../wrappers/cnvkit/env.yaml"
    shell:
        "(cnvkit.py heatmap {input.cns} -o {output.plot} {params.extra}) &> {log}"

# rule cnvkit_diagram:
#     input:
#         cns="variant_calls/{sample}/cnvkit/CNV_calls.cns",
#         cnr="variant_calls/{sample}/cnvkit/fixed_cov.cnr",
#     output:
#         pdf="variant_calls/{sample}/cnvkit/cnvkit_diagram.pdf",
#     params:
#         extra="",
#     log:
#         "logs/{sample}/cnvkit_diagram.log"
#     threads: 10
#     resources: mem=10
#     conda:
#         "../wrappers/cnvkit/env.yaml"
#     shell:
#         "(cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} {params.extra}) &> {log}"






###########################################################################################################################################
# # new function defined elsewhere
# def caller_bam_input(wildcards):
#     if config["analysis_mode"] == "somatic_paired":
#         return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="tumor", "sample_name"])[1]
#     elif config["analysis_mode"] == "reference_cohort":
#         return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="normal", "sample_name"])[1]
#     elif config["analysis_mode"] == "baseline":
#         return expand("mapped/{all_bams}.bam", all_bams=SAMPLES_LIST)[1]
#     else:
#         return expand("mapped/{all_bams}.bam", all_bams=SAMPLES_LIST)[1]

# # old function
# def get_bam_input(wildcards):
#     if config["tumor_normal_paired"] == True:
#         return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample), "sample"])[0],
#     else:
#         return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample)[0]