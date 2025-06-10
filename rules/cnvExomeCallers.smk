
def caller_bam_input(wildcards):
    if config["analysis_mode"] == "somatic_paired":
        return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="tumor", "sample_name"])[1]
    elif config["analysis_mode"] == "reference_cohort":
        return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="normal", "sample_name"])[1]
    elif config["analysis_mode"] == "baseline":
        return expand("mapped/{all_bams}.bam", all_bams=SAMPLES_LIST)[1]
    else:
        return expand("mapped/{all_bams}.bam", all_bams=SAMPLES_LIST)[1]


rule cnMOPS:
    input:  bam = caller_bam_input,
            ref = config["genome"],
            cohort = "cohort_data/cnMOPS_customCohort.Rdata"
    output: tsv = "variant_calls/{sample}/cnMOPS/cnMOPS_CNV_{sample}.tsv"
    log: "logs/{sample}/{sample}_cnMOPS.log"
    threads: 5
    resources: mem=8
    params:
        sampleName="{sample}"
    conda: "../wrappers/run_cnMOPS/env.yaml"
    script: "../wrappers/run_cnMOPS/script.py"

rule exomeDepth:
    input:  bam = caller_bam_input,
            ref = config["genome"],
            cohort = "cohort_data/exomeDepth_customCohort.Rdata"
    output: tsv = "variant_calls/{sample}/exomeDepth/exomeDepth_CNV_{sample}.tsv"
    log: "logs/{sample}/{sample}_exomeDepth.log"
    threads: 5
    resources: mem=8
    params:
        sampleName="{sample}"
    conda: "../wrappers/run_exomeDepth/env.yaml"
    script: "../wrappers/run_exomeDepth/script.py"

rule panelcnMOPS:
    input:  bam = caller_bam_input,
            ref = config["genome"],
            cohort = "cohort_data/panelcnMOPS_customCohort.Rdata"
    output: tsv = "variant_calls/{sample}/panelcnMOPS/panelcnMOPS_CNV_{sample}.tsv"
    log: "logs/{sample}/{sample}_panelcnMOPS.log"
    threads: 5
    resources: mem=8
    params:
        sampleName="{sample}"
    conda: "../wrappers/run_panelcnMOPS/env.yaml"
    script: "../wrappers/run_panelcnMOPS/script.py"