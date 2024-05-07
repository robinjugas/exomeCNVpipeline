
rule cnMOPS:
    input:  bam = "mapped/{sample}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            cohort="cohort_data/cnMOPS_customCohort.Rdata"
    output: tsv="variant_calls/{sample}/cnMOPS/cnMOPS_CNV_{sample}.tsv"
    log: "logs/{sample}/{sample}_cnMOPS.log"
    threads: 5
    resources: mem=8
    params:
        libROI=config["lib_ROI"],
        sampleName="{sample}"
    conda: "../wrappers/run_cnMOPS/env.yaml"
    script: "../wrappers/run_cnMOPS/script.py"

rule exomeDepth:
    input:  bam = "mapped/{sample}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            cohort="cohort_data/exomeDepth_customCohort.Rdata"
    output: tsv="variant_calls/{sample}/exomeDepth/exomeDepth_CNV_{sample}.tsv"
    log: "logs/{sample}/{sample}_exomeDepth.log"
    threads: 5
    resources: mem=8
    params:
        libROI=config["lib_ROI"],
        sampleName="{sample}"
    conda: "../wrappers/run_exomeDepth/env.yaml"
    script: "../wrappers/run_exomeDepth/script.py"



rule panelcnMOPS:
    input:  bam = "mapped/{sample}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            cohort="cohort_data/panelcnMOPS_customCohort.Rdata"
    output: tsv="variant_calls/{sample}/panelcnMOPS/panelcnMOPS_CNV_{sample}.tsv"
    log: "logs/{sample}/{sample}_panelcnMOPS.log"
    threads: 5
    resources: mem=8
    params:
        libROI=config["lib_ROI"],
        sampleName="{sample}"
    conda: "../wrappers/run_panelcnMOPS/env.yaml"
    script: "../wrappers/run_panelcnMOPS/script.py"