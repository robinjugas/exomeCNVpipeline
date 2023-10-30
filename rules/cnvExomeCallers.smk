rule cnMOPS:
    input:  bam = "mapped/{sample}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0]
    output: tsv="variant_calls/{sample}/cnMOPS/cnMOPS_CNV_{sample}.tsv"
    log: "logs/{sample}/cnMOPS/{sample}_cnMOPS.log"
    threads: 5
    resources: mem=8
    params:
        libROI=config["lib_ROI"],
        sampleName="{sample}"
    conda: "../wrappers/cnMOPS/env.yaml"
    script: "../wrappers/cnMOPS/script.py"

rule exomeDepth:
    input:  bam = "mapped/{sample}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0]
    output: tsv="variant_calls/{sample}/exomeDepth/exomeDepth_CNV_{sample}.tsv"
    log: "logs/{sample}/exomeDepth/{sample}_exomeDepth.log"
    threads: 5
    resources: mem=8
    params:
        libROI=config["lib_ROI"],
        sampleName="{sample}"
    conda: "../wrappers/exomeDepth/env.yaml"
    script: "../wrappers/exomeDepth/script.py"



rule panelcnMOPS:
    input:  bam = "mapped/{sample}.bam",
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0]
    output: tsv="variant_calls/{sample}/panelcnMOPS/panelcnMOPS_CNV_{sample}.tsv"
    log: "logs/{sample}/panelcnMOPS/{sample}_panelcnMOPS.log"
    threads: 5
    resources: mem=8
    params:
        libROI=config["lib_ROI"],
        sampleName="{sample}"
    conda: "../wrappers/panelcnMOPS/env.yaml"
    script: "../wrappers/panelcnMOPS/script.py"