

# SURVIVOR MERGE
rule survivor_merge:
    input:
        vcfs=expand("results/{cnv_caller}/{{SAMPLE}}/{{SAMPLE}}.vcf",cnv_caller=SVcallers)
    output:
        vcf="results/merged_survivor/{SAMPLE}.survivor_merged.vcf",
    params:
        sample_files="results/merged_survivor/ls_{SAMPLE}.txt",
        max_allowed_space=config["survivor_max_allowed_space"],#1000,
        min_callers=config["survivor_minCallers"],
        agree_on_type=1,
        agree_on_strand=0,
        estimate_sv_distance=0,
        min_length=config["survivor_min_sv_length"],
    log:
        "logs/survivor_merge/{SAMPLE}.log",
    threads: 2
    conda:
        "../wrappers/survivor/env.yaml"
    script:
        "../wrappers/survivor/script.py"
