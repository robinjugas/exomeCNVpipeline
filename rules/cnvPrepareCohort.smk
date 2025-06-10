def cohort_bam_inputs(wildcards):
    if config["analysis_mode"] == "somatic_paired":
        return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="normal", "sample_name"])
    elif config["analysis_mode"] == "reference_cohort":
        return expand("mapped/{normal_bams}.bam", normal_bams=sample_tab.loc[sample_tab["status"]=="control", "sample_name"])
    elif config["analysis_mode"] == "baseline":
        return expand("mapped/{all_bams}.bam", all_bams=SAMPLES_LIST)
    else:
        return expand("mapped/{all_bams}.bam", all_bams=SAMPLES_LIST)


# sample_tab.loc[index,"germinal"] = sample_tab.loc[(sample_tab["donor"] == donorvalue) & (sample_tab["status"]=="normal"),"sample_name"].to_string(index=False)

rule prepare_cnMOPS:
    input:  bams = cohort_bam_inputs,
            bed=config["bed_file"], #expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            ref=config["genome"] #expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output: Rdata="cohort_data/cnMOPS_customCohort.Rdata"
    log: "logs/prepare_cnMOPS.log"
    threads: 8
    resources: mem=64
    conda: "../wrappers/prepare_cnMOPS/env.yaml"
    script: "../wrappers/prepare_cnMOPS/script.py"


rule prepare_exomeDepth:
    input:  bams = cohort_bam_inputs,
            bed=config["bed_file"], #expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            ref=config["genome"] #expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output: Rdata="cohort_data/exomeDepth_customCohort.Rdata"
    log: "logs/prepare_exomeDepth.log"
    threads: 5
    resources: mem=64
    conda: "../wrappers/prepare_exomeDepth/env.yaml"
    script: "../wrappers/prepare_exomeDepth/script.py"



rule prepare_panelcnMOPS:
    input:  bams = cohort_bam_inputs,
            bed=config["bed_file"], 
            ref=config["genome"]
    output: Rdata="cohort_data/panelcnMOPS_customCohort.Rdata"
    log: "logs/prepare_panelcnMOPS.log"
    threads: 5
    resources: mem=64
    conda: "../wrappers/prepare_panelcnMOPS/env.yaml"
    script: "../wrappers/prepare_panelcnMOPS/script.py"