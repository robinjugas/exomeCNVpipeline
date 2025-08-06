import re
import os
import pandas as pd
import json
import yaml

#testhh
# configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]
# GLOBAL_REF_PATH = "/home/rj/4TB/CEITEC/"

##### Config processing #####
#conversion from json
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


#### Reference info processing



#### FOLDERS



# ####################################
# # VARIABLES FROM CONFIG
used_SV_callers = []

if config["use_control_cnMOPS"]:
    used_SV_callers.append("cnMOPS")
if config["use_control_exomeDepth"]:
    used_SV_callers.append("exomeDepth")
if config["use_control_panelcnMOPS"]:
    used_SV_callers.append("panelcnMOPS")
    
# if config["use_cnvkit"]:
#     used_SV_callers.append("cnvkit")

# if config["use_gatk_cnv"]:
#     used_SV_callers.append("gatk_cnv")
# if config["use_control_freec"]:
#     used_SV_callers.append("control_freec")


wildcard_constraints:
    status = "tumor|normal|control"


####################################
# SEPARATE RULES
include: "rules/cnvPrepareCohort.smk"
include: "rules/cnvExomeCallers.smk"
include: "rules/cnvkit.smk"
include: "rules/cnvMergeWhole.smk"
# include: "rules/cnvMergeTargetRegions.smk"

####################################
# RULE ALL
# SAMPLES_LIST = sample_tab.loc[sample_tab["status"]=="tumor"),"sample_name"].to_string(index=False)

if config["analysis_mode"] == "somatic_paired":
    SAMPLES_LIST = sample_tab.loc[sample_tab["status"]=="tumor", "sample_name"]
elif config["analysis_mode"] == "reference_cohort":
    SAMPLES_LIST = sample_tab.loc[sample_tab["status"]=="normal", "sample_name"] 
elif config["analysis_mode"] == "baseline":
    SAMPLES_LIST = sample_tab.sample_name
else:
    SAMPLES_LIST = sample_tab.sample_name

# print("SAMPLES OF INTEREST \n")
# print(SAMPLES_LIST)



######
# for all samples:
ALL_SAMPLES = sample_tab.sample_name
STATUS = sample_tab.status


# print("All samples \n")
# print(ALL_SAMPLES)
# print("STATUS of all samples \n")
# print(STATUS)

sample_tab.to_csv('sample_tab.tsv', sep="\t")

rule all:
    input:
        # TargetRegions_merged=expand("CNV_TargetRegions/{sample}.callers_merged.tsv", sample=SAMPLES_LIST),
        # TargetRegions_annotated=expand("CNV_TargetRegions/{sample}_called_TargetRegions.tsv", sample=SAMPLES_LIST),
        merged_tsv=expand("mergedCNVs_final/{sample}_final_CNVs_annotated.tsv",sample=SAMPLES_LIST),
        merged_xlsx=expand("mergedCNVs_final/{sample}_final_CNVs_annotated.xlsx",sample=SAMPLES_LIST),

        #cnvkit
        tsv=expand("cnvkit_final/{sample}_final_CNVs_annotated.tsv",sample=ALL_SAMPLES,status=STATUS),
        xlsx=expand("cnvkit_final/{sample}_final_CNVs_annotated.xlsx",sample=ALL_SAMPLES,status=STATUS),

        targetcoverage=expand("variant_calls/{sample}/cnvkit/{status}.targetcoverage.cnn", sample=ALL_SAMPLES,status=STATUS),
        antitargetcoverage=expand("variant_calls/{sample}/cnvkit/{status}.antitargetcoverage.cnn", sample=ALL_SAMPLES,status=STATUS),
        calls=expand("variant_calls/{sample}/cnvkit/CNV_calls.cns", sample=SAMPLES_LIST),
        TSVs=expand("variant_calls/{sample}/cnvkit/cnvkit_CNV_{sample}.tsv", sample=SAMPLES_LIST),
        GENES=expand("variant_calls/{sample}/cnvkit/trusted-genes.txt", sample=SAMPLES_LIST),
        SCATTER=expand("variant_calls/{sample}/cnvkit/{sample}_cnvkit_scatter.png", sample=SAMPLES_LIST),
        heatmap="variant_calls/all_samples/cnvkit_heatmap.png"