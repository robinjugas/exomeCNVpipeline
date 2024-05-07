import os
import pandas as pd


#testhh
# configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]
# GLOBAL_REF_PATH = "/home/rj/4TB/CEITEC/"

##### Config processing #####
#conversion from json
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


# for CNV pipeline to get corresponding germinal sample column
# for index, row in sample_tab.iterrows():
#     donorvalue = sample_tab.loc[index,"donor"]
#     sample_tab.loc[index,"germinal"] = sample_tab.loc[(sample_tab["donor"] == donorvalue) & (sample_tab["tumor_normal"]=="normal"),"sample_name"].to_string(index=False)

#### Reference info processing

#### Setting reference from lib_ROI
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]

#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"))
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]




#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])


# ####################################
# # VARIABLES FROM CONFIG
used_SV_callers = []


if config["use_control_cnMOPS"]:
    used_SV_callers.append("cnMOPS")
if config["use_control_exomeDepth"]:
    used_SV_callers.append("exomeDepth")
if config["use_control_panelcnMOPS"]:
    used_SV_callers.append("panelcnMOPS")
# if config["use_gatk_cnv"]:
#     used_SV_callers.append("gatk_cnv")
# if config["use_control_freec"]:
#     used_SV_callers.append("control_freec")
# if config["use_cnvkit"]:
#     used_SV_callers.append("cnvkit")


wildcard_constraints:
    tumor_normal = "tumor|normal"


####################################
# SEPARATE RULES
include: "rules/cnvExomeprepareCohort.smk"
include: "rules/cnvExomeCallers.smk"
include: "rules/cnvExomeMerge.smk"
include: "rules/cnvExomeExonMerge.smk"



####################################
# RULE ALL
# SAMPLES_LIST = sample_tab.loc[sample_tab["tumor_normal"]=="tumor"),"sample_name"].to_string(index=False)

if config["tumor_normal_paired"] == True:
    SAMPLES_LIST = sample_tab.loc[sample_tab["tumor_normal"]=="tumor", "sample_name"]
else:
    SAMPLES_LIST = sample_tab.sample_name

print(SAMPLES_LIST)


rule all:
    input:
        # merged=expand("CNV_exon_merged/{sample}.callers_merged.tsv", sample=sample_tab.sample_name),
        Exons=expand("CNV_exon_merged/{sample}_final.tsv", sample=SAMPLES_LIST),
        CNVs=expand("CNV_merged/{sample}.callers_merged.tsv", sample=SAMPLES_LIST)


