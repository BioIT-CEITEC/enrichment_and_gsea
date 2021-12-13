import os
import re
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

if not "ref_from_trans_assembly" in config:
    config["ref_from_trans_assembly"] = "F"

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
#config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["species"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = [re.sub(r" \(.*\)","",organism_name).lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# setting references for enrichment
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","GO_reference.json"),)
reference_GO = json.load(f)
f.close()
if config["species"] in reference_GO.keys():
    config["organism_go"] = reference_GO[config["species"]]
    config["onthology"] = True

else:
    raise ValueError("There is no "+config["species"]+" in GO references!")
    config["onthology"] = False


f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","kegg_reference.json"),)
reference_kegg = json.load(f)
f.close()
config["organism_kegg"] = reference_kegg[config["species"]]

f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","wp_reference.json"),)
reference_wp = json.load(f)
f.close()
config["organism_wp"] = reference_wp[config["species"]]

f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reactome_reference.json"),)
reference_wp = json.load(f)
f.close()
config["organism_reactome"] = reference_wp[config["species"]]


# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     analysis_type= "feature_count|RSEM",
     #data_type= "tsv|RData"


##### Target rules #####

def input_all(wildcards):
    input = {}
    if config["feature_count"]:
        input["feature_count"] = "results/feature_count_final_report.html"
    if config["RSEM"]:
        input["RSEM"] = "results/RSEM_final_report.html"
    return input

rule all:
    input:  unpack(input_all)

##### Modules #####

include: "rules/enrichment_analysis.smk"