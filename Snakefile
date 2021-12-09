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
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
#config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["species"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = [re.sub(r" \(.*\)","",organism_name).lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

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

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

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
    if (sample_tab.condition != "").all() and (sample_tab.tag != "").all():
        if config['conditions_to_compare'] == "all":
            condition_list = sorted(sample_tab.condition)
        else:
            condition_list = config['conditions_to_compare'].split(",")

        comparison_dir_list = list()
        for condition1 in condition_list:
            if ':' in condition1:
                conditions = condition1.split(":")
                comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
            else:
                for condition2 in condition_list[condition_list.index(condition1):]:
                    if ':' not in condition2 and condition2 != condition1:
                        comparison_dir_list.append(condition2 + "_vs_" + condition1)

        biotype_dir_list = config['biotypes'].split(",")

        if config["feature_count"]:
            input['tsv'] = expand("results/DE_feature_count/{comparison}/{biotype}/DESeq2.tsv",comparison=comparison_dir_list,biotype=biotype_dir_list)
        if config["RSEM"]:
            input['tsv'] = expand("results/DE_RSEM/{comparison}/{biotype}/DESeq2.tsv",comparison=comparison_dir_list,biotype=biotype_dir_list)

        # input['tsv'] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv",comparison=comparison_dir_list,biotype=biotype_dir_list)

        if config["ref_from_trans_assembly"] != False:
            input['trans_ids_map'] = \
            expand("{ref_dir}/annot/{ref}.transdecoder_ids_map",ref_dir=reference_directory,ref=config["reference"])[0]
    else:
        raise ValueError("There is no conditions or tag for samples!")
    return input

rule all:
    input:  unpack(input_all)


##### Modules #####

include: "rules/enrichment_analysis.smk"