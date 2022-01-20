import os
import re
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

# if not "ref_from_trans_assembly" in config:
#     config["ref_from_trans_assembly"] = "F"

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

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


# setting references for enrichment
# GO
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","GO_reference.json"),)
reference_GO = json.load(f)
f.close()
if config["species"] in reference_GO.keys():
    config["organism_go"] = reference_GO[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in GO references!")
    config["onthology"] = False

# KEGG
ff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","kegg_reference.json"),)
reference_kegg = json.load(ff)
ff.close()
if config["species"] in reference_kegg.keys():
    config["organism_kegg"] = reference_kegg[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in KEGG references!")
    config["kegg"] = False

# WIKIPATHWAYS
fff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","wp_reference.json"),)
reference_wp = json.load(fff)
fff.close()
if config["species"] in reference_wp.keys():
    config["organism_wp"] = reference_wp[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in WikiPathways references!")
    config["wikipathways"] = False

# REACTOME
ffff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reactome_reference.json"),)
reference_reactome = json.load(ffff)
ffff.close()
if config["species"] in reference_reactome.keys():
    config["organism_reactome"] = reference_reactome[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in REACTOME references!")
    config["reactome"] = False

#
wildcard_constraints:
    sample = "|".join(sample_tab.sample_name) + "|all_samples",
    analysis_type = "feature_count|RSEM",

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
        if config["onthology"]:
            input["feature_count_go_e"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/enrichment_GO/GO_enrich_CC.png", comparison=comparison_dir_list, biotype=biotype_dir_list),
            input["feature_count_go_g"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_CC.png",comparison=comparison_dir_list,biotype=biotype_dir_list)
        if config["kegg"]:
            input["feature_count_kegg_e"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/enrichment_KEGG/KEGG_enrich.png", comparison=comparison_dir_list, biotype=biotype_dir_list),
            input["feature_count_kegg_g"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/GSEA_KEGG/GSEA_KEGG.png", comparison=comparison_dir_list, biotype=biotype_dir_list)
        if config["wikipathways"]:
            input["feature_count_wp_e"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/enrichment_WP/WP_enrich.png", comparison=comparison_dir_list, biotype=biotype_dir_list),
            input["feature_count_wp_g"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/GSEA_WP/GSEA_WP.png", comparison=comparison_dir_list, biotype=biotype_dir_list)
        if config["reactome"]:
            input["feature_count_reactome_e"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/enrichment_REACTOME/REACTOME_enrich.png", comparison=comparison_dir_list, biotype=biotype_dir_list),
            input["feature_count_reactome_g"] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/GSEA_REACTOME/GSEA_REACTOME.png", comparison=comparison_dir_list, biotype=biotype_dir_list)

    if config["RSEM"]:
        if config["onthology"]:
            input["RSEM_go"] = [
                "results/DE_RSEM/{comparison}/{biotype}/enrichment_GO/GO_enrich_CC.png",
                "results/DE_RSEM/{comparison}/{biotype}/GSEA_GO/GSEA_GO_CC.png"]
        if config["kegg"]:
            input["RSEM_kegg"] = [
                "results/DE_RSEM/{comparison}/{biotype}/enrichment_KEGG/KEGG_enrich.png",
                "results/DE_RSEM/{comparison}/{biotype}/GSEA_KEGG/GSEA_KEGG.png"]
        if config["wikipathways"]:
            input["RSEM_wp"] = ["results/DE_RSEM/{comparison}/{biotype}/enrichment_WP/WP_enrich.png",
                "results/DE_RSEM/{comparison}/{biotype}/GSEA_WP/GSEA_WP.png"]
        if config["reactome"]:
            input["RSEM_reactome"] = [
                "results/DE_RSEM/{comparison}/{biotype}/enrichment_REACTOME/REACTOME_enrich.png",
                "results/DE_RSEM/{comparison}/{biotype}/GSEA_REACTOME/GSEA_REACTOME.png"]
    return input

rule all:
    input:  unpack(input_all)

##### Modules #####

include: "rules/enrichment_analysis.smk"