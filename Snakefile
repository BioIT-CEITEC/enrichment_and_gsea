import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

# config = BR.load_organism()

# setting references for enrichment
# GO
# f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","GO_reference.json"))
# reference_GO = json.load(f)
# f.close()
# if config["onthology"]:
#     if config["species"] in reference_GO.keys():
#         config["organism_go"] = reference_GO[config["species"]]
#     else:
#         raise ValueError("There is no "+config["species"]+" in GO references!")
#         config["onthology"] = False
#         config["organism_go"] = ""
# else:
#     config["organism_go"] = ""
#
# # KEGG
# ff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","kegg_reference.json"))
# reference_kegg = json.load(ff)
# ff.close()
# if config["kegg"]:
#     if config["species"] in reference_kegg.keys():
#         config["organism_kegg"] = reference_kegg[config["species"]]
#     else:
#         raise ValueError("There is no "+config["species"]+" in KEGG references!")
#         config["kegg"] = False
#         config["organism_kegg"] = ""
# else:
#     config["organism_kegg"] = ""
#
#     # WIKIPATHWAYS
# fff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","wp_reference.json"))
# reference_wp = json.load(fff)
# fff.close()
# if config["wikipathways"]:
#     if config["species"] in reference_wp.keys():
#         config["organism_wp"] = reference_wp[config["species"]]
#     else:
#         raise ValueError("There is no "+config["species"]+" in WikiPathways references!")
#         config["wikipathways"] = False
#         config["organism_wp"] = ""
# else:
#     config["organism_wp"] = ""
#
#     # REACTOME
# ffff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reactome_reference.json"))
# reference_reactome = json.load(ffff)
# ffff.close()
# if config["reactome"]:
#     if config["species"] in reference_reactome.keys():
#         config["organism_reactome"] = reference_reactome[config["species"]]
#     else:
#         raise ValueError("There is no "+config["species"]+" in REACTOME references!")
#         config["reactome"] = False
#         config["organism_reactome"] = ""
# else:
#     config["organism_reactome"] = ""

#set analysis selected analysis types from config and rise exception if no selected
analysis = []
if config["featureCount"]:
    count_over_list = config['count_over'].split(",")
    if ("exon" in count_over_list):
        config["featureCount_exon"] = True
        analysis.append("featureCount_exon")
    if ("gene" in count_over_list):
        config["featureCount_gene"] = True
        analysis.append("featureCount_gene")
    if ("transcript" in count_over_list):
        config["featureCount_transcript"] = True
        analysis.append("featureCount_transcript")
    if ("three_prime_UTR" in count_over_list):
        config["featureCount_3pUTR"] = True
        analysis.append("featureCount_3pUTR")
    if ("five_prime_UTR" in count_over_list):
        config["featureCount_5pUTR"] = True
        analysis.append("featureCount_5pUTR")

if config["RSEM"]:
    analysis.append("RSEM")
if config["salmon_align"]:
    analysis.append("salmon_align")
if config["salmon_map"]:
    analysis.append("salmon_map")
if config["kallisto"]:
    analysis.append("kallisto")
if len(analysis) == 0:
    raise ValueError("There was no RSEM or featureCount used in previous analysis!")

if config["kegg"] == True and config["organism_kegg"] != "-":
    config["kegg_path"] = os.path.join(GLOBAL_REF_PATH, "tools","KEGG", config["organism_kegg"])

def get_comparison_dir_list(condition_list):
    comparison_dir_list = list()
    for condition1 in condition_list:
        if ':' in condition1:
            conditions = condition1.split(":")
            comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
        else:
            for condition2 in condition_list[condition_list.index(condition1):]:
                if ':' not in condition2 and condition2 != condition1:
                    comparison_dir_list.append(condition1 + "_vs_" + condition2)
    return comparison_dir_list


## create list of conditions
if config['conditions_to_compare'] == "all":
    condition_list = sorted(sample_tab.condition.unique())
    condition_list_first = [condition for condition in condition_list if
                            not re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list_second = [condition for condition in condition_list if
                             re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list = condition_list_first + condition_list_second
    comparison_dir_list = get_comparison_dir_list(condition_list)
else:
    comparison_dir_list = get_comparison_dir_list(config['conditions_to_compare'].split(","))
    condition_list = set(config['conditions_to_compare'].replace(':',',').split(","))
    # if not config['keep_not_compared_samples_for_normalization']:
    #     sample_tab = sample_tab[sample_tab['condition'].isin(condition_list)]

biotype_dir_list = config['biotypes'].split(",")

config["analysis_type"] = "|".join(analysis)
config["biotype_list"] = "|".join(biotype_dir_list)
config["comparison"] = "|".join(comparison_dir_list)

os.makedirs("enrichment_gsea",exist_ok=True)

f=open("enrichment_gsea/config_enrichment_gsea.json","w")
json.dump(config,f,indent = 4)
f.close()

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name) + "|all_samples",
    lib_name = "[^\.\/]+",
    analysis_type= "featureCount_exon|featureCount_gene|featureCount_transcript|featureCount_3pUTRn|featureCount_5pUTR|RSEM|salmon_map|salmon_align|kallisto",
    condition_list = "|".join(condition_list),
    biotype = "|".join(biotype_dir_list),
    comparison = "|".join(comparison_dir_list)

##### Target rules #####

rule all:
    input:  report = "enrichment_gsea/enrichment_GSEA_final_report.html"

##### Modules #####

include: "rules/enrichment_analysis.smk"
