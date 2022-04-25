import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

# if not "ref_from_trans_assembly" in config:
#     config["ref_from_trans_assembly"] = "F"

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"))
reference_dict = json.load(f)
f.close()

#config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["species"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species"].split(" (")[0].lower().replace(" ","_")
print(config["species"])
print(config["organism"])

##### Config processing #####
# Folders
#
#reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

# setting references for enrichment
# GO
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","GO_reference.json"))
reference_GO = json.load(f)
f.close()
if config["species"] in reference_GO.keys():
    config["organism_go"] = reference_GO[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in GO references!")
    config["onthology"] = False

# KEGG
ff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","kegg_reference.json"))
reference_kegg = json.load(ff)
ff.close()
if config["species"] in reference_kegg.keys():
    config["organism_kegg"] = reference_kegg[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in KEGG references!")
    config["kegg"] = False

# WIKIPATHWAYS
fff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","wp_reference.json"))
reference_wp = json.load(fff)
fff.close()
if config["species"] in reference_wp.keys():
    config["organism_wp"] = reference_wp[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in WikiPathways references!")
    config["wikipathways"] = False

# REACTOME
ffff = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reactome_reference.json"))
reference_reactome = json.load(ffff)
ffff.close()
if config["species"] in reference_reactome.keys():
    config["organism_reactome"] = reference_reactome[config["species"]]
else:
    raise ValueError("There is no "+config["species"]+" in REACTOME references!")
    config["reactome"] = False

#
analysis = []
if config["feature_count"]:
    analysis.append("feature_count")
if config["RSEM"]:
    analysis.append("RSEM")

condition_list = sorted(sample_tab.condition.unique()) if config['conditions_to_compare'] == "all" else config['conditions_to_compare'].split(",")
biotype_dir_list = config['biotypes'].split(",")
comparison_dir_list = list()
for condition1 in condition_list:
    if ':' in condition1:
        conditions = condition1.split(":")
        comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
    else:
        for condition2 in condition_list[condition_list.index(condition1):]:
            if ':' not in condition2 and condition2 != condition1:
                comparison_dir_list.append(condition2 + "_vs_" + condition1)

config["analysis_type"] = analysis
config["biotype_list"] = biotype_dir_list
config["comparison"] = comparison_dir_list

os.makedirs("enrichment_gsea",exist_ok=True)

f=open("enrichment_gsea/config_enrichment_gsea.json","w")
json.dump(config,f,indent = 4)
f.close()

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name) + "|all_samples",
    lib_name = "[^\.\/]+",
    analysis_type = "|".join(analysis),
    condition_list = "|".join(condition_list),
    biotype = "|".join(biotype_dir_list),
    comparison = "|".join(comparison_dir_list)

##### Target rules #####

rule all:
    input:  report = "enrichment_gsea/enrichment_GSEA_final_report.html"

##### Modules #####

include: "rules/enrichment_analysis.smk"
