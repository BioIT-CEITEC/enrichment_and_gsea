import os
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
config["organism"] = [species_name.split(" (")[0].lower().replace(" ","_") for species_name in config["species"]][0]

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

rule all:
    input:  html = "results/enrichment_GSEA_final_report.html"

##### Modules #####

include: "rules/enrichment_analysis.smk"