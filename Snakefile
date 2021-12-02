import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Conditions
#
def gene_enrichment_input(wildcards):
  input = {}
  condition_list = config["conditions"]
  comparison_dir_list = list()
  for condition1 in condition_list:
    if ':' in condition1:
      conditions = condition1.split(":")
      comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
    else:
      for condition2 in condition_list[condition_list.index(condition1):]:
        if ':' not in condition2 and condition2 != condition1:
          comparison_dir_list.append(condition2 + "_vs_" + condition1)
  input['tsv'] = expand("mRNA_DE_{de_analysis}/{comparison}/all/DESeq2.tsv",comparison=comparison_dir_list)



wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?",

##### Target rules #####
rule all:
    input:
        expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
        expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name),
        expand("mapped/{sample}.bigWig",sample = sample_tab.sample_name),
        expand("mapped/{sample}.bedGraph",sample = sample_tab.sample_name),
        "qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html"

##### Modules #####

include: "rules/enrichment_GO.smk"