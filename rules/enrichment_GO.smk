def DE_computation_input(wildcards):
    input={}
    input["cfg_tab"] = "config.json"
    input["biotype_groups"] = os.path.join(GLOBAL_REF_PATH,"general/default/annot/biotypes_list_mod.txt")
    input["sqlite"] = expand("{ref_dir}/annot/{ref}.sqlite.gz",ref_dir=reference_directory,ref=config["reference"])[0]
    if wildcards.analysis_type == "feature_count":
        input["expression_tab"] = "results/analysis_{analysis_type}_table/complete.{analysis_type}.tsv"
    if wildcards.analysis_type == "RSEM":
        input["expression_tab"] = "results/analysis_{analysis_type}_table/complete.{analysis_type}.RData"
    return input


rule enrichment_GO:
    input: unpack(DE_computation_input)
    output: table = "results/DE_{analysis_type}/{comparison}/{biotype}/edgeR.tsv",
    params: count_type="{analysis_type}",
            organism = config["organism"],
            use_tag_to_pair_samples = config["use_tag_to_pair_samples"],
            ref_from_trans_assembly = config["ref_from_trans_assembly"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.log"
    conda:  "../wrappers/DE_computation/env.yaml"
    script: "../wrappers/DE_computation/script.py"
