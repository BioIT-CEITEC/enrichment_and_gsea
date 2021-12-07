def DE_computation_input(wildcards):
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

        input['tsv'] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv", comparison=comparison_dir_list, biotype=biotype_dir_list)

        if config["ref_from_trans_assembly"] != False:
            input['trans_ids_map'] = expand("{ref_dir}/annot/{ref}.transdecoder_ids_map", ref_dir=reference_directory,ref=config["reference"])[0]
    else:
        raise ValueError("There is no conditions or tag for samples!")
    return input


rule enrichment_GO:
    input: unpack(DE_computation_input)
    output: table = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_CC.png",
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir= "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO",
            organism_go = config["organism_go"],
            cutoff_log2fc = config["cutoff_log2fc"],
            cutoff_padj = config["cutoff_padj"],
            n_up = config["n_up"],
            colors = config["colors"],
            enrich_padj = config["enrich_padj"],
            enrich_padjmethod = config["enrich_padjmethod"],
            enrich_minGSSize = config["enrich_minGSSize"],
            enrich_maxGSSize = config["enrich_maxGSSize"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.enrichment_GO.log"
    conda:  "../wrappers/enrichment_GO/env.yaml"
    script: "../wrappers/enrichment_GO/script.py"

rule GSEA_GO:
    input: unpack(DE_computation_input)
    output: table = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_CC.png",
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO",
            organism_go = config["organism_go"],
            cutoff_log2fc = config["cutoff_log2fc"],
            cutoff_padj = config["cutoff_padj"],
            n_up = config["n_up"],
            n_down= config["n_down"],
            colors = config["colors"],
            gsea_padj = config["gsea_padj"],
            gsea_padjmethod = config["gsea_padjmethod"],
            gsea_minGSSize = config["gsea_minGSSize"],
            gsea_maxGSSize = config["gsea_maxGSSize"],
            gsea_eps = config["gsea_eps"],
            gsea_nPermSimple = config["gsea_nPermSimple"],
            gsea_by = config["gsea_by"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.GSEA_GO.log"
    conda:  "../wrappers/GSEA_GO/env.yaml"
    script: "../wrappers/GSEA_GO/script.py"
