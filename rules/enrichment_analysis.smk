def final_input(wildcards):
    input = {}
    if config["onthology"]:
        input["gobpE"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_BP.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
        input["gomfE"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_MF.svg", analysis_type=analysis,comparison=comparison_dir_list,biotype=biotype_dir_list)
        input["goccE"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_CC.svg", analysis_type=analysis,comparison=comparison_dir_list,biotype=biotype_dir_list)
        input["gobpG"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_BP.svg", analysis_type=analysis,comparison=comparison_dir_list,biotype=biotype_dir_list)
        input["gomfG"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_MF.svg", analysis_type=analysis,comparison=comparison_dir_list,biotype=biotype_dir_list)
        input["goccG"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_CC.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
    if config["kegg"]:
        input["keggE"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_KEGG/KEGG_enrich.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
        input["keggG"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_KEGG/GSEA_KEGG.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
    if config["reactome"]:
        input["reactE"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_REACTOME/REACTOME_enrich.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
        input["reactG"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_REACTOME/GSEA_REACTOME.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
    if config["wikipathways"]:
        input["wpE"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_WP/WP_enrich.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)
        input["wpG"] = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_WP/GSEA_WP.svg", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list)

    return input


rule final_report:
    input:  txtfile = "enrichment_gsea/config_enrichment_gsea.txt"
    output: html = "enrichment_gsea/enrichment_GSEA_final_report.html"
    #params: config = "enrichment_gsea/config_enrichment_gsea.json"
    conda: "../wrappers/final_report/env.yaml"
    script: "../wrappers/final_report/enrichment_GSEA_final_report.Rmd"
    #shell: "touch {output.html}"

rule completion:
    input:  unpack(final_input),
            enrich = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/gene_for_enrichment.tsv", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list),
            gsea = expand("enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/gene_for_gsea.tsv", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list),
    output: txtfile = "enrichment_gsea/config_enrichment_gsea.txt"
    params: config = "./config.json"
    #conda: "../wrappers/final_report/env.yaml"
    #script: "../wrappers/final_report/enrichment_GSEA_final_report.Rmd"
    shell: "touch {output.txtfile}"

rule sampling:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: enrich = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/gene_for_enrichment.tsv",
            gsea = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/gene_for_gsea.tsv",
    params: universe = "enrichment_gsea/gene_universe.tsv",
            organism_go = config["organism_go"],
            cutoff_log2fc_enrich = config["cutoff_log2fc_enrich"],
            cutoff_padj_enrich = config["cutoff_padj_enrich"],
            cutoff_log2fc_gsea = config["cutoff_log2fc_gsea"],
            cutoff_padj_gsea =config["cutoff_padj_gsea"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.sampling.log"
    conda:  "../wrappers/sampling/env.yaml"
    script: "../wrappers/sampling/script.py"

rule enrichment_GO:
    input:  tsv = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/gene_for_enrichment.tsv"
    output: plotBP = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_BP.svg",
            plotMF = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_MF.svg",
            plotCC = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_CC.svg"
    params: outdir= "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO",
            organism_go = config["organism_go"],
            n_up = config["n_up"],
            colors = config["colors"],
            enrich_padj = config["enrich_padj"],
            enrich_padjmethod = config["enrich_padjmethod"],
            enrich_minGSSize = config["enrich_minGSSize"],
            enrich_maxGSSize = config["enrich_maxGSSize"],
            universe = "enrichment_gsea/gene_universe.tsv"
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.enrichment_GO.log"
    conda:  "../wrappers/enrichment_GO/env.yaml"
    script: "../wrappers/enrichment_GO/script.py"

rule GSEA_GO:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plotBP = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_BP.svg",
            plotMF = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_MF.svg",
            plotCC = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_CC.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO",
            organism_go = config["organism_go"],
            cutoff_log2fc = config["cutoff_log2fc_gsea"],
            cutoff_padj = config["cutoff_padj_gsea"],
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

rule enrichment_kegg:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_KEGG/KEGG_enrich.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_KEGG",
            organism_kegg = config["organism_kegg"],
            cutoff_log2fc = config["cutoff_log2fc_enrich"],
            cutoff_padj = config["cutoff_padj_enrich"],
            n_up = config["n_up"],
            colors = config["colors"],
            enrich_padj = config["enrich_padj"],
            enrich_padjmethod = config["enrich_padjmethod"],
            enrich_minGSSize = config["enrich_minGSSize"],
            enrich_maxGSSize = config["enrich_maxGSSize"],
            organism_go = config["organism_go"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.enrichment_KEGG.log"
    conda:  "../wrappers/enrichment_kegg/env.yaml"
    script: "../wrappers/enrichment_kegg/script.py"

rule GSEA_kegg:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_KEGG/GSEA_KEGG.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_KEGG",
            organism_kegg = config["organism_kegg"],
            cutoff_log2fc = config["cutoff_log2fc_gsea"],
            cutoff_padj = config["cutoff_padj_gsea"],
            n_up = config["n_up"],
            n_down= config["n_down"],
            colors = config["colors"],
            gsea_padj = config["gsea_padj"],
            gsea_padjmethod = config["gsea_padjmethod"],
            gsea_minGSSize = config["gsea_minGSSize"],
            gsea_maxGSSize = config["gsea_maxGSSize"],
            gsea_eps = config["gsea_eps"],
            gsea_nPermSimple = config["gsea_nPermSimple"],
            gsea_by = config["gsea_by"],
            organism_go = config["organism_go"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.GSEA_KEGG.log"
    conda:  "../wrappers/GSEA_kegg/env.yaml"
    script: "../wrappers/GSEA_kegg/script.py"

rule enrichment_reactome:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_REACTOME/REACTOME_enrich.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_REACTOME",
            organism_reactome = config["organism_reactome"],
            cutoff_log2fc = config["cutoff_log2fc_enrich"],
            cutoff_padj = config["cutoff_padj_enrich"],
            n_up = config["n_up"],
            colors = config["colors"],
            enrich_padj = config["enrich_padj"],
            enrich_padjmethod = config["enrich_padjmethod"],
            enrich_minGSSize = config["enrich_minGSSize"],
            enrich_maxGSSize = config["enrich_maxGSSize"],
            organism_go = config["organism_go"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.enrichment_REACTOME.log"
    conda:  "../wrappers/enrichment_reactome/env.yaml"
    script: "../wrappers/enrichment_reactome/script.py"

rule GSEA_reactome:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_REACTOME/GSEA_REACTOME.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_REACTOME",
            organism_reactome = config["organism_reactome"],
            cutoff_log2fc = config["cutoff_log2fc_gsea"],
            cutoff_padj = config["cutoff_padj_gsea"],
            n_up = config["n_up"],
            n_down= config["n_down"],
            colors = config["colors"],
            gsea_padj = config["gsea_padj"],
            gsea_padjmethod = config["gsea_padjmethod"],
            gsea_minGSSize = config["gsea_minGSSize"],
            gsea_maxGSSize = config["gsea_maxGSSize"],
            gsea_eps = config["gsea_eps"],
            gsea_nPermSimple = config["gsea_nPermSimple"],
            gsea_by = config["gsea_by"],
            organism_go = config["organism_go"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.GSEA_REACTOME.log"
    conda:  "../wrappers/GSEA_reactome/env.yaml"
    script: "../wrappers/GSEA_reactome/script.py"

rule enrichment_wp:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_WP/WP_enrich.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/enrichment_WP",
            organism_wp = config["organism_wp"],
            cutoff_log2fc = config["cutoff_log2fc_enrich"],
            cutoff_padj = config["cutoff_padj_enrich"],
            n_up = config["n_up"],
            colors = config["colors"],
            enrich_padj = config["enrich_padj"],
            enrich_padjmethod = config["enrich_padjmethod"],
            enrich_minGSSize = config["enrich_minGSSize"],
            enrich_maxGSSize = config["enrich_maxGSSize"],
            organism_go = config["organism_go"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.enrichment_WP.log"
    conda:  "../wrappers/enrichment_wp/env.yaml"
    script: "../wrappers/enrichment_wp/script.py"

rule GSEA_wp:
    input:  tsv = "results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_WP/GSEA_WP.svg"
    params: workdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "enrichment_gsea/DE_{analysis_type}/{comparison}/{biotype}/GSEA_WP",
            organism_wp = config["organism_wp"],
            cutoff_log2fc = config["cutoff_log2fc_gsea"],
            cutoff_padj = config["cutoff_padj_gsea"],
            n_up = config["n_up"],
            n_down= config["n_down"],
            colors = config["colors"],
            gsea_padj = config["gsea_padj"],
            gsea_padjmethod = config["gsea_padjmethod"],
            gsea_minGSSize = config["gsea_minGSSize"],
            gsea_maxGSSize = config["gsea_maxGSSize"],
            gsea_eps = config["gsea_eps"],
            gsea_nPermSimple = config["gsea_nPermSimple"],
            gsea_by = config["gsea_by"],
            organism_go = config["organism_go"]
    log:    "logs/all_samples/{comparison}.{biotype}.DE_{analysis_type}.GSEA_WP.log"
    conda:  "../wrappers/GSEA_wp/env.yaml"
    script: "../wrappers/GSEA_wp/script.py"