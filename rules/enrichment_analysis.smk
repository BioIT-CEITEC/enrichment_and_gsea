def de_computation_input(wildcards):
    input = {}
    path = r"results"
    deseq_files = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if "DESeq2.tsv" in name:
                deseq_files.append(os.path.join(root,name))

    analysis_type = [d.rsplit("/",4)[1].split("DE_")[1] for d in deseq_files]
    comparison_dir_list = [d.rsplit("/", 4)[2] for d in deseq_files]
    biotype_dir_list = [d.rsplit("/", 4)[3] for d in deseq_files]

    # comparison_dir_list = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path,d)) and "_vs_" in d]
    # biotype_dir_list = []
    # for comp in comparison_dir_list:
    #     biotype_dir_list.extend([d for d in os.listdir(os.path.join(path,comp))])

    input['tsv'] = expand("results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv", comparison=comparison_dir_list, biotype=biotype_dir_list)

    if len(comparison_dir_list) > 0:
        print(comparison_dir_list)
    else:
        raise ValueError("There is no differential analysis results!")

    return input


rule enrichment_GO:
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO/GO_enrich_CC.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir= "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_GO",
            organism_go = config["organism_go"],
            cutoff_log2fc = config["cutoff_log2fc_enrich"],
            cutoff_padj = config["cutoff_padj_enrich"],
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
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO/GSEA_GO_CC.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_GO",
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
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_KEGG/KEGG_enrich.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_KEGG",
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
    conda:  "../wrappers/enrichment_KEGG/env.yaml"
    script: "../wrappers/enrichment_KEGG/script.py"

rule GSEA_kegg:
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_KEGG/GSEA_KEGG.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_KEGG",
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
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_REACTOME/REACTOME_enrich.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_REACTOME",
            organism_kegg = config["organism_reactome"],
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
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_REACTOME/GSEA_REACTOME.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_REACTOME",
            organism_kegg = config["organism_reactome"],
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
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_WP/WP_enrich.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/enrichment_WP",
            organism_kegg = config["organism_wp"],
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
    input:  tsv = "results/DE_{{analysis_type}}/{comparison}/{biotype}/DESeq2.tsv"
    output: plot = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_WP/GSEA_WP.png"
    params: workdir = "results/DE_{analysis_type}/{comparison}/{biotype}",
            outdir = "results/DE_{analysis_type}/{comparison}/{biotype}/GSEA_WP",
            organism_kegg = config["organism_wp"],
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