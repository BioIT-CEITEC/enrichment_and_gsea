run_all <- function(args){

  input_genes <- args[1]
  output_enrich <- args[2]
  output_gsea <- args[3]
  output_universe <- args[4]
  organism <- args[5]
  cutoff_log2fc_enrich <- as.numeric(args[6])
  cutoff_padj_enrich <- as.numeric(args[7])
  cutoff_log2fc_gsea <- as.numeric(args[8])
  cutoff_padj_gsea <- as.numeric(args[9])

  output_enrich_up <- paste0(sub(".tsv", "", output_enrich), "_up.tsv")
  output_enrich_down <- paste0(sub(".tsv", "", output_enrich), "_down.tsv")

  library("data.table")

  if(!require(organism, character.only = T)) {BiocManager::install(organism, update = F)}
  library(organism, character.only = T)
  database <- get(organism)
  if(organism == "org.At.tair.db"){
    KEYID <- "TAIR"
  }else{
    KEYID <- "ENSEMBL"
  }

  # read results of DE analysis
  de <- fread(input_genes,header = T)
  if(names(de)[1]!="Geneid"){names(de)[1]<-"Geneid"} # In case we need to rename

  ## filter genes
  deseq_cutoff <- function(deseq = deseq2_tab, LOG2FC = 0, PADJ = 1){
    x <- deseq[is.na(padj) == F & is.na(pvalue) == F,]
    x <- x[abs(log2FoldChange) >= LOG2FC & padj <= PADJ,]
    return(x)
  }
  enrich_tab <- deseq_cutoff(de, cutoff_log2fc_enrich, cutoff_padj_enrich)
  gsea_tab <- deseq_cutoff(de, cutoff_log2fc_gsea, cutoff_padj_gsea)

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  universe <- as.data.table(select(database, keys=keys(database), columns = c(KEYID,'ENTREZID','SYMBOL')))

  enrich_tab <- merge(enrich_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)
  gsea_tab <- merge(gsea_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)

  enrich_tab_up <- enrich_tab[log2FoldChange > 0,]
  enrich_tab_down <- enrich_tab[log2FoldChange < 0,]

  if("gene_name" %in% names(de)){
    fwrite(enrich_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = output_enrich, sep="\t")
    fwrite(enrich_tab_up[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = output_enrich_up, sep="\t")
    fwrite(enrich_tab_down[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = output_enrich_down, sep="\t")

    fwrite(gsea_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = output_gsea, sep="\t")
  }else{
    fwrite(enrich_tab[,.(Geneid, gene_name=Feature_name, ENTREZID, log2FoldChange, padj)], file = output_enrich, sep="\t")
    fwrite(enrich_tab_up[,.(Geneid, gene_name=Feature_name, ENTREZID, log2FoldChange, padj)], file = output_enrich_up, sep="\t")
    fwrite(enrich_tab_down[,.(Geneid, gene_name=Feature_name, ENTREZID, log2FoldChange, padj)], file = output_enrich_down, sep="\t")

    fwrite(gsea_tab[,.(Geneid, gene_name=Feature_name, ENTREZID, log2FoldChange, padj)], file = output_gsea, sep="\t")
  }

  fwrite(universe, file = output_universe, sep="\t")

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(9)
# args[1] <- "DESeq2.tsv" # input_genes
# args[2] <- "gene_for_enrichment.tsv" # output_enrich
# args[3] <- "gene_for_gsea.tsv" # output_gsea
# args[4] <- "gene_universe.tsv" # output_universe
# args[5] <- "org.Hs.eg.db" # organism
# args[6] <- "1" # cutoff_log2fc_enrich
# args[7] <- "0.05" # cutoff_padj_enrich
# args[8] <- "0" # cutoff_log2fc_gsea
# args[9] <- "1" # cutoff_padj_gsea

run_all(args)
