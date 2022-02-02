run_all <- function(args){

  WORKDIR <- args[1]
  input_genes <- args[2]
  organism <- args[3]
  cutoff_log2fc_enrich <- as.numeric(args[4])
  cutoff_padj_enrich <- as.numeric(args[5])
  cutoff_log2fc_gsea <- as.numeric(args[6])
  cutoff_padj_gsea <- as.numeric(args[7])

  library("data.table")
  library("clusterProfiler")

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
  colnames(de)[colnames(de) == 'V1'] <- 'Geneid' # In case we need to rename
  ## filter genes
  deseq_cutoff <- function(deseq = deseq2_tab, LOG2FC = 0, PADJ = 1){
    x <- deseq[is.na(padj) == F & is.na(pvalue) == F,]
    x <- x[abs(log2FoldChange) >= LOG2FC & padj <= PADJ,]
    return(x)
  }
  enrich_tab <- deseq_cutoff(de, cutoff_log2fc_enrich, cutoff_padj_enrich)
  gsea_tab <- deseq_cutoff(de, cutoff_log2fc_gsea, cutoff_padj_gsea)

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  universe <- select(database, keys=keys(database), columns = c(KEYID,'ENTREZID','SYMBOL'))

  enrich_tab <- merge(enrich_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)
  gsea_tab <- merge(gsea_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)

  fwrite(enrich_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = paste0(WORKDIR,"/gene_for_enrichment.tsv"), sep="\t")

  fwrite(gsea_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = paste0(WORKDIR,"/gene_for_gsea.tsv"), sep="\t")

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(16)
# args[1] <- "E:/OneDrive - MUNI/TF_Daniel" # WORKDIR
# args[2] <- "DESeq2.tsv" # input_genes
# args[3] <- "org.Hs.eg.db" # organism
# args[4] <- "1" # cutoff_log2fc_enrich
# args[5] <- "0.05" # cutoff_padj_enrich
# args[6] <- "0" # cutoff_log2fc_gsea
# args[7] <- "1" # cutoff_padj_gsea

run_all(args)
