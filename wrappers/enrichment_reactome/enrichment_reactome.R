run_all <- function(args){

  WORKDIR <- args[1]
  input_genes <- args[2]
  OUTPUT_DIR <- args[3]
  organism_reactome <- args[4]
  cutoff_log2fc <- as.numeric(args[5])
  cutoff_padj <- as.numeric(args[6])
  n_up <- as.integer(args[7])
  COLORS <- unlist(strsplit(args[8],split=":"))[1]
  enrich_padj <- as.numeric(args[9])
  enrich_padjmethod <- args[10]
  enrich_minGSSize <- as.numeric(args[11])
  enrich_maxGSSize <- as.numeric(args[12])
  organism_go <- args[13]

  library("data.table")
  library("clusterProfiler")
  library("ReactomePA")
  library("ggplot2")

  if(!require(organism_go, character.only = T)) {BiocManager::install(organism_go, update = F)}
  library(organism_go, character.only = T)
  database <- get(organism_go)
  if(organism_go == "org.At.tair.db"){
    KEYID <- "TAIR"
  }else{
    KEYID <- "ENSEMBL"
  }

  #setwd(WORKDIR)

  # read results of DE analysis
  de <- fread(input_genes,header = T)
  colnames(de)[colnames(de) == 'V1'] <- 'Geneid' # In case we need to rename
  ## filter genes
  deseq_cutoff <- function(deseq = deseq2_tab, LOG2FC = 1, PADJ = 0.05){
    x <- deseq[is.na(padj) == F & is.na(pvalue) == F,]
    x <- x[abs(log2FoldChange) >= LOG2FC & padj <= PADJ,]
    return(x)
  }
  deseq2_tab <- deseq_cutoff(de, cutoff_log2fc, cutoff_padj)

  if(dir.exists(OUTPUT_DIR)==F){
    dir.create(OUTPUT_DIR, recursive = T)
  }

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  universe <- select(database, keys=keys(database), columns = c(KEYID,'ENTREZID','SYMBOL'))

  deseq2_tab <- merge(deseq2_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)

  fwrite(deseq2_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = paste0(OUTPUT_DIR,"/Gene_ID.tsv"), sep="\t")

  ereact <- enrichPathway(gene        = deseq2_tab$ENTREZID,
                        universe      = universe$ENTREZID,
                        organism      = organism_reactome,
                        pAdjustMethod = enrich_padjmethod,
                        pvalueCutoff  = enrich_padj,
                        minGSSize     = enrich_minGSSize,
                        maxGSSize     = enrich_maxGSSize)

  dtereact <- as.data.table(ereact)
  fwrite(dtereact, file = paste0(OUTPUT_DIR,"/REACTOME_enrich.tsv"), sep="\t")

  # Plot enrichment plot
  myEnrichPlot <- function(go.table = dtereact,
                         nUp = 10,
                         PADJ = 0.05,
                         mycol = "firebrick",
                         ploTitle = "REACTOME pathways"){
    go.table <- go.table[p.adjust <= PADJ,]
    setorder(go.table, p.adjust)

    nUp <- ifelse(nUp < 0, 10, nUp)
    nUp <- ifelse(nUp > length(go.table[, Count]),
                length(go.table[, Count]),nUp)

    filtRes <- head(go.table, n = nUp)

    g <- ggplot(filtRes, aes(reorder(Description, -p.adjust), Count)) +
      geom_col(fill = mycol) +
      coord_flip() +
      labs(x="", y="Count",
           title=ploTitle) +
      theme_minimal() +
      theme(legend.position = "none")

    return(g)
    }

  REACTOME_plot <- myEnrichPlot(dtereact,
                             nUp = n_up,
                             PADJ = enrich_padj,
                             mycol = COLORS,
                             ploTitle = "REACTOME pathways")
  ggsave(REACTOME_plot, filename = paste0(OUTPUT_DIR,"/REACTOME_enrich.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(REACTOME_plot, filename = paste0(OUTPUT_DIR,"/REACTOME_enrich.svg",sep=""),
       width = 10, height = 7, device = svg, bg='transparent')
  ggsave(REACTOME_plot, filename = paste0(OUTPUT_DIR,"/REACTOME_enrich.png",sep=""),
       width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(13)
# args[1] <- "E:/OneDrive - MUNI/TF_Daniel" # WORKDIR
# args[2] <- "DESeq2.tsv" # input_genes
# args[3] <- "enrichment_REACTOME" # OUTPUT_DIR
# args[4] <- "human" # organism_reactome
# args[5] <- "1" # cutoff_log2fc
# args[6] <- "0.05" # cutoff_padj
# args[7] <- 15 # n_up
# args[8] <- "firebrick:white:royalblue" # colors
# args[9] <- "0.1" # enrich_padj = padj from the enrich result table
# args[10] <- "BH" # enrich_padjmethod (BH,BY,fdr,holm,hochberg,hommel,bonferroni,none)
# args[11] <- "2" # enrich_minGSSize
# args[12] <- "Inf" # enrich_maxGSSize
# args[13] <- "org.Hs.eg.db" # organism_go

run_all(args)
