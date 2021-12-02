run_all <- function(args){

  WORKDIR <- args[1]
  input_genes <- args[2]
  OUTPUT_DIR <- args[3]
  organism <- args[4]
  cutoff_log2fc <- as.numeric(args[5])
  cutoff_padj <- as.numeric(args[6])
  n.up <- as.integer(args[7])
  COLORS <- unlist(strsplit(args[8],split=":"))[1]
  enrich_padj <- as.numeric(args[9])
  enrich_padjmethod <- args[10]
  enrich_minGSSize <- as.numeric(args[11])
  enrich_maxGSSize <- as.numeric(args[12])

  library("data.table")
  library("clusterProfiler")
  library("ggplot2")

  if(!require(organism, character.only = T)) {BiocManager::install(organism, update = F)}
  library(organism, character.only = T)
  database <- get(organism)
  if(organism == "org.At.tair.db"){
    KEYID <- "TAIR"
  }else{
    KEYID <- "ENSEMBL"
  }

  setwd(WORKDIR)

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

  egoBP <- enrichGO(gene          = deseq2_tab$Geneid,
                    universe      = universe[,KEYID],
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "BP", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    readable      = FALSE,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

  dtegoBP <- as.data.table(egoBP)
  fwrite(dtegoBP, file = paste0(OUTPUT_DIR,"/GO_enrich_BP.tsv"), sep="\t")

  egoMF <- enrichGO(gene          = deseq2_tab$Geneid,
                    universe      = universe[,KEYID],
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "MF", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    readable      = FALSE,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

  dtegoMF <- as.data.table(egoMF)
  fwrite(dtegoMF, file = paste0(OUTPUT_DIR,"/GO_enrich_MF.tsv"), sep="\t")

  egoCC <- enrichGO(gene          = deseq2_tab$Geneid,
                    universe      = universe[,KEYID],
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "CC", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    readable      = FALSE,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

  dtegoCC <- as.data.table(egoCC)
  fwrite(dtegoCC, file = paste0(OUTPUT_DIR,"/GO_enrich_CC.tsv"), sep="\t")

  # Plot enrichment plot
  myEnrichPlot <- function(go.table = dtegoBP,
                         nUp = 10,
                         GOPADJ = 0.05,
                         mycol = "firebrick",
                         ploTitle = "GO - Biological Processes"){
    go.table <- go.table[p.adjust <= GOPADJ,]
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

  GO_BP_plot <- myEnrichPlot(dtegoBP,
                             nUp = n.up,
                             GOPADJ = enrich_padj,
                             mycol = COLORS,
                             ploTitle = "GO - Biological Process")
  ggsave(GO_BP_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_BP.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(GO_BP_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_BP.svg",sep=""),
       width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GO_BP_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_BP.png",sep=""),
       width = 10, height = 7, device = "png", bg='transparent')

  GO_MF_plot <- myEnrichPlot(dtegoMF,
                           nUp = n.up,
                           GOPADJ = GOPADJ,
                           mycol = COLOR,
                           ploTitle = "GO - Molecular Function")
  ggsave(GO_MF_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_MF.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(GO_MF_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_MF.svg",sep=""),
       width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GO_MF_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_MF.png",sep=""),
       width = 10, height = 7, device = "png", bg='transparent')

  GO_CC_plot <- myEnrichPlot(dtegoCC,
                           nUp = n.up,
                           GOPADJ = GOPADJ,
                           mycol = COLOR,
                           ploTitle = "GO - Cellular Componenet")
  ggsave(GO_CC_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_CC.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(GO_CC_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_CC.svg",sep=""),
       width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GO_CC_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_CC.png",sep=""),
       width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(12)
# args[1] <- "E:/OneDrive - MUNI/TF_Daniel" # WORKDIR
# args[2] <- "DESeq2.tsv" # input_genes
# args[3] <- "clusterProfiler" # OUTPUT_DIR
# args[4] <- "org.Hs.eg.db" # organism
# args[5] <- "1" # cutoff_log2fc
# args[6] <- "0.05" # cutoff_padj
# args[7] <- 10 # n.up
# args[8] <- "firebrick:white:royalblue" # colors
# args[9] <- "0.1" # enrich_padj = padj from the enrich result table
# args[10] <- "BH" # enrich_padjmethod (BH,BY,fdr,holm,hochberg,hommel,bonferroni,none)
# args[11] <- "2" # enrich_minGSSize
# args[12] <- "Inf" # enrich_maxGSSize

run_all(args)
