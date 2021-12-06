run_all <- function(args){

  WORKDIR <- args[1]
  input_genes <- args[2]
  OUTPUT_DIR <- args[3]
  organism <- args[4]
  cutoff_log2fc <- as.numeric(args[5])
  cutoff_padj <- as.numeric(args[6])
  n_up <- as.integer(args[7])
  n_down <- as.integer(args[8])
  COLORS <- unlist(strsplit(args[9],split=":"))
  gsea_padj <- as.numeric(args[10])
  gsea_padjmethod <- args[11]
  gsea_minGSSize <- as.numeric(args[12])
  gsea_maxGSSize <- as.numeric(args[13])
  gsea_eps <- as.numeric(args[14])
  gsea_nPermSimple <- as.numeric(args[15])
  gsea_by <- args[16]

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
  deseq_cutoff <- function(deseq = deseq2_tab, LOG2FC = 0, PADJ = 1){
    x <- deseq[is.na(padj) == F & is.na(pvalue) == F,]
    x <- x[abs(log2FoldChange) >= LOG2FC & padj <= PADJ,]
    return(x)
  }
  deseq2_tab <- deseq_cutoff(de, cutoff_log2fc, cutoff_padj)

  ## select just entrez id and stat/logFC
  genes <- deseq2_tab[,.(Geneid, logFC = log2FoldChange)]
  ## remove NA values
  genes <- na.omit(genes)
  ## order by decreasing logFC
  setorder(genes, -logFC)

  ## change into named vector
  rankGenes <- genes$logFC
  names(rankGenes) <- genes$Geneid

    if(dir.exists(OUTPUT_DIR)==F){
    dir.create(OUTPUT_DIR, recursive = T)
  }


  if(dir.exists(OUTPUT_DIR)==F){
    dir.create(OUTPUT_DIR, recursive = T)
  }

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  universe <- select(database, keys=keys(database), columns = c(KEYID,'ENTREZID','SYMBOL'))

  deseq2_tab <- merge(deseq2_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)
  fwrite(deseq2_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = paste0(OUTPUT_DIR,"/Gene_ID.tsv"), sep="\t")

  gseGOBP <- gseGO(gene          = rankGenes,
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "BP", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    readable      = FALSE,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseGOBP <- as.data.table(gseGOBP)
  fwrite(dtgseGOBP, file = paste0(OUTPUT_DIR,"/GSEA_GO_BP.tsv"), sep="\t")

  gseGOMF <- gseGO(gene          = rankGenes,
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "MF", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    readable      = FALSE,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseGOMF <- as.data.table(gseGOMF)
  fwrite(dtgseGOMF, file = paste0(OUTPUT_DIR,"/GSEA_GO_MF.tsv"), sep="\t")

  gseGOCC <- gseGO(gene          = rankGenes,
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "CC", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    readable      = FALSE,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseGOCC <- as.data.table(gseGOCC)
  fwrite(dtgseGOCC, file = paste0(OUTPUT_DIR,"/GSEA_GO_CC.tsv"), sep="\t")

  # Plot enrichment plot
  myGSEAPlot <- function(fgsea.table = dtgseGOBP,
                         nUp = 10,
                         nDown = 10,
                         GOPADJ = 0.05,
                         gradient = c("firebrick","white","royalblue"),
                         ploTitle = "GO - Biological Processes"){
    fgsea.table <- fgsea.table[p.adjust <= Padj,]
    setorder(fgsea.table, p.adjust)
    fgsea.table[, Enrichment := ifelse(NES > 0, "Up-regulated", "Down-regulated")]

    nUp <- ifelse(nUp < 0, 0, nUp)
    nUp <- ifelse(nUp > length(fgsea.table[NES > 0, NES]),
                  length(fgsea.table[NES > 0, NES]),nUp)

    nDown <- ifelse(nDown < 0, 0, nDown)
    nDown <- ifelse(nDown > length(fgsea.table[NES < 0, NES]),
                    length(fgsea.table[NES < 0, NES]),nDown)

    filtRes <- rbind(head(fgsea.table, n = nUp),
                     tail(fgsea.table, n = nDown ))

    g <- ggplot(filtRes, aes(reorder(Description, NES), NES)) +
      geom_col( aes(fill = NES )) +
      coord_flip() +
      labs(x="", y="Normalized Enrichment Score",
           title=ploTitle) +
      theme_minimal() +
      scale_fill_gradient2(high = gradient[1],
                           mid = gradient[2],
                           low = gradient[3],
                           midpoint = 0) +
      theme(legend.position = "none")

    return(g)
    }

  GSEA_BP_plot <- myEnrichPlot(dtgseGOBP,
                               nUp = n_up,
                               nDown = n_dn,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - Biological Process")
  ggsave(GSEA_BP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_BP.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_BP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_BP.svg",sep=""),
         width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GSEA_BP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_BP.png",sep=""),
         width = 10, height = 7, device = "png", bg='transparent')

  GSEA_MF_plot <- myEnrichPlot(dtgseGOMF,
                               nUp = n_up,
                               nDown = n_dn,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - Molecular Function")
  ggsave(GSEA_MF_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_MF.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_MF_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_MF.svg",sep=""),
         width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GSEA_MF_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_MF.png",sep=""),
         width = 10, height = 7, device = "png", bg='transparent')

  GSEA_CC_plot <- myEnrichPlot(dtgseGOCC,
                               nUp = n_up,
                               nDown = n_dn,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - Cellular Component")
  ggsave(GSEA_CC_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_CC.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_CC_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_CC.svg",sep=""),
         width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GSEA_CC_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_CC.png",sep=""),
         width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(16)
# args[1] <- "E:/OneDrive - MUNI/TF_Daniel" # WORKDIR
# args[2] <- "DESeq2.tsv" # input_genes
# args[3] <- "GSEA_GO" # OUTPUT_DIR
# args[4] <- "org.Hs.eg.db" # organism
# args[5] <- "0" # cutoff_log2fc
# args[6] <- "1" # cutoff_padj
# args[7] <- 10 # n_up
# args[8] <- 10 # n_down
# args[9] <- "firebrick:white:royalblue" # colors
# args[10] <- "0.1" # gsea_padj = padj from the enrich result table
# args[11] <- "BH" # gsea_padjmethod (BH,BY,fdr,holm,hochberg,hommel,bonferroni,none)
# args[12] <- "2" # gsea_minGSSize
# args[13] <- "Inf" # gsea_maxGSSize
# args[14] <- "0" # gsea_eps
# args[15] <- 10000 # gsea_nPermSimple
# args[16] <- "fgsea" # gsea_by ("fgsea","DOSE")

run_all(args)
