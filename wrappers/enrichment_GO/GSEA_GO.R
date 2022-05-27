run_all <- function(args){

  input_genes <- args[1]
  OUTPUT_DIR <- args[2]
  organism <- args[3]
  n_up <- as.integer(args[4])
  n_down <- as.integer(args[5])
  COLORS <- unlist(strsplit(args[6],split=":"))
  gsea_padj <- as.numeric(args[7])
  gsea_padjmethod <- args[8]
  gsea_minGSSize <- as.numeric(args[9])
  gsea_maxGSSize <- as.numeric(args[10])
  gsea_eps <- as.numeric(args[11])
  gsea_nPermSimple <- as.integer(args[12])
  gsea_by <- args[13]
  input_universe <- args[14]

  library("data.table")
  library("clusterProfiler")
  library("ggplot2")
  library("stringr")

  if(!require(organism, character.only = T)) {BiocManager::install(organism, update = F)}
  library(organism, character.only = T)
  database <- get(organism)
  if(organism == "org.At.tair.db"){
    KEYID <- "TAIR"
  }else{
    KEYID <- "ENSEMBL"
  }

  deseq2_tab <- fread(input_genes,header = T)
  deseq2_tab$ENTREZID <- as.character(deseq2_tab$ENTREZID)

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  universe <- fread(input_universe)
  universe$ENTREZID <- as.character(universe$ENTREZID)

  convert_geneid <- function(dt, deseq_tab = deseq2_tab, is.gsea = FALSE, is.entrez = FALSE){
    if (is.gsea == FALSE){
      tabl <- setDT(dt)[, strsplit(as.character(geneID), "/", fixed=TRUE),
                          by = .(ID, Description, pvalue, p.adjust, qvalue, geneID)
      ][,.(ID, Description, pvalue, p.adjust, qvalue, geneID = V1)]

      if (is.entrez == FALSE){
        tabl <- merge(tabl[, ENSEMBL := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENSEMBL", all.x=T)
      }
      else{
        tabl <- merge(tabl[, ENTREZID := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENTREZID", all.x=T)
      }
      tabl <- tabl[, .(ID, Description, pvalue, p.adjust, qvalue, ENSEMBL, gene_name, ENTREZID)]
      setorder(tabl, p.adjust, pvalue, ID, ENSEMBL)
    }
    else{
      tabl <- setDT(dt)[, strsplit(as.character(core_enrichment), "/", fixed=TRUE),
                          by = .(ID, Description, NES, pvalue, p.adjust, qvalues, core_enrichment)
      ][,.(ID, Description, NES, pvalue, p.adjust, qvalues, geneID = V1)]

      if (is.entrez == FALSE){
        tabl <- merge(tabl[, ENSEMBL := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENSEMBL", all.x=T)
      }
      else{
        tabl <- merge(tabl[, ENTREZID := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENTREZID", all.x=T)
      }
      tabl <- tabl[, .(ID, Description, NES, pvalue, p.adjust, qvalues, ENSEMBL, gene_name, ENTREZID)]
      setorder(tabl, p.adjust, pvalue, ID, ENSEMBL)
    }
    return(tabl)
  }

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

  gseaGOBP <- gseGO(gene          = rankGenes,
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "BP", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseaGOBP <- as.data.table(gseaGOBP)
  if(length(dtgseaGOBP$ID) > 0){
    dtgseaGOBPex <- convert_geneid(dtgseaGOBP, deseq2_tab, is.gsea = T, is.entrez = F)
    fwrite(dtgseaGOBPex, file = paste0(OUTPUT_DIR,"/GSEA_GO_BP_extended.tsv"), sep="\t")
  }
  fwrite(dtgseaGOBP, file = paste0(OUTPUT_DIR,"/GSEA_GO_BP.tsv"), sep="\t")

  gseaGOMF <- gseGO(gene          = rankGenes,
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "MF", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseaGOMF <- as.data.table(gseaGOMF)
  if(length(dtgseaGOMF$ID) > 0){
    dtgseaGOMFex <- convert_geneid(dtgseaGOMF, deseq2_tab, is.gsea = T, is.entrez = F)
    fwrite(dtgseaGOMFex, file = paste0(OUTPUT_DIR,"/GSEA_GO_MF_extended.tsv"), sep="\t")
  }
  fwrite(dtgseaGOMF, file = paste0(OUTPUT_DIR,"/GSEA_GO_MF.tsv"), sep="\t")

  gseaGOCC <- gseGO(gene          = rankGenes,
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "CC", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseaGOCC <- as.data.table(gseaGOCC)
  if(length(dtgseaGOCC$ID) > 0){
    dtgseaGOCCex <- convert_geneid(dtgseaGOCC, deseq2_tab, is.gsea = T, is.entrez = F)
    fwrite(dtgseaGOCCex, file = paste0(OUTPUT_DIR,"/GSEA_GO_CC_extended.tsv"), sep="\t")
  }
  fwrite(dtgseaGOCC, file = paste0(OUTPUT_DIR,"/GSEA_GO_CC.tsv"), sep="\t")

  # Plot enrichment plot
  myGSEAPlot <- function(fgsea.table = dtgseaGOBP,
                         nUp = 10,
                         nDown = 10,
                         Padj = 0.05,
                         gradient = c("firebrick","white","royalblue"),
                         ploTitle = "GO - Biological Processes"){
    fgsea.table <- fgsea.table[p.adjust <= Padj,]
    setorder(fgsea.table, -NES)
    fgsea.table[, Enrichment := ifelse(NES > 0, "Up-regulated", "Down-regulated")]

    nUp <- ifelse(nUp < 0, 0, nUp)
    nUp <- ifelse(nUp > length(fgsea.table[NES > 0, NES]),
                  length(fgsea.table[NES > 0, NES]),nUp)

    nDown <- ifelse(nDown < 0, 0, nDown)
    nDown <- ifelse(nDown > length(fgsea.table[NES < 0, NES]),
                    length(fgsea.table[NES < 0, NES]),nDown)

    filtRes <- rbind(head(fgsea.table, n = nUp),
                     tail(fgsea.table, n = nDown ))

    if(length(filtRes$ID) == 1){
      if(filtRes$NES > 0){
        gradient[2] <- gradient[1]
      }
      if(filtRes$NES < 0){
        gradient[2] <- gradient[3]
      }
    }

    if(length(filtRes$Description)>0){
      filtRes$nDescription <- str_wrap(filtRes$Description, width = 100)
    }else{
      filtRes[, nDescription := Description]
    }

    g <- ggplot(filtRes, aes(reorder(nDescription, NES), NES)) +
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

  GSEA_BP_plot <- myGSEAPlot(dtgseaGOBP,
                               nUp = n_up,
                               nDown = n_down,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - Biological Process")
  ggsave(GSEA_BP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_BP.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_BP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_BP.svg",sep=""),
         width = 10, height = 7, device = "svg")
  # ggsave(GSEA_BP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_BP.png",sep=""),
  #        width = 10, height = 7, device = "png", bg='transparent')

  GSEA_MF_plot <- myGSEAPlot(dtgseaGOMF,
                               nUp = n_up,
                               nDown = n_down,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - Molecular Function")
  ggsave(GSEA_MF_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_MF.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_MF_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_MF.svg",sep=""),
         width = 10, height = 7, device = "svg")
  # ggsave(GSEA_MF_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_MF.png",sep=""),
  #        width = 10, height = 7, device = "png", bg='transparent')

  GSEA_CC_plot <- myGSEAPlot(dtgseaGOCC,
                               nUp = n_up,
                               nDown = n_down,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - Cellular Component")
  ggsave(GSEA_CC_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_CC.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_CC_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_CC.svg",sep=""),
         width = 10, height = 7, device = "svg")
  # ggsave(GSEA_CC_plot, filename = paste0(OUTPUT_DIR,"/GSEA_GO_CC.png",sep=""),
  #        width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(14)
# args[1] <- "gene_for_enrichment.tsv" # input_genes
# args[2] <- "GSEA_GO" # OUTPUT_DIR
# args[3] <- "org.Hs.eg.db" # organism
# args[4] <- 10 # n_up
# args[5] <- 10 # n_down
# args[6] <- "firebrick:white:royalblue" # colors
# args[7] <- "0.1" # gsea_padj = padj from the enrich result table
# args[8] <- "BH" # gsea_padjmethod (BH,BY,fdr,holm,hochberg,hommel,bonferroni,none)
# args[9] <- "2" # gsea_minGSSize
# args[10] <- "Inf" # gsea_maxGSSize
# args[11] <- "0" # gsea_eps
# args[12] <- "10000" # gsea_nPermSimple
# args[13] <- "fgsea" # gsea_by ("fgsea","DOSE")
# args[14] <- "gene_universe.tsv"

run_all(args)
