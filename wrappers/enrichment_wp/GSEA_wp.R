run_all <- function(args){

  input_genes <- args[1]
  OUTPUT_DIR <- args[2]
  organism_wp <- args[3]
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

  deseq2_tab <- fread(input_genes,header = T)
  deseq2_tab$ENTREZID <- as.character(deseq2_tab$ENTREZID)
  deseq2_tab <- unique(deseq2_tab)

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
      }else{
        tabl <- merge(tabl[, ENTREZID := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENTREZID", all.x=T)
      }
      tabl <- tabl[, .(ID, Description, pvalue, p.adjust, qvalue, ENSEMBL, gene_name, ENTREZID)]
      setorder(tabl, p.adjust, pvalue, ID, ENSEMBL)
    }else{
      tabl <- setDT(dt)[, strsplit(as.character(core_enrichment), "/", fixed=TRUE),
                          by = .(ID, Description, NES, pvalue, p.adjust, qvalues, core_enrichment)
      ][,.(ID, Description, NES, pvalue, p.adjust, qvalues, geneID = V1)]

      if (is.entrez == FALSE){
        tabl <- merge(tabl[, ENSEMBL := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENSEMBL", all.x=T)
      }else{
        tabl <- merge(tabl[, ENTREZID := geneID], deseq_tab[, .(ENSEMBL = Geneid, gene_name, ENTREZID)],
                      by="ENTREZID", all.x=T)
      }
      tabl <- tabl[, .(ID, Description, NES, pvalue, p.adjust, qvalues, ENSEMBL, gene_name, ENTREZID)]
      setorder(tabl, p.adjust, pvalue, ID, ENSEMBL)
    }
    return(tabl)
  }

  ## select just entrez id and stat/logFC
  genes <- deseq2_tab[,.(ENTREZID, logFC = log2FoldChange)]
  ## remove NA values
  genes <- na.omit(genes)
  ## order by decreasing logFC
  setorder(genes, -logFC)

  ## change into named vector
  rankGenes <- genes$logFC
  names(rankGenes) <- genes$ENTREZID

  if(dir.exists(OUTPUT_DIR)==F){
    dir.create(OUTPUT_DIR, recursive = T)
  }

  gseaWP <- gseWP(gene            = rankGenes,
                    organism      = organism_wp,
                    pAdjustMethod = gsea_padjmethod,
                    pvalueCutoff  = gsea_padj,
                    minGSSize     = gsea_minGSSize,
                    maxGSSize     = gsea_maxGSSize,
                    nPermSimple   = gsea_nPermSimple,
                    eps           = gsea_eps,
                    by            = gsea_by)

  dtgseaWP <- as.data.table(gseaWP)
  if(length(dtgseaWP$ID) > 0){
    dtgseaWPex <- convert_geneid(dtgseaWP, deseq2_tab, is.gsea = T, is.entrez = T)
    fwrite(dtgseaWPex, file = paste0(OUTPUT_DIR,"/GSEA_WP_extended.tsv"), sep="\t")
  }
  fwrite(dtgseaWP, file = paste0(OUTPUT_DIR,"/GSEA_WP.tsv"), sep="\t")

  # Plot enrichment plot
  myGSEAPlot <- function(fgsea.table = dtgseaWP,
                         nUp = 10,
                         nDown = 10,
                         Padj = 0.05,
                         gradient = c("firebrick","white","royalblue"),
                         ploTitle = "GSEA - WikiPathways"){
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

  GSEA_WP_plot <- myGSEAPlot(dtgseaWP,
                               nUp = n_up,
                               nDown = n_down,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - WikiPathways")
  ggsave(GSEA_WP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_WP.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_WP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_WP.svg",sep=""),
         width = 10, height = 7, device = "svg")
  # ggsave(GSEA_WP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_WP.png",sep=""),
  #        width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(14)
# args[1] <- "gene_for_enrichment.tsv" # input_genes
# args[2] <- "GSEA_GO" # OUTPUT_DIR
# args[3] <- "Mus musculus" # organism_wp
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
