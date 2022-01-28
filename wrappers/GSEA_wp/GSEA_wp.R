run_all <- function(args){

  WORKDIR <- args[1]
  input_genes <- args[2]
  OUTPUT_DIR <- args[3]
  organism_wp <- args[4]
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
  gsea_nPermSimple <- as.integer(args[15])
  gsea_by <- args[16]
  organism_go <- args[17]

  library("data.table")
  library("clusterProfiler")
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
  deseq_cutoff <- function(deseq = deseq2_tab, LOG2FC = 0, PADJ = 1){
    x <- deseq[is.na(padj) == F & is.na(pvalue) == F,]
    x <- x[abs(log2FoldChange) >= LOG2FC & padj <= PADJ,]
    return(x)
  }
  deseq2_tab <- deseq_cutoff(de, cutoff_log2fc, cutoff_padj)

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  universe <- select(database, keys=keys(database), columns = c(KEYID,'ENTREZID','SYMBOL'))

  deseq2_tab <- merge(deseq2_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)
  fwrite(deseq2_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = paste0(OUTPUT_DIR,"/Gene_ID.tsv"), sep="\t")

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

  GSEA_WP_plot <- myGSEAPlot(dtgseaWP,
                               nUp = n_up,
                               nDown = n_down,
                               Padj = gsea_padj,
                               gradient = COLORS,
                               ploTitle = "GSEA - WikiPathways")
  ggsave(GSEA_WP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_WP.pdf",sep=""),
         width = 10, height = 7, device = "pdf")
  ggsave(GSEA_WP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_WP.svg",sep=""),
         width = 10, height = 7, device = svg, bg='transparent')
  ggsave(GSEA_WP_plot, filename = paste0(OUTPUT_DIR,"/GSEA_WP.png",sep=""),
         width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(17)
# args[1] <- "E:/OneDrive - MUNI/TF_Daniel" # WORKDIR
# args[2] <- "DESeq2.tsv" # input_genes
# args[3] <- "GSEA_WP" # OUTPUT_DIR
# args[4] <- "Homo sapiens" # organism_wp
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
# args[17] <- "org.Hs.eg.db" # organism_go

run_all(args)
