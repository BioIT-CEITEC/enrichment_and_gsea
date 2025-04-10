run_all <- function(args){

  input_genes <- args[1]
  OUTPUT_DIR <- args[2]
  organism_wp <- args[3]
  n_up <- as.integer(args[4])
  COLORS <- unlist(strsplit(args[5],split=":"))[1]
  enrich_padj <- as.numeric(args[6])
  enrich_padjmethod <- args[7]
  enrich_minGSSize <- as.numeric(args[8])
  enrich_maxGSSize <- as.numeric(args[9])
  input_universe <- args[10]

  library("data.table")
  library("clusterProfiler")
  library("ggplot2")
  library("stringr")

  deseq2_tab <- fread(input_genes)
  deseq2_tab$ENTREZID <- as.character(deseq2_tab$ENTREZID)
  deseq2_tab <- unique(deseq2_tab)

  if(dir.exists(OUTPUT_DIR)==F){
    dir.create(OUTPUT_DIR, recursive = T)
  }

  emptytable<-data.table(ID=character(),Description=character(),GeneRatio=character(),BgRatio=character(),pvalue=numeric(),p.adjust=numeric(),qvalue=numeric(),geneID=character(),Count=integer())

  if(length(deseq2_tab$ENTREZID) == 0){
    # create an empty table
    dtewp<-emptytable
  }else{
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

    ewp <- enrichWP(gene        = deseq2_tab$ENTREZID,
                    universe      = universe$ENTREZID,
                    organism      = organism_wp,
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

    dtewp <- as.data.table(ewp)
    if(length(dtewp$ID) > 0){
      dtewpex <- convert_geneid(dtewp, deseq2_tab, is.gsea = F, is.entrez = T)
      fwrite(dtewpex, file = paste0(OUTPUT_DIR,"/enrich_WP_extended.tsv"), sep="\t")
    }else{
        dtewp<-emptytable
    }
  }
  fwrite(dtewp, file = paste0(OUTPUT_DIR,"/enrich_WP.tsv"), sep="\t")

  # Plot enrichment plot
  myEnrichPlot <- function(go.table = dtewp,
                         nUp = 10,
                         PADJ = 0.05,
                         mycol = "firebrick",
                         ploTitle = "WikiPathways"){
    go.table <- go.table[p.adjust <= PADJ,]
    setorder(go.table, p.adjust)

    nUp <- ifelse(nUp < 0, 10, nUp)
    nUp <- ifelse(nUp > length(go.table[, Count]),
                length(go.table[, Count]),nUp)

    filtRes <- head(go.table, n = nUp)

    if(length(filtRes$Description)>0){
      filtRes$nDescription <- str_wrap(filtRes$Description, width = 100)
    }else{
      filtRes[, nDescription := Description]
    }

    g <- ggplot(filtRes, aes(reorder(nDescription, -p.adjust), Count)) +
      geom_col(fill = mycol) +
      coord_flip() +
      labs(x="", y="Count",
           title=ploTitle) +
      theme_minimal() +
      theme(legend.position = "none")

    return(g)
    }

  WP_plot <- myEnrichPlot(dtewp,
                             nUp = n_up,
                             PADJ = enrich_padj,
                             mycol = COLORS,
                             ploTitle = "WikiPathways")
  ggsave(WP_plot, filename = paste0(OUTPUT_DIR,"/enrich_WP.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(WP_plot, filename = paste0(OUTPUT_DIR,"/enrich_WP.svg",sep=""),
       width = 10, height = 7, device = "svg")
  # ggsave(WP_plot, filename = paste0(OUTPUT_DIR,"/enrich_WP.png",sep=""),
  #      width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(10)
# args[1] <- "gene_for_enrichment.tsv" # input_genes
# args[2] <- "enrichment_WP" # OUTPUT_DIR
# args[3] <- "Mus musculus" # organism_wp
# args[4] <- 10 # n_up
# args[5] <- "firebrick:white:royalblue" # colors
# args[6] <- "0.1" # enrich_padj = padj from the enrich result table
# args[7] <- "BH" # enrich_padjmethod (BH,BY,fdr,holm,hochberg,hommel,bonferroni,none)
# args[8] <- "2" # enrich_minGSSize
# args[9] <- "Inf" # enrich_maxGSSize
# args[10] <- "gene_universe.tsv" # input_universe

run_all(args)
