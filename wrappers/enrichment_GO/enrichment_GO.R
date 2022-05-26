run_all <- function(args){

  input_genes <- args[1]
  OUTPUT_DIR <- args[2]
  organism_go <- args[3]
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

  if(!require(organism_go, character.only = T)) {BiocManager::install(organism_go, update = F)}
  library(organism_go, character.only = T)
  database <- get(organism_go)
  if(organism_go == "org.At.tair.db"){
    KEYID <- "TAIR"
  }else{
    KEYID <- "ENSEMBL"
  }

  # # read results of DE analysis
  # de <- fread(input_genes,header = T)
  # colnames(de)[colnames(de) == 'V1'] <- 'Geneid' # In case we need to rename
  # ## filter genes
  # deseq_cutoff <- function(deseq = deseq2_tab, LOG2FC = 1, PADJ = 0.05){
  #   x <- deseq[is.na(padj) == F & is.na(pvalue) == F,]
  #   x <- x[abs(log2FoldChange) >= LOG2FC & padj <= PADJ,]
  #   return(x)
  # }
  # deseq2_tab <- deseq_cutoff(de, cutoff_log2fc, cutoff_padj)
  deseq2_tab <- fread(input_genes)
  deseq2_tab$ENTREZID <- as.character(deseq2_tab$ENTREZID)

  if(dir.exists(OUTPUT_DIR)==F){
    dir.create(OUTPUT_DIR, recursive = T)
  }

  ## lookup gene symbol and unigene ID for the 1st 6 keys
  #universe <- select(database, keys=keys(database), columns = c(KEYID,'ENTREZID','SYMBOL'))
  universe <- fread(input_universe)
  universe$ENTREZID <- as.character(universe$ENTREZID)

  #deseq2_tab <- merge(deseq2_tab, universe, by.x = "Geneid", by.y = KEYID, all.x=T)
  #fwrite(deseq2_tab[,.(Geneid, gene_name, ENTREZID, log2FoldChange, padj)], file = paste0(OUTPUT_DIR,"/Gene_ID.tsv"), sep="\t")

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
  egoBP <- enrichGO(gene          = deseq2_tab$Geneid,
                    universe      = universe[,get(KEYID)],
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "BP", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    readable      = FALSE,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

  dtegoBP <- as.data.table(egoBP)
  dtegoBPex <- convert_geneid(dtegoBP, deseq2_tab, F, F)
  fwrite(dtegoBPex, file = paste0(OUTPUT_DIR,"/GO_enrich_BP_extended.tsv"), sep="\t")
  fwrite(dtegoBP, file = paste0(OUTPUT_DIR,"/GO_enrich_BP.tsv"), sep="\t")

  egoMF <- enrichGO(gene          = deseq2_tab$Geneid,
                    universe      = universe[,get(KEYID)],
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "MF", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    readable      = FALSE,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

  dtegoMF <- as.data.table(egoMF)
  dtegoMFex <- convert_geneid(dtegoMF, deseq2_tab, F, F)
  fwrite(dtegoMFex, file = paste0(OUTPUT_DIR,"/GO_enrich_MF_extended.tsv"), sep="\t")
  fwrite(dtegoMF, file = paste0(OUTPUT_DIR,"/GO_enrich_MF.tsv"), sep="\t")

  egoCC <- enrichGO(gene          = deseq2_tab$Geneid,
                    universe      = universe[,get(KEYID)],
                    OrgDb         = database,
                    keyType       = KEYID,
                    ont           = "CC", # "MF", "BP", and "CC", "ALL" (?)
                    pAdjustMethod = enrich_padjmethod,
                    pvalueCutoff  = enrich_padj,
                    readable      = FALSE,
                    minGSSize     = enrich_minGSSize,
                    maxGSSize     = enrich_maxGSSize)

  dtegoCC <- as.data.table(egoCC)
  dtegoCCex <- convert_geneid(dtegoCC, deseq2_tab, F, F)
  fwrite(dtegoCCex, file = paste0(OUTPUT_DIR,"/GO_enrich_CC_extended.tsv"), sep="\t")
  fwrite(dtegoCC, file = paste0(OUTPUT_DIR,"/GO_enrich_CC.tsv"), sep="\t")

  # Plot enrichment plot
  myEnrichPlot <- function(go.table = dtegoBP,
                         nUp = 10,
                         PADJ = 0.05,
                         mycol = "firebrick",
                         ploTitle = "GO - Biological Processes"){
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

  GO_BP_plot <- myEnrichPlot(dtegoBP,
                             nUp = n_up,
                             PADJ = enrich_padj,
                             mycol = COLORS,
                             ploTitle = "GO - Biological Process")
  ggsave(GO_BP_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_BP.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(GO_BP_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_BP.svg",sep=""),
       width = 10, height = 7, device = "svg")
  # ggsave(GO_BP_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_BP.png",sep=""),
  #      width = 10, height = 7, device = "png", bg='transparent')

  GO_MF_plot <- myEnrichPlot(dtegoMF,
                           nUp = n_up,
                           PADJ = enrich_padj,
                           mycol = COLORS,
                           ploTitle = "GO - Molecular Function")
  ggsave(GO_MF_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_MF.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(GO_MF_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_MF.svg",sep=""),
       width = 10, height = 7, device = "svg")
  # ggsave(GO_MF_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_MF.png",sep=""),
  #      width = 10, height = 7, device = "png", bg='transparent')

  GO_CC_plot <- myEnrichPlot(dtegoCC,
                           nUp = n_up,
                           PADJ = enrich_padj,
                           mycol = COLORS,
                           ploTitle = "GO - Cellular Component")
  ggsave(GO_CC_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_CC.pdf",sep=""),
       width = 10, height = 7, device = "pdf")
  ggsave(GO_CC_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_CC.svg",sep=""),
       width = 10, height = 7, device = "svg")
  # ggsave(GO_CC_plot, filename = paste0(OUTPUT_DIR,"/GO_enrich_CC.png",sep=""),
  #     width = 10, height = 7, device = "png", bg='transparent')

}

# run as Rscript
args <- commandArgs(trailingOnly = T)

# args <- character(10)
# args[1] <- "gene_for_enrichment.tsv" # input_genes
# args[2] <- "enrichment_GO" # OUTPUT_DIR
# args[3] <- "org.Hs.eg.db" # organism_go
# args[4] <- 10 # n_up
# args[5] <- "firebrick:white:royalblue" # colors
# args[6] <- "0.1" # enrich_padj = padj from the enrich result table
# args[7] <- "BH" # enrich_padjmethod (BH,BY,fdr,holm,hochberg,hommel,bonferroni,none)
# args[8] <- "2" # enrich_minGSSize
# args[9] <- "Inf" # enrich_maxGSSize
# args[10] <- "gene_universe.tsv" # input_universe

run_all(args)
