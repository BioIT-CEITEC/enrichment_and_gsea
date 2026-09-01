# Gene Enrichment and GSEA Analysis Pipeline

A Snakemake-based workflow for functional enrichment analysis and Gene Set Enrichment Analysis (GSEA) of RNA-seq differential expression results.

## Overview

This pipeline takes differential expression results from various quantification methods and performs:

- **Over-representation Analysis (ORA)** - Testing if specific gene sets are enriched among differentially expressed genes
- **Gene Set Enrichment Analysis (GSEA)** - Pre-ranked GSEA to detect subtle coordinated changes

### Supported Enrichment Databases

| Database | Description |
|----------|-------------|
| **Gene Ontology (GO)** | Biological Process, Molecular Function, Cellular Component |
| **KEGG** | Kyoto Encyclopedia of Genes and Genomes pathways |
| **Reactome** | Curated pathway database |
| **WikiPathways** | Community-curated pathway resource |
| **Enrichr** | Web-based gene set enrichment analysis |

### Supported Quantification Methods

- featureCount (exon, gene, transcript, 3' UTR, 5' UTR)
- HTSeq-count (exon, gene, transcript, 3' UTR, 5' UTR)
- RSEM
- Salmon (alignment-based and mapping-based)
- kallisto

## Quick Start

### 1. Setup Configuration

```bash
cp config.example.json config.json
```

Edit `config.json` to define your samples and analysis parameters.

### 2. Run the Pipeline

```bash
# Dry run to verify configuration
snakemake --cores 4 --dryrun

# Execute the full pipeline
snakemake --cores 8

# Generate HTML report
snakemake --cores 4 enrichment_gsea/enrichment_GSEA_final_report.html
```

## Configuration

### Required Settings

```json
{
  "entity_name": "project_name",
  "onthology": true,
  "kegg": true,
  "reactome": true,
  "wikipathways": true,
  "conditions_to_compare": "all",
  "samples": {
    "sample_id": {
      "sample_name": "display_name",
      "condition": "experimental_condition",
      "tag": "replicate_number"
    }
  }
}
```

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `cutoff_log2fc_enrich` | Log2FC threshold for gene selection in enrichment | 1 |
| `cutoff_padj_enrich` | Adjusted p-value threshold for enrichment | 0.05 |
| `cutoff_log2fc_gsea` | Log2FC threshold for GSEA ranking | 0 |
| `cutoff_padj_gsea` | Adjusted p-value threshold for GSEA | 1 |
| `n_up` / `n_down` | Number of genes to analyze per direction | 15 |
| `enrich_padj` | Significance cutoff for enriched terms | 0.1 |
| `gsea_nPermSimple` | Permutations for GSEA p-values | 10000 |
| `colors` | Color scheme for plots (up:white:down) | firebrick:white:royalblue |

### Enabling Quantification Methods

Set the desired methods to `true` in your configuration:

```json
{
  "featureCount": true,
  "HTSeqCount": false,
  "RSEM": false,
  "salmon_map": false,
  "salmon_align": false,
  "kallisto": false
}
```

For featureCount/HTSeqCount, also specify counting level:

```json
{
  "count_over": "gene"  // or: exon, transcript, three_prime_UTR, five_prime_UTR
}
```

## Output Structure

```
enrichment_gsea/
├── DE_{analysis_type}/{comparison}/
│   ├── gene_for_enrichment_*.tsv      # Gene lists for enrichment
│   ├── gene_for_gsea.tsv              # Ranked gene list for GSEA
│   ├── enrichment_GO_{enrich}/        # GO enrichment results + plots
│   ├── enrichment_KEGG_{enrich}/      # KEGG results + plots
│   ├── enrichment_REACTOME_{enrich}/  # Reactome results + plots
│   ├── enrichment_WP_{enrich}/        # WikiPathways results + plots
│   ├── GSEA_GO/                       # GO GSEA results + plots
│   ├── GSEA_KEGG/                     # KEGG GSEA results + plots
│   ├── GSEA_REACTOME/                 # Reactome GSEA results + plots
│   ├── GSEA_WP/                       # WikiPathways GSEA results + plots
│   └── enrichr_{enrich}/              # Enrichr results + plots
└── enrichment_GSEA_final_report.html  # Summary HTML report
```

Each enrichment output includes:
- `.tsv` - Tab-separated results table with p-values, odds ratios, gene lists
- `.svg` - Publication-ready bar plots of enriched terms

## Dependencies

The pipeline uses Conda environments defined in each wrapper directory. Snakemake will automatically create these environments when running.

### Core Dependencies
- Snakemake >= 5.18.0
- R with clusterProfiler, fgsea, DOSE, ggplot2
- Python with pandas, enrichr

### External Dependencies
- BioIT-CEITEC/bioroots_utilities (loaded via GitHub)
- Local reference files for organism-specific databases

## Requirements

- **Input**: Differential expression results in `DE_{analysis_type}/{comparison}/DESeq2.tsv` format
- **References**: Organism-specific annotation files must be available in the configured reference path

## License

Internal use at CEITEC Bioinformatics.

## Support

For issues or questions, contact the BioIT-CEITEC team.
