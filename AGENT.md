# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Snakemake pipeline for gene enrichment analysis and GSEA (Gene Set Enrichment Analysis). It processes differential expression results from various quantification methods (featureCount, HTSeq-count, RSEM, Salmon, kallisto) and performs pathway enrichment analysis across multiple databases: Gene Ontology (GO), KEGG, Reactome, WikiPathways, and Enrichr.

The workflow depends on `BioIT-CEITEC/bioroots_utilities` GitHub module for sample loading and shared utilities.

## Running the Pipeline

```bash
# Run full pipeline (default target: final HTML report)
snakemake --cores <n>

# Run with specific config
snakemake --cores <n> --configfile config.json

# Dry run to check workflow
snakemake --cores <n> --dryrun

# View available rules
snakemake --cores <n> --list
```

## Configuration

1. Copy `config.example.json` to `config.json` and customize:
   - `entity_name`: Project identifier
   - `onthology`, `kegg`, `reactome`, `wikipathways`: Toggle enrichment databases
   - `conditions_to_compare`: "all" or comma-separated condition pairs (e.g., "KO:WT")
   - `biotypes`: Gene biotype filter (default: "all")
   - `cutoff_log2fc_enrich/gsea`: Log2FC thresholds for gene selection
   - `n_up/n_down`: Number of genes to enrich per direction
   - `colors`: Plot color scheme (format: "up:white:down")
   - `samples`: Sample metadata (sample_name, condition, replicate)

2. Enable quantification methods in `workflow.config.json`:
   - `featureCount`, `HTSeqCount`, `RSEM`, `salmon_map`, `salmon_align`, `kallisto`
   - For featureCount/HTSeq: specify `count_over` (exon, gene, transcript, 3pUTR, 5pUTR)

## Architecture

### Main Components

- **[Snakefile](snakefile)**: Entry point, loads config, defines wildcard constraints, includes rules
- **[rules/enrichment_analysis.smk](rules/enrichment_analysis.smk)**: Core enrichment rules (sampling, GO/KEGG/Reactome/WP enrichment + GSEA, Enrichr)
- **wrappers/**: Per-database wrappers containing:
  - `script*.py`: Snakemake shell wrappers that invoke R scripts
  - `*.R`: Actual enrichment/GSEA analysis using clusterProfiler/fgsea
  - `env.yaml`: Conda environment per wrapper

### Data Flow

```
DE results (DE_{analysis_type}/{comparison}/DESeq2.tsv)
    ↓ sampling rule (filters genes by cutoffs)
gene_for_enrichment_{up/down/all}.tsv + gene_for_gsea.tsv
    ↓ enrichment rules
enrichment_{DB}_{enrich}/enrich_{DB}.{tsv/svg}
    ↓ final_report rule
enrichment_gsea/enrichment_GSEA_final_report.html
```

### Analysis Types Supported

| Method | Wildcard Pattern |
|--------|------------------|
| featureCount | `featureCount_exon`, `featureCount_gene`, `featureCount_transcript`, `featureCount_3pUTR`, `featureCount_5pUTR` |
| HTSeq-count | `HTSeqCount_exon`, `HTSeqCount_gene`, `HTSeqCount_transcript`, `HTSeqCount_3pUTR`, `HTSeqCount_5pUTR` |
| RSEM | `RSEM` |
| Salmon (aligned) | `salmon_align` |
| Salmon (mapping-based) | `salmon_map` |
| kallisto | `kallisto` |

### Enrichment Databases

Each database has paired rules:
- `rule enrichment_*`: Over-representation analysis (ORA) for up/down/all genes
- `rule GSEA_*`: Pre-ranked GSEA using fgsea

Databases: GO (BP/MF/CC), KEGG, Reactome, WikiPathways, plus `gseapy_enrichr` for Enrichr API.

## Key Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `cutoff_log2fc_enrich` | Log2FC threshold for enrichment gene selection | 1 |
| `cutoff_padj_enrich` | Adjusted p-value threshold for enrichment | 0.05 |
| `cutoff_log2fc_gsea` | Log2FC threshold for GSEA ranking | 0 |
| `cutoff_padj_gsea` | Adjusted p-value threshold for GSEA | 1 |
| `enrich_padj` | Enrichment result adj. p-value cutoff | 0.1 |
| `gsea_nPermSimple` | Permutations for GSEA p-value estimation | 10000 |
| `gsea_by` | GSEA ranking method ("fgsea" or "DOSE") | fgsea |

## Debugging Tips

- Check logs in `logs/all_samples/{comparison}.DE_{analysis_type}.*.log`
- Sampling script filters DE results; verify input DESeq2.tsv exists
- KEGG requires local reference path set via `GLOBAL_REF_PATH` in config
- Ensure `bioRoot` samples sheet has correct condition labels matching `conditions_to_compare`

## Modifying the Workflow

To add a new enrichment database:
1. Create `wrappers/enrichment_newdb/` with `env.yaml`, `script_enrich.py`, `script_gsea.py`, `enrichment_newdb.R`, `GSEA_newdb.R`
2. Add rules to [rules/enrichment_analysis.smk](rules/enrichment_analysis.smk) following existing patterns
3. Update `final_input()` function to include new database outputs
4. Add organism references to `reference.json` or equivalent
