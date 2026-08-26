#############################################################
# wrapper for rule: enrichment_enrichr
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: enrichment_enrichR \n##\n")
f.close()

# Create output directory
command = "mkdir -p " + snakemake.params.outdir
f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)

# Generate gene list
command = "cut -f 2 " + snakemake.input.tsv + " | tail -n +2 > " + snakemake.params.gene_list
f = open(log_filename, 'a+')
f.write("## COMMAND: " + command + "\n")
shell(command)

# Check if gene list is empty
if os.path.getsize(snakemake.params.gene_list) > 0:
    # Run gseapy enrichr if gene list is not empty
    command = "gseapy enrichr -f svg" + \
              " -i " + snakemake.params.gene_list + \
              " -o " + snakemake.params.outdir + \
              " --org " + snakemake.params.enrichr_org + \
              " -g " + snakemake.params.enrichr_db + \
              " -t " + snakemake.params.n_up + " " + \
              " --cut " + snakemake.params.enrich_padj
    f = open(log_filename, 'a+')
    f.write("## COMMAND: " + command + "\n")

    # Run gseapy and capture stderr to check for ValueError
    # Use "true ||" prefix to prevent shell from failing on non-zero exit code
    error_log = snakemake.params.outdir + "/gseapy_error.tmp"
    command_with_capture = "( " + command + " 2> " + error_log + " ) || true"
    shell(command_with_capture)

    # Check if ValueError occurred (no enrich terms after cutoff) by reading the error log
    has_error = False
    if os.path.exists(error_log):
        with open(error_log, 'r') as errf:
            error_content = errf.read()
        if 'ValueError' in error_content or 'No enrich terms' in error_content:
            has_error = True

    if has_error:
        f.write("## WARNING: gseapy enrichr found no enrich terms after cutoff. Creating empty outputs.\n")
        # Create empty plot file
        shell("touch " + snakemake.output.plot)
        # Create empty TAB-delimited output file with header
        header = "Gene_set\tTerm\tOverlap\tP-value\tAdjusted P-value\tOld P-value\tOld Adjusted P-value\tOdds Ratio\tCombined Score\tGenes\n"
        with open(snakemake.output.table, 'w') as outf:
            outf.write(header)
        f.write("## Created empty output files.\n")
        shell("rm " + error_log)
    else:
        # Append logs and clean up
        shell("cat " + snakemake.params.outdir + "/gseapy.enrichr.*.log >> " + log_filename)
        shell("rm " + snakemake.params.outdir + "/gseapy.enrichr.*.log")
        shell("rm -f " + error_log)
else:
    # Create an empty plot file if gene list is empty
    command = "touch " + snakemake.output.plot
    f = open(log_filename, 'a+')
    f.write("## COMMAND: " + command + "\n")
    shell(command)

    # Create an empty TAB-delimited output file with the specified header
    header = "Gene_set\tTerm\tOverlap\tP-value\tAdjusted P-value\tOld P-value\tOld Adjusted P-value\tOdds Ratio\tCombined Score\tGenes\n"
    with open(snakemake.output.table, 'w') as f:
        f.write(header)

    f = open(log_filename, 'a+')
    f.write("## Created empty TAB-delimited output file with header: " + header)
    f.close()