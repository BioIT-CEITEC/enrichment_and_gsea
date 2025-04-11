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
    shell(command)

    # Append logs and clean up
    command = "cat " + snakemake.params.outdir + "/gseapy.enrichr.*.log >> " + log_filename
    shell(command)

    command = "rm " + snakemake.params.outdir + "/gseapy.enrichr.*.log"
    shell(command)
else:
    # Create an empty plot file if gene list is empty
    command = "touch " + snakemake.output.plot
    f = open(log_filename, 'a+')
    f.write("## COMMAND: " + command + "\n")
    shell(command)
