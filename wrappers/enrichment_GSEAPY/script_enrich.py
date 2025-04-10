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

command = " gseapy enrichr -f svg" +\
            " -i " + snakemake.input.tsv +\
            " -o " + snakemake.params.outdir +\
            " --org " + snakemake.params.enrichr_org +\
            " -g " + snakemake.params.enrichr_db +\
            " -t " snakemake.params.n_up + " " +\
            " --cut " snakemake.params.enrich_padj

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " cat " + snakemake.params.outdir + "/gseapy.enrichr.*.log  >> " + snakemake.log
shell(command)
