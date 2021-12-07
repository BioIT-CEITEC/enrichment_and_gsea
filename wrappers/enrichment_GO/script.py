#############################################################
# wrapper for rule: enrichment_GO
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: enrichment_GO \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/enrichment_GO.R " +\
            snakemake.params.workdir + " " +\
            snakemake.input.tsv + " " +\
            snakemake.params.outdir + " " +\
            snakemake.params.organism_go + " " +\
            snakemake.params.cutoff_log2fc + " " +\
            snakemake.params.cutoff_padj + " " +\
            snakemake.params.n_up + " " +\
            snakemake.params.colors + " " +\
            snakemake.params.enrich_padj + " " +\
            snakemake.params.enrich_padjmethod + " " +\
            snakemake.params.enrich_minGSSize + " " +\
            snakemake.params.enrich_maxGSSize +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
