#############################################################
# wrapper for rule: enrichment_reactome
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: enrichment_reactome \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/enrichment_reactome.R " +\
            snakemake.input.tsv + " " +\
            snakemake.params.outdir + " " +\
            snakemake.params.organism_reactome + " " +\
            snakemake.params.n_up + " " +\
            snakemake.params.colors + " " +\
            snakemake.params.enrich_padj + " " +\
            snakemake.params.enrich_padjmethod + " " +\
            snakemake.params.enrich_minGSSize + " " +\
            snakemake.params.enrich_maxGSSize + " " +\
            snakemake.params.universe +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
