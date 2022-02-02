#############################################################
# wrapper for rule: sampling
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: sampling \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/sampling.R " +\
            snakemake.input.tsv + " " +\
            snakemake.output.enrich + " " +\
            snakemake.output.gsea + " " +\
            snakemake.output.universe + " " +\
            snakemake.params.organism_go + " " +\
            snakemake.params.cutoff_log2fc_enrich + " " +\
            snakemake.params.cutoff_padj_enrich + " " +\
            snakemake.params.cutoff_log2fc_gsea + " " +\
            snakemake.params.cutoff_padj_gsea +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
