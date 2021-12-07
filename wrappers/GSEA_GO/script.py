#############################################################
# wrapper for rule: GSEA_GO
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: GSEA_GO \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/GSEA_GO.R "+\
            snakemake.params.workdir + " " +\
            snakemake.input.tsv + " " +\
            snakemake.params.outdir + " " +\
            snakemake.params.organism_go + " " +\
            snakemake.params.cutoff_log2fc + " " +\
            snakemake.params.cutoff_padj + " " +\
            snakemake.params.n_up + " " +\
            snakemake.params.n_down + " " +\
            snakemake.params.colors + " " +\
            snakemake.params.gsea_padj + " " +\
            snakemake.params.gsea_padjmethod + " " +\
            snakemake.params.gsea_minGSSize + " " +\
            snakemake.params.gsea_maxGSSize + " " +\
            snakemake.params.gsea_eps + " " +\
            snakemake.params.gsea_nPermSimple + " " +\
            snakemake.params.gsea_by +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
