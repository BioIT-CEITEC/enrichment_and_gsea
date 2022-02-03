#############################################################
# wrapper for rule: final_report
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: final_report \n##\n")
f.close()

command = " cp '"+os.path.abspath(os.path.dirname(__file__))+"/enrichment_GSEA_final_report.Rmd' enrichment_gsea/"

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = """ Rscript -e "rmarkdown::render('enrichment_gsea/enrichment_GSEA_final_report.Rmd', params=list(config = '""" + snakemake.params.config + """' ))" """ +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " ls " + snakemake.output.html

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
