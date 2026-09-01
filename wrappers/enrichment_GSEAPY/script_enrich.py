#############################################################
# wrapper for rule: enrichment_enrichr
#############################################################
import os
import subprocess
import sys
import time
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
command = "awk -F'\t' '$1 != $2' " + snakemake.input.tsv + " | cut -f 2 | tail -n +2 | sed -e 's/__.*$//' | sort | uniq | grep -v '^$' > " + snakemake.params.gene_list
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
              " -v -g " + snakemake.params.enrichr_db + \
              " -t " + snakemake.params.n_up + " " + \
              " --cut " + snakemake.params.enrich_padj
    f = open(log_filename, 'a+')
    f.write("## COMMAND: " + command + "\n")

    # Run gseapy, capturing stderr to a per-run temp file. The exit code is
    # captured via subprocess.run (no raise on non-zero) instead of being
    # masked with "|| true", so real failures are detected, not just the
    # "no enrich terms" ValueError.
    error_log = snakemake.params.outdir + "/gseapy_error.tmp"
    command_with_capture = command + " 2> " + error_log

    # Transient Enrichr API failures ("Error sending gene list", timeouts,
    # HTTP 5xx) are retried with backoff before giving up.
    transient_markers = ["Error sending gene list", "timed out", "Timeout",
                         "ConnectionError", "Max retries", "HTTPError",
                         "500", "502", "503", "504"]
    max_attempts = 3
    failed = False
    error_content = ""
    for attempt in range(1, max_attempts + 1):
        # shell.has_fail is not available in all Snakemake versions, so run
        # the command directly; stderr is already redirected to error_log by
        # the shell redirection in the command itself. check=False keeps the
        # exit code for inspection instead of raising.
        proc = subprocess.run(["/bin/bash", "-c", command_with_capture], check=False)
        failed = proc.returncode != 0
        error_content = ""
        if os.path.exists(error_log):
            with open(error_log, 'r') as errf:
                error_content = errf.read()
        f = open(log_filename, 'a+')
        if not failed:
            f.write("## gseapy enrichr succeeded on attempt %d\n" % attempt)
            f.close()
            break
        is_transient = any(m in error_content for m in transient_markers)
        if is_transient and attempt < max_attempts:
            wait = 30 * attempt
            f.write("## WARNING: transient Enrichr API failure on attempt %d/%d, retrying in %ds:\n%s\n"
                    % (attempt, max_attempts, wait, error_content))
            f.close()
            time.sleep(wait)
            continue
        f.write("## gseapy enrichr failed on attempt %d/%d\n" % (attempt, max_attempts))
        f.close()

    # "No enrich terms" after cutoff -> expected empty-result case
    no_terms = ('ValueError' in error_content or 'No enrich terms' in error_content)

    if failed and no_terms:
        f = open(log_filename, 'a+')
        f.write("## WARNING: gseapy enrichr found no enrich terms after cutoff. Creating empty outputs.\n")
        # Create empty plot file
        shell("touch " + snakemake.output.plot)
        # Create empty TAB-delimited output file with header
        header = "Gene_set\tTerm\tOverlap\tP-value\tAdjusted P-value\tOld P-value\tOld Adjusted P-value\tOdds Ratio\tCombined Score\tGenes\n"
        with open(snakemake.output.table, 'w') as outf:
            outf.write(header)
        f.write("## Created empty output files.\n")
        f.close()
        shell("rm -f " + error_log)
    elif failed:
        # Genuine failure (API error after retries, bad parameters, crash):
        # keep the traceback in the log and fail loudly so Snakemake reports
        # the real cause instead of a cryptic missing-output error.
        f = open(log_filename, 'a+')
        f.write("## ERROR: gseapy enrichr failed:\n" + error_content + "\n")
        f.close()
        shell("cat " + snakemake.params.outdir + "/gseapy.enrichr.*.log >> " + log_filename + " 2> /dev/null || true")
        shell("rm -f " + snakemake.params.outdir + "/gseapy.enrichr.*.log")
        shell("rm -f " + error_log)
        sys.exit(1)
    else:
        # Success: append logs and clean up
        shell("cat " + snakemake.params.outdir + "/gseapy.enrichr.*.log >> " + log_filename + " 2> /dev/null || true")
        shell("cat " + error_log + " >> " + log_filename + " 2> /dev/null || true")
        shell("rm -f " + snakemake.params.outdir + "/gseapy.enrichr.*.log")
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