
# Run Snakemake with the specified Snakefile and resources
snakemake -s gwas_snakefile --use-conda --cores 92 --resources mem_mb=300000 -j 2 -k --latency-wait 5 --rerun-incomplete --rerun-triggers mtime #--summary > pipeline_benchmark.txt

# --snakefile or -s: specifies the Snakefile to use
# --forcerun <rule_name> or -R <rule_name>: trigger the execution of a specific rule
# --forceall or -F: trigger the execution of all rules (BE CAREFUL)
# --use-conda: enables the use of conda environments for rules
# --cores: specifies the number of cores to use
# --resources: specifies the resources available for the workflow
# -j: specifies the maximum number of jobs to run simultaneously
# -k: keeps going even if some jobs fail
# --keep-going: allows the workflow to continue running other jobs even if some fail
# --keep-incomplete: keeps incomplete jobs in the output directory
# --latency-wait: sets the time to wait for files to appear before retrying
# --rerun-incomplete: reruns incomplete jobs
# --rerun-triggers: specifies triggers for rerunning jobs based on modification times of input and output files
# --summary: generates a summary of the workflow execution
# > pipeline_benchmark.txt: redirects the output to a file for benchmarking