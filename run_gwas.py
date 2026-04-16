import argparse
from subprocess import run
import os

def main():
    parser = argparse.ArgumentParser(description="Starting GWAS pipeline script.",
                                    epilog="""
                                    This script requires:
                                    1. Annotation files:
                                    Genome '.fna',
                                    Annotated '.gff'.
                                    2. Phenotypes table '.tsv".
                                    3. Reference Genome annotation:
                                    Reference '.fna' and '.gff'.
                                    For more reference go to 'https://github.com/hector-furini/Integrated-GWAS-Pipeline-for-Bacteria'
                                    """)


    # POSITIONAL ARGUMENT
    parser.add_argument("--anno",
                        type=str,
                        help="Path to annotation genome files directory. Attention: In this directory must exist subdirectories for each genome.",
                        required=True,
                        metavar="['ANNOT_DIR']")

    parser.add_argument("--pheno",
                        type=str,
                        help="Tab separated phenotypes table (.tsv) path",
                        required=True,
                        metavar="['TSV_TABLE']")

    parser.add_argument("--ref_fna",
                        type=str,
                        help="Reference genome .fna file path",
                        required=True,
                        metavar="['REF_FNA']")

    parser.add_argument("--ref_gff",
                        type=str,
                        help="Reference genome .gff file path",
                        required=True,
                        metavar="['REF_GFF']")
    
    # OPTINAL ARGUMENT
    parser.add_argument("--max_threads",
                        default=8,
                        type=int,
                        help="Number of threads to be used. Default: 8")

    parser.add_argument("--max_mem",
                        default=8000,
                        type=int,
                        help="Max amount of memory to be used by the pipeline in megabytes. Default: 8000")

    parser.add_argument("--jobs",
                        default=1,
                        type=int,
                        help="Max simultaneos jobs to run. Default: 1")

    args = parser.parse_args()
    writeConfigFile(checkFiles(args))
    runSnakemake(args.max_threads, args.max_mem, args.jobs)

def checkFiles(args):
    for key, value in vars(args).items():
        if key in ["anno"]:
            if not os.path.exists(value):
                print(f"ERROR: --{key} not a valid dir. Given: {value}")
                return
        elif key in ["pheno", "ref_fna", "ref_gff"]:
            if not os.path.isfile(value) or not value.rsplit(".", 1)[-1] in ["tsv", "fna", "gff"]:
                print(f"ERROR: --{key} not a valid file. Given: {value}")
                return
        elif key in ["max_threads", "max_mem", "jobs"]:
            if value <= 0:
                print(f"--{key} value must be bigger than 0. Given value: {value}")
                return
    print(f"All files verified sucessfully!")
    return args

def writeConfigFile(args):
    print("Writing to Config File...")
    with open("config.yaml", "w") as config_file:
        config_file.write(base_config_yaml)
        config_file.write(
f"""
# -- Input Parameters --
prokka_out: {args.anno}
ref_fna: {args.ref_fna}
ref_gff: {args.ref_gff}
phen_table: {args.pheno}
max_threads: {args.max_threads}
""")

def runSnakemake(threads, mem, jobs):

    snakemake_summary_run = [f"snakemake", "-s", "gwas_snakefile", "-n", "--cores",
    f"{threads}", "--resources", f"mem_mb={mem}", "--rerun-triggers", "mtime"]

    snakemake_run = ["snakemake", "-s", "gwas_snakefile", "--use-conda", "--cores",
    f"{threads}", "--resources", f"mem_mb={mem}", "-j", f"{jobs}", "-k", "--latency-wait", "5",
    "--rerun-incomplete", "--rerun-triggers", "mtime"]

    print("Generating pipeline summary file at root directory")
    with open("GWAS_summary.txt", "w") as f:
        run(snakemake_summary_run, stdout=f)

    print(f"Running snakemake with {jobs} simultaneos jobs. Max {threads} threads and {mem}mb memory.")
    run(snakemake_run)

base_config_yaml ="""# -- Root Directories --
pyseer_out: "output/gwas/pyseer"
gwas_tree: "output/gwas/tree"
gwas_data: "output/gwas"
gwas_snippy: "output/gwas/snippy"
gwas_unitig: "output/gwas/unitig"
gwas_pyseer: "output/gwas/pyseer_results"

# -- Scripts --
r_effect_bubble : "scripts/plot_effect_bubble.R"
gwas_summarize : "scripts/summarize_gwas_results.py"
gwas_filter_hits : "scripts/filter_significant_hits.py"
gwas_phandango : "scripts/phandango_plot.py"
pyseer_phylogeny : "scripts/pyseer_scripts/phylogeny_distance.py"
pyseer_qqplot: "scripts/pyseer_scripts/qq_plot.py"
pyseer_patterns: "scripts/pyseer_scripts/count_patterns.py"

# -- Conda enviroments --
conda_fasttree: "envs/fasttree.yml"
conda_unitig: "envs/unitig_caller.yml"
conda_pyseer: "envs/pyseer.yml"
conda_snippy: "envs/snippy.yml"
conda_roary: "envs/roary.yml"
conda_r: "envs/r_env.yml"
"""

if __name__ == "__main__":
    main()