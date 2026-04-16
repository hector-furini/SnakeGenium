# Integrated-GWAS-Pipeline-for-Bacteria

Requirements: Anaconda package manager. Annotated genome files ".fna" and ".gff" (Tested with prokka annotation).

1. If you haven't configure your conda installion to use bioconda and conda-forge channels, run:

'''
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
'''

2. Install the conda environment containing snakemake depencies using the file "./envs/snakemake.yml"
'conda env create -f environment.yml'

3. Activate the conda environment using
'conda activate snakemake'

4. To run the pipeline, start the script run_gwas.py using the following command. Substiture parameters and input files as needed:
'python run_gwas.py --anno annotations_directory/ --pheno phenotypes.tsv --ref_fna reference_genome.fna --ref_gff reference_genome.gff --max-threads 32 --max-mem 30000(mb)'

## References:
FastTree
Pyseer
Roary
Snippy
Unitig_caller
