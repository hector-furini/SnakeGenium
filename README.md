# Integrated-GWAS-Pipeline-for-Bacteria

### **Requirements:** 
* #### Anaconda environment manager
* #### Phenotypes Table (.tsv)
* #### Annotated genome files ".fna" and ".gff" for both reference genome and target genomes

---
## About the Phenotypes Table:

Mandatory phenotypes table separated by tabulation (.tsv) <ins>__must__</ins> follow this organization:
- __(1)__ Column must be the Genome Identifier Number (unique per genome).
- __(2)__ Other columns must contain the resistant/suceptible phenotype by antimicrobial using a binary presence\abscence identifier. Using 1 for resistant phenotype, 0 for suceptible and -1 for empty values as empty values can't be processed by the pipeline.
- __(3)__ Avoid any spaces in columns names or genome identifiers.
- __(4)__ The first column MUST be the genome identifier, however, you can add as much columns for antimicrobial phenotypes as you like beyond the first column.

### **EXAMPLE:**

| Genome Identifier  | Antimicrobial_01_Phenotype | Antimicrobial_02_Phenotype | Antimicrobial_03_Phenotype |
| ------------- | ------------- |  ------------- | ------------- |
| GENOME_01  | 0  | 1 | 0 |
| GENOME_02  | 1  | 1 | 0 |
| GENOME_03  | -1  | 0 | 1 |
| GENOME_04  | 1  | -1 | 0 |

---

## About the Genomes Directories structuring:

Mandatory genomes annotation directory must be structure as it follows:
- (1) The main directory need to be parsed as a argument when running the script.
- (2) The main directory must include subdirectories named <ins>__EQUALLY__</ins> to the genome identifiers found on the phenotypes table.
- (3) Each of the genomes identifier subdirectory must include inside already annotated genome files such as the genome file in ".fna" format and the annotated genomic data ".gff" file.
- (4) Its recommended that the ".fna" and ".gff" files are named equally to the genome identifier.

### **EXAMPLE:**
<pre>
/annotation_directory
        ├── GENOME_01/
        |       ├── /GENOME_01.fna
        |       └── /GENOME_01.gff
        ├── GENOME_02/
        |       ├── /GENOME_02.fna
        |       └── /GENOME_02.gff
        ├── GENOME_03/
        |       ├── /GENOME_03.fna
        |       └── /GENOME_03.gff
        └── GENOME_03/
                ├── /GENOME_03.fna
                └── /GENOME_03.gff
</pre>

## Getting Started:

1. Install Anaconda from its official site: [Anaconda Download](https://www.anaconda.com/download)

2. Activate Conda base environment after installation on your preferred terminal:

`conda activate`

3. If you haven't configure your conda installation to use bioconda and conda-forge channels, from your conda environment, run:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

4. Clone this repository, and install Snakemake using anaconda package manager. You can install snakemake either from the repositories or,
as we strongly recommend, install using the 'snakemake.yml' file in the 'envs' directory. This ensures correct versioning for the packages.
<br>-- Setting conda channel priority for flexible or disabled might be necessary --

`conda env create -f envs/snakemake.yml`

5. Activate the conda environment using:

`conda activate snakemake`

6. To run the pipeline, start the script run_gwas.py using the following command. Change parameters and input files as needed:

```
python run_gwas.py --anno /annotations_directory --pheno phenotypes.tsv --ref_fna reference_genome.fna --ref_gff reference_genome.gff --max-threads 32 --max-mem 30000(mb)'`
```

**ATTENTION:** The first time the pipeline is run will take several minuters at its beginning due to dependencies and packages installation. 

5. Output will be generated inside the root folder in the `./output` directory.

## References:
FastTree v2.1.11: [github.io/fasttree](https://morgannprice.github.io/fasttree/)
<br>Pyseer v1.4.0: [github.io/pyseer](https://github.com/mgalardini/pyseer)
<br>└── Pyseer Scripts: [github.io/pyseer/scripts](https://github.com/mgalardini/pyseer/tree/master/scripts)
<br>Roary v3.13.0: [github.io/roary](https://github.com/sanger-pathogens/Roary)
<br>Snakemake v9.9.0: [github.io/snakemake](https://github.com/snakemake/snakemake)
<br>Snippy v4.6.0: [github.io/snippy](https://github.com/tseemann/snippy)
<br>Unitig-Caller v1.3.0: [github.io/unitig-caller](https://github.com/bacpop/unitig-caller)