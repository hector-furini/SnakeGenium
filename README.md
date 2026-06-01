# **SnakeGenium v1.01**
## Integrated GWAS Pipeline for Bacteria

* Read the [Latest Changelog](CHANGELOG.md)

### **Requirements:**
* #### Linux64 <sub>(Tested, might work on other unix systems)</sub>
* #### Anaconda environment manager
* #### Phenotypes Table (.tsv)
* #### Annotated genome files ".fna", ".gff" and ".gbk" for both reference genome and target genomes

---
## About the Phenotypes Table:

Mandatory phenotypes table separated by tabulation (.tsv) <ins>__must__</ins> follow this organization:
- __(1)__ The First column must be the Genome Identifier Number (unique per genome).
- __(2)__ Other columns must contain the resistant/suceptible phenotype by antimicrobial using a binary presence\abscence identifier. Using 1 for resistant phenotype, 0 for susceptible and -1 for empty values as empty values can't be processed by the pipeline.
- __(3)__ Avoid any spaces in columns names or genome identifiers.
- __(4)__ The first column MUST be the genome identifier, however, you can add as much columns for antimicrobial phenotypes as you like beyond the first column.
- __(5)__ At least 10 phenotypes for both resistant and susceptible are required for each antimicrobial to run the analysis. Otherwise the pipeline might fail or crash.

### **EXAMPLE:**

| Genome_Identifier  | Antimicrobial_01_Phenotype | Antimicrobial_02_Phenotype | Antimicrobial_03_Phenotype |
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
- (3) Each of the genomes identifier subdirectory must include inside already annotated genome files such as the genome file in ".fna" format and the annotated genomic data ".gff" and ".gbk" files.
- (4) The ".fna", ".gff" and ".gbk" files must be named equally to its respective genome identifier.

### **EXAMPLE:**
<pre>
/annotation_directory
        ├── GENOME_01/
        |       ├── /GENOME_01.fna
        |       ├── /GENOME_01.gbk
        |       └── /GENOME_01.gff
        ├── GENOME_02/
        |       ├── /GENOME_02.fna
        |       ├── /GENOME_02.gbk
        |       └── /GENOME_02.gff
        ├── GENOME_03/
        |       ├── /GENOME_03.fna
        |       ├── /GENOME_03.gbk
        |       └── /GENOME_03.gff
        └── GENOME_04/
                ├── /GENOME_04.fna
                ├── /GENOME_04.gbk
                └── /GENOME_04.gff
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
as we strongly recommend, install using the 'snakemake.yml' file in the 'envs' directory. This ensures correct versioning for the packages.<br>
**-- Setting conda channel priority for flexible or disabled might be necessary --**

`conda env create -f envs/snakemake.yml`

5. Activate the conda environment using:

`conda activate snakemake`

6. To run the pipeline, start the script run_gwas.py using the following command. Change parameters and input files as needed.
**Running the pipeline may overwrite existing files in base directory "./output"**

```
python run_gwas.py --anno /annotations_directory --pheno phenotypes.tsv --ref_fna reference_genome.fna --ref_gff reference_genome.gff --ref_gff reference_genome.gbk --max-threads 32 --max-mem 30000(mb) --jobs 1
```

#### Arguments description:
``` --anno: Annotations directory structured as shown in the examples ```

``` --pheno: Tab Separated Table with Phenotypes ```

``` --ref_fna: Reference Genome FASTA or .FNA  ```

``` --ref_gff: Reference Genome Annotated .GFF ```

``` --ref_gbk: Reference Genome Annotated .GBK ```

``` --max_threads: Max threads to be used in the analisys. Default: [8]```

``` --max-mem: Max memory in megabytes to be used in the analisys. Default: [8000] ```

``` --jobs: Number of simultaneos jobs to be run by the pipeline. Default: [1] ```

**ATTENTION:** The first time the pipeline is run will take several minuters at its beginning due to dependencies and packages installation. 

5. Output will be generated inside the root folder in the `./output` directory.

<pre>
/output
      └── gwas/
            ├── phenotypes/                                     # Phenotypes tables separated by antimicrobial
                    ├── antimicrobial_1_name.tsv
                    └── antimicrobial_2_name.tsv
            ├── pyseer_results/                                 # Pyseer runs results
                    ├── genes_results/
                    ├── snps_results/
                    ├── summary/                                # All of Pyseers results summarized
                    └── unitigs_results/
            ├── roary/                                          # Roary matrix results
                    └── gene_presence_absecence.csv             # Gene presence/abscence matrix
            ├── samples/                                        # Samples directories list separated by antimicrobial
                    ├── antimicrobial_1_sample_list.txt
                    └── antimicrobial_2_sample_list.txt
            ├── snippy/
                    ├── core.aln                                # Core aligment file
                    └── core.vcf                                # Variants file
            ├── tree/
                    ├── core_genome.treefile                    # Phylogenetic tree
                    └── kinship_matrix                          # Kinship matrix for phylogenetic distances
            ├── unitig/
                    └── unitigs.pyseer.gz                       # Unitigs generated from sampled genomes
            └── fasta_list.txt                                  # Genomes Directory List
</pre>

## Please Cite

If you have used SnakeGenium for your results, please cite the tool and any other packages and datasets you may have used.

If a reference needs to be updated please let me know!

## References:
FastTree v2.1.11: [github.io/fasttree](https://morgannprice.github.io/fasttree/)
<br>Pyseer v1.4.0: [github.io/pyseer](https://github.com/mgalardini/pyseer)
<br>└── Pyseer Scripts: [github.io/pyseer/scripts](https://github.com/mgalardini/pyseer/tree/master/scripts)
<br>Roary v3.13.0: [github.io/roary](https://github.com/sanger-pathogens/Roary)
<br>Snakemake v9.9.0: [github.io/snakemake](https://github.com/snakemake/snakemake)
<br>Snippy v4.6.0: [github.io/snippy](https://github.com/tseemann/snippy)
<br>Unitig-Caller v1.3.0: [github.io/unitig-caller](https://github.com/bacpop/unitig-caller)

## Name Etimology:
SNAKE - **SNAKE**make/Python origin<br>
GEN - **GEN**ome Analysis<br>
IUM - Bacter**IUM**

## License:
**SnakeGenium** is free software, released under the [GPL (version 3)](https://www.gnu.org/licenses/gpl-3.0).

## Author

Hector H. Furini
