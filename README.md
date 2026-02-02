# ❄️ Lake Untersee BGC mining

This rep contains code used for the analysis described in:

**"Long-read metagenomic sequencing recovers diverse biosynthetic gene clusters in Antarctic microbial mats from Lake Untersee, Dronning Maud Land"**

The workflow enables processing of long-read metagenomic data and recovery of biosynthetic gene clusters (BGCs) and metagenome-assembled genomes (MAGs).


## 🧬 Data Availability

The following datasets from the study are publicly available:

- **Raw sequencing data (FASTQ)**: [NCBI SRA, accession pending]
- **Metagenome assemblies**: [NCBI GenBank, accession pending]
- **High-quality (HQ) and medium-quality (MQ) MAGs**: [NCBI GenBank, accession pending]


## 🚀 Running the Analysis

### Requirements

Make sure you have the following tools installed:

- `conda` ≥ *[insert version]*  
- `snakemake` ≥ *[insert version]*  

These are used to manage the workflow and software environments. All tools required by the analysis are handled automatically by Snakemake through conda environment yaml files specified within the workflow (in envs folder).

A database for Kraken2 is required to run this analysis and is available [here](https://benlangmead.github.io/aws-indexes/k2) - the version used in this analysis is the December 2022 PPF. The AntiSMASH database is also required to run this analysis. This can be downloaded here following the instructions on the [antismash website](https://docs.antismash.secondarymetabolites.org/install/). 

### Usage

1. **Clone this repository**:

   ```bash
   git clone https://github.com/<your-org>/untersee-bgc.git
   cd untersee-bgc
   ```

2. **Download the data from NCBI SRA**

3. **Run the analysis**:

    ```bash
    snakemake --use-conda --cores <number_of_threads>
    ```
This command will execute the full workflow using Snakemake and create all necessary conda environments.
