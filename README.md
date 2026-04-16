[![DOI](https://zenodo.org/badge/1203660686.svg)](https://doi.org/10.5281/zenodo.19482094)
# APA localization - Analysis and Pipelines

This repository contains the computational workflows and downstream analysis notebooks related to the analysis of alternative polyadenylation isoforms in subcellular compartments.

The repository is optimized for running BOTH the workflows and analysis in jupyter notebook on HPC cluster.

On sciCORE HPC, running jupyter notebook on a computational node is nicely enabled by [OnDemand service](https://docs.scicore.unibas.ch/HPC%20Cluster/interactivecomputing/#open-ondemand-ood).

We utilize a hybrid approach: **Snakemake** for robust, scalable data processing on HPC clusters (sciCORE), and **Jupyter Notebooks** for interactive downstream analysis and visualization.

# Current state
Currently, we've reanalyzed the bulk RNA-seq data from the study [System-wide analysis of RNA and protein subcellular localization dynamics](https://www.nature.com/articles/s41592-023-02101-9)
We've run basic processing, alignment, and quantification of gene expression using customly prepared .gtf file and [FeatureCounts](https://subread.sourceforge.net/featureCounts.html) utility.
We further quantified relative usage of polyadenylation sites (PASs) with a [PAQR2 workflow](https://link.springer.com/article/10.1186/s13059-018-1415-3). 
At that, we've used the modified version of the workflow from the paper [Leveraging multi-omics data to infer regulators of mRNA 3’ end processing in glioblastoma](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2024.1363933/full).
As input for PAQR2, we used human [PolyASite Atlas v3.0](https://academic.oup.com/nar/article/53/D1/D197/7893321) filtered for 62% stringency level (see the paper for explanation of the optimality of that particular threshold value).

We further focused on **NFYA** gene and alternative polyadenylation at its terminal exon, in a collaborative project with the [group of Prof. Dr. Paolo Gandellini](https://www.unimi.it/it/ugov/person/paolo-gandellini).

We are also working on **nanopore data processing**.

**Action plan** 
In general, different sub-projects related to subcellular localization of APA isoforms, will be corresponding to different jupyter notebooks.

For now, NFYA project code and data were finalized. ONT nanopore data processing and analysis is ongoing.

## Repository Structure

```text
.
├── NFYA_project.ipynb                  # a Jupyter notebook dedicated to NFYA project, includes analysis and workflow configuration
├── ONT_analysis.ipynb                  # a Jupyter notebook dedicated to nanopore data processing and analysis
├── APA_localization.template.env       # Template for required environment variables/paths
└── WF/                                 # Snakemake Workflow Engine
    ├── Snakefile-prepare-faster        # Pipeline Step 1: RNA-seq data processing (alignment, FastQC)
    ├── Snakefile-quantification-faster # Pipeline Step 2: Quantification of gene expression with FeatureCounts, separating .bam files by chromosomes for efficiency, preparation of coverages for PAQR quantification
    ├── Snakefile-PAQR-quantify         # Pipeline Step 3: Running final PAQR quantification to obtain PAS-vs-sample count matrix
    ├── config.template.yaml            # Template configuration for Snakemake parameters
    ├── envs/                           # Conda environments isolated for specific Snakemake rules
    ├── profile/                        # SLURM execution profile for the HPC
    └── scripts/                        # Python and R scripts utilized by both Snakemake and Jupyter
```

## Quick Start & Setup

To ensure strict reproducibility and security, this project uses `.env` files to manage all absolute paths (data directories, genome annotations, etc.). **Do not hardcode paths into the Python or Snakemake files.**

### 1. Clone the Repository
Clone this repository into your local user space (`$HOME`):
```bash
git clone https://github.com/zavolanlab/APA_localization.git
cd APA_localization
```

### 2. Configure Environment Paths
You must map the project to your local HPC paths. 
**First**, copy the template, rename it, and fill in your absolute paths, for example like that:
  ```bash
  cp APA_localization.template.env APA_localization.scicore.env
  # Open .env and edit the "Base Directories" section to match your system
  ```
* **Recommended if you are a Zavolan group member on sciCORE:** move the `APA_localization.scicore.env` to Project GROUP folder and symlink into your local repository directory:
  ```bash
  ln -s <a file with specified sciCORE paths> APA_localization.scicore.env
  ```
This way `APA_localization.scicore.env` will be automatically accessible by group members but will not be tracked by git.
*(Note: `*.env` files are ignored by git to protect private cluster paths, except the `APA_localization.template.env` file).*
**(`APA_localization.scicore.env` does exist in the GROUP folder of the Project on Scicore. Look for README there.)

### 3. Install the conda environment with zavolab_pyutils
Analysis in the notebook is largely based on the functions from [zavolab_pyutils](https://github.com/zavolanlab/zavolab_pyutils/tree/dev) repository.
Follow the instruction from that repo "Developer Setup from source, with conda environment".
Use the created conda environment "zavolab_pyutils" to execute the Jupyter Notebook.

### 4. Essential for developpers! Install nbstripout
When in the `APA_localization` directory, run:
```bash
nbstripout --install
```
This will automatically hide the output of cells in juputer notebooks when pushed to github! Otherwise there is a risk of exposing your HPC cluster paths to public.

### 5. Use the juputer notebook to configure the workflow and input table preparation
Configuration of the workflows (i.e. creation of input .tsv with sample specification and .yaml config is done **inside** the jupyter notebook)

### 6. Executing the Workflows

The heavy lifting is divided into (currently, three) separate Snakemake workflows located in the `WF/` directory.

`Bash` commands are also prepared inside the jupyter notebook. They should be further copied into command line and executed.

**On an HPC cluster like sciCORE**, workflows should be executed on a **login** node. Snakemake further automatically submits jobs to computational nodes.

### 7. Downstream Analysis
Once the Snakemake workflows are complete, all results are routed to the shared group directories defined in your `.env` file. 

Use respective sections of the Jupyter Notebook to analyze the outputs. 

The notebook automatically loads your `.env` paths using `python-dotenv`, allowing it to dynamically locate all workflow results, figures, and metadata regardless of where you cloned this repository.
