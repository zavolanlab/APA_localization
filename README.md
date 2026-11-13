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
├── NFYA_project.ipynb                  # NFYA project: analysis and workflow configuration
├── ONT_analysis.ipynb                  # nanopore (ONT) data: custom genome with ERCC spike-ins, input tables and configuration of the nanoflowz runs, bulk RNA-seq of the same samples and its comparison with ONT
├── ONT_organelles_vs_Total.ipynb       # ONT subcellular fractions (Free, Org2, ER, TIS, Golgi, PEX, End/TJ) vs Total: read QC, nuclear contamination, log2 CPM slopes, library complexity, input-amount control (exported as an HTML report)
├── ONT_organelles.fracSanity.ipynb     # ONT subcellular fractions vs Total (siCtrl): fracSanity (alpha_g, F), ISMARA gene-set over-/under-representation, cleavage-site usage (one CS per read 3' end position), tandem APA in PAQR3 terminal exons, alpha of cleavage sites (exported as an HTML report)
├── ONT_fracSanity.sbatch               # SLURM job: WF/scripts/ONT_fracSanity_pipeline.py, then execution of ONT_organelles.fracSanity.ipynb
├── APEX_seq_vs_Control.ipynb           # public APEX-seq (Fazal et al. 2019, GSE116008, HEK293T) vs unlabelled controls: read QC, nuclear RNA, log2 CPM slopes, library complexity, APEX vs fractionation (exported as an HTML report)
├── APEX_seq.fracSanity.ipynb           # APEX-seq vs controls: fracSanity, ISMARA gene sets, tandem PAS usage (PAQR), APA, alpha of PAS (exported as an HTML report)
├── APEX_vs_ONT_organelles.fracSanity.ipynb # comparison of the fracSanity patterns (r_g, gene-set activities, APA) of the ONT fractions and the APEX-seq constructs (exported as an HTML report)
├── ONT_organelles_and_APEX.ReadMapping.ipynb # read numbers and gene-class / gene-biotype content of all ONT libraries (fractions and Totals) and APEX-seq libraries, biotype enrichment vs the references (as pbody_RNA_Spang.ReadMapping.ipynb; exported as an HTML report)
├── APEX_fracSanity.sbatch              # SLURM job: WF/scripts/APEX_fracSanity_pipeline.py, then execution of the three APEX notebooks
├── LoRNA_U2OS.fracSanity.ipynb         # U2OS density-gradient fractionation (LoRNA, Villanueva et al. 2024; run U2OS_subcellular_fractions): fracSanity tested against the SIRV spike-ins and the mass balance of the fractions, comparison with the ONT fractions and APEX-seq (exported as an HTML report)
├── LoRNA_fracSanity.sbatch             # SLURM job: WF/scripts/LoRNA_fracSanity_pipeline.py, then execution of LoRNA_U2OS.fracSanity.ipynb
├── CD47_analysis.ipynb                 # ONT data of CD47 pgRNA clones: nanoflowz configuration, gene- and cleavage-site-level expression analysis
├── APA_localization.template.env       # Template for required environment variables/paths
└── WF/                                 # Snakemake Workflow Engine
    ├── Snakefile-prepare-faster        # Pipeline Step 1: RNA-seq data processing (alignment, FastQC), paired-end reads
    ├── Snakefile-prepare-faster-se     # Pipeline Step 1 for single-end reads
    ├── Snakefile-quantification-faster # Pipeline Step 2: Quantification of gene expression with FeatureCounts, separating .bam files by chromosomes for efficiency, preparation of coverages for PAQR quantification
    ├── Snakefile-PAQR-quantify         # Pipeline Step 3: Running final PAQR quantification to obtain PAS-vs-sample count matrix
    ├── Snakefile-basic-TakaraPicoUMI-pe # paired-end Takara SMART-seq Pico UMI libraries: UMI extraction, adapter trimming, STAR mapping, UMI deduplication, FeatureCounts
    ├── config.yaml                     # Snakemake parameters (written by the notebooks)
    ├── cluster.json                    # default SLURM resources
    ├── jobstats-config.ini             # settings of the jobstats utility
    ├── envs/                           # Conda environments isolated for specific Snakemake rules
    ├── profile/                        # SLURM execution profile for the HPC
    └── scripts/                        # Python and R scripts utilized by both Snakemake and Jupyter
        ├── ONT_fracSanity_pipeline.py  # heavy steps of ONT_organelles.fracSanity.ipynb: counts and cleavage sites from the nanoflowz outputs, Sanity, fracSanity, CS usage and APA tests (Sanity_utils)
        ├── APEX_fracSanity_pipeline.py # heavy steps of the APEX-seq notebooks: counts and tandem PAS from the bulk workflow outputs, fragment duplicates, Sanity, fracSanity, PAS usage and APA tests
        ├── LoRNA_fracSanity_pipeline.py # heavy steps of LoRNA_U2OS.fracSanity.ipynb: counts, spike-in RNA content per fraction, mass balance, Sanity, fracSanity with three priors on F
        ├── count_spike_ins.sbatch      # SLURM array: SIRV / ERCC read pairs per library (STAR on the spike-in sequences in $EXTERNAL_DATA_DIR/spike_ins)
        ├── download_sra_fastq.sh       # SRA toolkit download (prefetch, fasterq-dump, pigz) of paired-end runs; run on the login / vscode node
        ├── run_bulk_workflow.sh        # runs Snakefile-prepare-faster, -quantification-faster and -PAQR-quantify for one config (as printed in NFYA_project.ipynb)
        ├── execute_notebook.py         # executes a notebook in place (e.g. in an sbatch job)
        └── make_html_reports.py        # HTML reports of the analysis notebooks without code (see section 8)
```

Nanopore data are processed with the Nextflow workflow [nanoflowz](https://github.com/zavolanlab/nanoflowz); its runs are configured in `ONT_analysis.ipynb` and `CD47_analysis.ipynb`.

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

### 3. Install the conda environments for the notebooks
The notebooks run in one of two conda environments. The environment of each notebook is stored in its kernel metadata:

| Environment | Notebooks |
|---|---|
| `zavolab_pyutils` | `NFYA_project.ipynb`, `ONT_analysis.ipynb`, `CD47_analysis.ipynb` |
| `sanity_utils` | `ONT_organelles_vs_Total.ipynb`, `ONT_organelles.fracSanity.ipynb` (with `ONT_fracSanity.sbatch`), the three APEX notebooks (with `APEX_fracSanity.sbatch`), `LoRNA_U2OS.fracSanity.ipynb` (with `LoRNA_fracSanity.sbatch`), HTML export (`WF/scripts/make_html_reports.py`) |

**`zavolab_pyutils`.** Most of the analysis is based on the functions from the [zavolab_pyutils](https://github.com/zavolanlab/zavolab_pyutils/tree/dev) repository.
Follow the instruction from that repo "Developer Setup from source, with conda environment".

**`sanity_utils`.** [Sanity_utils](https://github.com/nimwegenLab/Sanity_utils) is a private repository: cloning it needs an SSH key with access to `nimwegenLab`. Its environment also installs `zavolab_pyutils`, `ipykernel` and `python-dotenv`:
```bash
git clone git@github.com:nimwegenLab/Sanity_utils.git
cd Sanity_utils
PYTHONNOUSERSITE=1 conda env create -f environment.yml
conda activate sanity_utils
sanity-utils-install-sanity      # the Sanity binary; on sciCORE add --arch-flag=-march=x86-64-v2
```
`PYTHONNOUSERSITE=1` keeps packages from `pip install --user` (`~/.local`) from shadowing the environment's own.

**Jupyter kernels.** Both environments need `ipykernel` (`conda install -n zavolab_pyutils -c conda-forge ipykernel` if it is missing). On sciCORE OnDemand, conda environments with `ipykernel` show up as kernels automatically ("conda-env-…-py"). Elsewhere, register them once:
```bash
conda activate sanity_utils
python -m ipykernel install --user --name sanity_utils --display-name "Python (sanity_utils)" --env PYTHONNOUSERSITE 1
conda activate zavolab_pyutils
python -m ipykernel install --user --name zavolab_pyutils --display-name "Python (zavolab_pyutils)"
```
and select the kernel in Jupyter (Kernel > Change Kernel).

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

**Public data (APEX-seq, Fazal et al. 2019).** The runs are listed in `$METADATA_DIR/APEX_seq_Fazal2019.selected_runs.txt` (sample sheet: `start_samples.APEX_seq_Fazal2019.tsv`, annotation: `APEX_seq_Fazal2019.sample_annotation.tsv`) and downloaded on the login / vscode node (compute nodes have no internet), then processed with the three workflows (config `wf_runs/APEX_seq_Fazal2019/config.yaml`, as for the NFYA project but stranded PAQR coverages; the STAR index of `U2OS_subcellular_fractions` is linked):
```bash
nohup setsid WF/scripts/download_sra_fastq.sh $METADATA_DIR/APEX_seq_Fazal2019.selected_runs.txt $FASTQ_DIR/APEX_seq_Fazal2019 3 > $FASTQ_DIR/APEX_seq_Fazal2019/download.log 2>&1 &
nohup setsid WF/scripts/run_bulk_workflow.sh $WF_RUNS_DIR/APEX_seq_Fazal2019/config.yaml $FASTQ_DIR/APEX_seq_Fazal2019/download.log > $WF_RUNS_DIR/APEX_seq_Fazal2019/workflow_head.log 2>&1 &
sbatch APEX_fracSanity.sbatch    # afterwards: pipeline and the three APEX notebooks
```

**Density-gradient fractionation (LoRNA, U2OS; run `U2OS_subcellular_fractions` of `NFYA_project.ipynb`) with spike-ins.** The SIRV sequences and annotation (Lexogen SIRV-Set, release 170612a / Lot00141 annotation, from a public mirror; converted to per-SIRV coordinates) and ERCC92 are in `$EXTERNAL_DATA_DIR/spike_ins` (`SIRV_ERCC92.fa`, `SIRV_ERCC92.gtf`, `STAR_index/`):
```bash
S=$METADATA_DIR/start_samples.U2OS_subcellular_fractions.tsv
sbatch --array=1-$(($(wc -l < $S) - 1)) WF/scripts/count_spike_ins.sbatch $S $WF_RUNS_DIR/U2OS_subcellular_fractions_spike_ins
sbatch --dependency=afterok:<array job id> LoRNA_fracSanity.sbatch    # pipeline and LoRNA_U2OS.fracSanity.ipynb
```

### 7. Downstream Analysis
Once the Snakemake workflows are complete, all results are routed to the shared group directories defined in your `.env` file. 

Use respective sections of the Jupyter Notebook to analyze the outputs. 

The notebook automatically loads your `.env` paths using `python-dotenv`, allowing it to dynamically locate all workflow results, figures, and metadata regardless of where you cloned this repository.

### 8. HTML reports for collaborators
The notebooks listed in `NOTEBOOKS` of `WF/scripts/make_html_reports.py` (the ONT organelle and APEX-seq notebooks) double as reports: each saved figure is shown with its file name and links to the full-size PNG/PDF, and cells tagged `report-skip` (setup, configuration) are left out. After running a notebook (Kernel > Restart & Run All; for `ONT_organelles.fracSanity.ipynb`: `sbatch ONT_fracSanity.sbatch`, for the APEX notebooks: `sbatch APEX_fracSanity.sbatch`, which also runs its pipeline and needs the human gene sets in `$EXTERNAL_DATA_DIR/genesets_human`, downloaded once on a node with internet access with `sanity-utils-fetch-genesets all --organism human --dest <dir>`), export it without code from the project folder:
```bash
conda activate sanity_utils
python WF/scripts/make_html_reports.py --out <report_folder>                              # all notebooks listed in NOTEBOOKS of the script
python WF/scripts/make_html_reports.py --out <report_folder> ONT_organelles_vs_Total.ipynb  # selected notebooks only
```
The report folder contains one HTML file per notebook, an `index.html`, and `figures/` with the full-size files that the reports link to. Share the folder as a whole (zip, SWITCHdrive).
