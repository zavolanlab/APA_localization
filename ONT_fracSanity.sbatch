#!/bin/bash
#SBATCH --job-name=ONT_fracSanity
#SBATCH --cpus-per-task=48
#SBATCH --mem=160G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=ONT_fracSanity_%j.log

# fracSanity, cleavage-site usage and APA of the ONT subcellular fractions (siCtrl) vs Total
# (WF/scripts/ONT_fracSanity_pipeline.py), then execution of ONT_organelles.fracSanity.ipynb.
# Steps whose outputs exist are skipped, so the job can be resubmitted after a failure.
# Optional: STEPS=<comma-separated pipeline steps> to run only these; NO_NOTEBOOK=1 to skip the notebook
#   e.g. sbatch --export=ALL,NO_NOTEBOOK=1 ONT_fracSanity.sbatch

source ~/miniforge3_Dec2024/etc/profile.d/conda.sh
conda activate sanity_utils
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=6

cd ~/Projects/APA_localization
set -a; source APA_localization.scicore.env; set +a
OUT=${WF_RUNS_DIR}Sanity_utils_ONT_organelles_siCtrl

python WF/scripts/ONT_fracSanity_pipeline.py --outdir "$OUT" --threads "$SLURM_CPUS_PER_TASK" ${STEPS:+--steps "$STEPS"} || exit 1
[[ -n "${NO_NOTEBOOK:-}" ]] && exit 0
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8   # gene-set fits (SVD); more threads were slower
python WF/scripts/execute_notebook.py ONT_organelles.fracSanity.ipynb
