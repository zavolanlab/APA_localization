#!/bin/bash
#SBATCH --job-name=APEX_fracSanity
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=APEX_fracSanity_%j.log

# fracSanity, tandem PAS usage and APA of the APEX-seq data (Fazal et al. 2019) vs the unlabelled controls
# (WF/scripts/APEX_fracSanity_pipeline.py), then execution of APEX_seq_vs_Control.ipynb, ONT_organelles_and_APEX.ReadMapping.ipynb,
# APEX_seq.fracSanity.ipynb and APEX_vs_ONT_organelles.fracSanity.ipynb. Needs the outputs of the bulk workflow run APEX_seq_Fazal2019
# (WF/scripts/run_bulk_workflow.sh). Steps whose outputs exist are skipped, so the job can be resubmitted after a failure.
# Optional: STEPS=<pipeline steps separated by +> (sbatch --export splits on commas); NO_NOTEBOOK=1 to skip the notebooks.

source ~/miniforge3_Dec2024/etc/profile.d/conda.sh
conda activate sanity_utils
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=6

cd ~/Projects/APA_localization
set -a; source APA_localization.scicore.env; set +a
OUT=${WF_RUNS_DIR}Sanity_utils_APEX_seq_Fazal2019

python WF/scripts/APEX_fracSanity_pipeline.py --outdir "$OUT" --threads "$SLURM_CPUS_PER_TASK" ${STEPS:+--steps "${STEPS//+/,}"} || exit 1
[[ -n "${NO_NOTEBOOK:-}" ]] && exit 0
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8   # gene-set fits (SVD); more threads were slower
for nb in APEX_seq_vs_Control.ipynb ONT_organelles_and_APEX.ReadMapping.ipynb APEX_seq.fracSanity.ipynb APEX_vs_ONT_organelles.fracSanity.ipynb; do
    python WF/scripts/execute_notebook.py "$nb" || exit 1
done
