#!/bin/bash
#SBATCH --job-name=LoRNA_fracSanity
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=LoRNA_fracSanity_%j.log

# fracSanity on the U2OS density-gradient fractionation (LoRNA, Villanueva et al. 2024) and its test against the SIRV
# spike-ins (WF/scripts/LoRNA_fracSanity_pipeline.py), then execution of LoRNA_U2OS.fracSanity.ipynb. Needs the bulk
# workflow run U2OS_subcellular_fractions and, for the steps spike_ins and fit_cohorts_spikeprior, the SIRV counts of
# WF/scripts/count_spike_ins.sbatch. Steps whose outputs exist are skipped, so the job can be resubmitted after a failure.
# Optional (sbatch --export=ALL,...): STEPS=<pipeline steps separated by + or ,>; NO_NOTEBOOK=1 to skip the notebook.

source ~/miniforge3_Dec2024/etc/profile.d/conda.sh
conda activate sanity_utils
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=6

cd ~/Projects/APA_localization
set -a; source APA_localization.scicore.env; set +a
OUT=${WF_RUNS_DIR}Sanity_utils_LoRNA_U2OS

python WF/scripts/LoRNA_fracSanity_pipeline.py --outdir "$OUT" --threads "$SLURM_CPUS_PER_TASK" ${STEPS:+--steps "${STEPS//+/,}"} || exit 1
[[ -n "${NO_NOTEBOOK:-}" ]] && exit 0
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
python WF/scripts/execute_notebook.py LoRNA_U2OS.fracSanity.ipynb || exit 1
