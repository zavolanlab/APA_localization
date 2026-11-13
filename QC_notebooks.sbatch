#!/bin/bash
#SBATCH --job-name=QC_notebooks
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=QC_notebooks_%j.log

# Executes notebooks in place (default: the two library-QC notebooks with the slope / compression analyses).
# Optional (sbatch --export=ALL,...): NOTEBOOKS=<notebooks separated by +>.

source ~/miniforge3_Dec2024/etc/profile.d/conda.sh
conda activate sanity_utils
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

cd ~/Projects/APA_localization
NOTEBOOKS=${NOTEBOOKS:-ONT_organelles_vs_Total.ipynb+APEX_seq_vs_Control.ipynb}
status=0
for nb in ${NOTEBOOKS//+/ }; do
    python WF/scripts/execute_notebook.py "$nb" || status=1
done
exit $status
