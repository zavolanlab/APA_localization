#!/bin/bash
# Runs the bulk RNA-seq workflow (Snakefile-prepare-faster, -quantification-faster, -PAQR-quantify) for one config,
# as printed in NFYA_project.ipynb; snakemake submits the jobs to SLURM. Run from the login / vscode node
# (detached, e.g. nohup setsid ... &), from the project folder, so that the conda envs in .snakemake/ are reused.
#
# Usage: WF/scripts/run_bulk_workflow.sh <config.yaml> [wait-for file] [steps, default "prepare-faster quantification-faster PAQR-quantify"]
#   wait-for file: start only once this file contains "ALL FINISHED" (e.g. the log of download_sra_fastq.sh)

set -uo pipefail
CONFIG=$1
WAIT_FOR=${2:-}
STEPS=${3:-"prepare-faster quantification-faster PAQR-quantify"}
PROJECT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$PROJECT"
set -a; source APA_localization.scicore.env; set +a
source ~/miniforge3_Dec2024/etc/profile.d/conda.sh
conda activate snakemake

if [[ -n "$WAIT_FOR" ]]; then
    until grep -q "ALL FINISHED" "$WAIT_FOR" 2>/dev/null; do sleep 120; done
    if grep -q "FAILED" "$WAIT_FOR"; then echo "$(date) downloads failed, see $WAIT_FOR"; exit 1; fi
fi

for step in $STEPS; do
    echo "$(date) === $step"
    snakemake \
        --snakefile "${WF_DIR}Snakefile-$step" \
        --scheduler greedy \
        --configfile "$CONFIG" \
        --printshellcmds \
        --software-deployment-method conda apptainer \
        --conda-frontend conda \
        --apptainer-args "--bind ${WF_DIR},${LAB_GROUP_DIR},${RAW_SEQUENCING_DATA_DIR}" \
        --executor slurm \
        --profile "${WF_DIR}profile" \
        --nolock \
        --rerun-incomplete \
        || { echo "$(date) $step FAILED"; exit 1; }
done
echo "$(date) WORKFLOW FINISHED"
