#!/bin/bash
#SBATCH --job-name=nanoflowz_head
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=nanoflowz_head_%j.log

source ~/miniforge3_Dec2024/etc/profile.d/conda.sh
conda activate nextflow
export NXF_OPTS="-Xms2g -Xmx16g"   # keep the Java heap below the SLURM --mem

cd ~/Projects/APA_localization
nextflow run /scicore/home/zavolan/mirono0000/libs/nanoflowz/main.nf \
    -params-file /scicore/home/zavolan/GROUP/APA_localization/wf_runs/v2p0_ONT_FACSorganelles_ER_TIS_PEX_Golgi_Aleksei_v1/run_params.json \
    -profile conda -resume
