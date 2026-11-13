#!/bin/bash
#SBATCH --job-name=spike_ins
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=/scicore/home/zavolan/GROUP/APA_localization/wf_runs/U2OS_subcellular_fractions_spike_ins/logs/%x_%A_%a.log

# SIRV read pairs per library of the U2OS density-gradient libraries (LoRNA, Villanueva et al. 2024; Lexogen SIRV
# isoform mix; ERCC92 is in the reference as a check): STAR on the spike-in sequences only, strict filters (>= 80% of the
# pair aligned, <= 4% mismatches). One array task per row of the start-samples file (header excluded).
# Usage: sbatch --array=1-<n> WF/scripts/count_spike_ins.sbatch <start_samples.tsv> <outdir>
# Output per library: <outdir>/<sample>/{Log.final.out,ReadsPerGene.out.tab}; skipped if it exists.

SAMPLES=$1; OUT=$2
INDEX=/scicore/home/zavolan/GROUP/APA_localization/external_data/spike_ins/STAR_index
read -r sample fq1 fq2 < <(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR == n + 1 {print $1, $4, $5}' "$SAMPLES")
[[ -s "$OUT/$sample/Log.final.out" ]] && { echo "$sample: done before"; exit 0; }
source ~/miniforge3_Dec2024/etc/profile.d/conda.sh; conda activate STAR
rm -rf "$OUT/$sample.tmp"; mkdir -p "$OUT/$sample.tmp"
STAR --runThreadN "$SLURM_CPUS_PER_TASK" --genomeDir "$INDEX" --readFilesIn "$fq1" "$fq2" --readFilesCommand zcat \
     --outFileNamePrefix "$OUT/$sample.tmp/" --outSAMtype None --quantMode GeneCounts \
     --outFilterMismatchNoverReadLmax 0.04 --outFilterScoreMinOverLread 0.8 --outFilterMatchNminOverLread 0.8 \
     --alignIntronMax 10000 --alignMatesGapMax 10000 --outTmpDir "/scratch/$USER/star_spike_${SLURM_JOB_ID}" || exit 1
mv "$OUT/$sample.tmp" "$OUT/$sample"
echo "$sample: $(grep 'Uniquely mapped reads number' "$OUT/$sample/Log.final.out")"
