#!/bin/bash
# Download paired-end runs from SRA with the SRA toolkit (prefetch + fasterq-dump) and gzip them.
# Run on the login / vscode node (compute nodes have no internet). Runs whose fastq files exist are skipped.
#
# Usage: download_sra_fastq.sh <run list (one SRR per line)> <output dir> [parallel downloads, default 3]
# Output: <output dir>/<SRR>_1.fastq.gz, <SRR>_2.fastq.gz

set -uo pipefail
RUNS=$1
OUT=$2
PAR=${3:-3}
SRA_BIN=~/libs/sratoolkit.3.0.2-ubuntu64/bin
PIGZ_BIN=~/miniforge3_Dec2024/envs/cutadapt/bin/pigz   # not PIGZ: pigz reads that variable as options
TMP=/scratch/$USER/sra_download
mkdir -p "$OUT" "$TMP"

one() {
    srr=$1
    if [[ -s $OUT/${srr}_1.fastq.gz && -s $OUT/${srr}_2.fastq.gz ]]; then echo "$srr: present, skipped"; return 0; fi
    w=$TMP/$srr
    mkdir -p "$w"
    # resumable: prefetch skips a complete .sra file; a finished fasterq-dump (log with "reads written") is kept;
    # a mate whose .fastq is gone was compressed completely (pigz removes the input only at the end)
    rm -f "$w/$srr/$srr.sra.lock"     # one worker per run: a lock left at the start is from an interrupted earlier attempt
    $SRA_BIN/prefetch --max-size 50G -O "$w" "$srr" > "$w/prefetch.log" 2>&1 || { echo "$srr: prefetch FAILED"; return 1; }
    if ! grep -q "reads written" "$w/fasterq.log" 2>/dev/null; then
        rm -rf "$w"/*.fastq "$w"/*.fastq.gz "$w/tmp"
        $SRA_BIN/fasterq-dump --split-3 -e 3 -t "$w/tmp" -O "$w" "$w/$srr" > "$w/fasterq.log" 2>&1 || { echo "$srr: fasterq-dump FAILED"; return 1; }
    fi
    for r in 1 2; do
        if [[ -f "$w/${srr}_$r.fastq" ]]; then
            $PIGZ_BIN -f -1 -p 2 "$w/${srr}_$r.fastq" || { echo "$srr: pigz FAILED"; return 1; }   # level 1: the login node has 6 cores
        fi
        [[ -s "$w/${srr}_$r.fastq.gz" ]] || { echo "$srr: ${srr}_$r.fastq.gz missing"; return 1; }
    done
    mv "$w/${srr}_1.fastq.gz" "$w/${srr}_2.fastq.gz" "$OUT/" && rm -rf "$w"
    echo "$srr: done $(date +%H:%M)"
}
export -f one
export OUT TMP SRA_BIN PIGZ_BIN

grep -v '^\s*$' "$RUNS" | xargs -P "$PAR" -I{} bash -c 'one {}'
echo "ALL FINISHED $(date)"
