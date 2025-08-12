#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# === USER ADJUSTABLE ===
RAW_DIR="/scratch/alpine/beyo2625/arthur_truffles/trimmed_reads/all_reads"
TRIM_DIR="/scratch/alpine/beyo2625/arthur_truffles/trimmed_reads/all_reads_trimmed"
CONT_DIR="/scratch/alpine/beyo2625/arthur_truffles/contam_clean_reads"
OUTFILE="read_summary.tsv"

# Patterns: customize for different naming
# Example assumptions:
#   raw:    ${sample}_R1.fastq[.gz] and ${sample}_R2.fastq[.gz]
#   trimmed: ${sample}_R1.trimmed.fastq[.gz] and ${sample}_R2.trimmed.fastq[.gz]
#   contaminant-removed: ${sample}_clean_R1.fastq[.gz] and ${sample}_clean_R2.fastq[.gz]
# === end user section ===

count_reads_in_file() {
    local f=$1
    if [[ ! -e "$f" ]]; then
        echo 0
        return
    fi
    if [[ "$f" == *.gz ]]; then
        local lines
        lines=$(gzip -dc "$f" | wc -l)
    else
        local lines
        lines=$(wc -l < "$f")
    fi
    echo $(( lines / 4 ))
}

# Build sample list (deduplicated base names)
SAMPLES=$(ls "${RAW_DIR}"/*.fastq* 2>/dev/null | \
    sed -E 's#.*/##; s/_R[12].*//' | sort -u)

if [[ -z "$SAMPLES" ]]; then
    echo "No raw fastq files found in $RAW_DIR. Check paths/patterns." >&2
    exit 1
fi

printf "sample\traw_reads\ttrimmed_reads\tcontam_removed_reads_150\tcontam_removed_reads_200\n" > "$OUTFILE"

for sample in $SAMPLES; do
    # Raw
    raw_r1=$(ls "${RAW_DIR}/${sample}_R1_paired.fastq"* 2>/dev/null | head -n1 || true)
    raw_r2=$(ls "${RAW_DIR}/${sample}_R2_paired.fastq"* 2>/dev/null | head -n1 || true)
    raw1_ct=$(count_reads_in_file "$raw_r1")
    raw2_ct=$(count_reads_in_file "$raw_r2")
    raw_total=$(( raw1_ct + raw2_ct ))

    # Trimmed
    trim_r1=$(ls "${TRIM_DIR}/${sample}_val_1.fq"* 2>/dev/null | head -n1 || true)
    trim_r2=$(ls "${TRIM_DIR}/${sample}_val_2.fq"* 2>/dev/null | head -n1 || true)
    trim1_ct=$(count_reads_in_file "$trim_r1")
    trim2_ct=$(count_reads_in_file "$trim_r2")
    trim_total=$(( trim1_ct + trim2_ct ))

    # Contaminant removed bit score 100
    cont_100_r1=$(ls "${CONT_DIR}/${sample}_fwd_150.paired.fastq.gz"* 2>/dev/null | head -n1 || true)
    cont_100_r2=$(ls "${CONT_DIR}/${sample}_rev_150.paired.fastq.gz"* 2>/dev/null | head -n1 || true)
    cont1_100_ct=$(count_reads_in_file "$cont_100_r1")
    cont2_100_ct=$(count_reads_in_file "$cont_100_r2")
    cont_100_total=$(( cont1_100_ct + cont2_100_ct ))
    
    # Contaminant removed bit score 200
    cont_200_r1=$(ls "${CONT_DIR}/${sample}_fwd_200.paired.fastq.gz"* 2>/dev/null | head -n1 || true)
    cont_200_r2=$(ls "${CONT_DIR}/${sample}_rev_200.paired.fastq.gz"* 2>/dev/null | head -n1 || true)
    cont1_200_ct=$(count_reads_in_file "$cont_200_r1")
    cont2_200_ct=$(count_reads_in_file "$cont_200_r2")
    cont_200_total=$(( cont1_200_ct + cont2_200_ct ))

    printf "%s\t%d\t%d\t%d\t%d\n" "$sample" "$raw_total" "$trim_total" "$cont_100_total" "$cont_200_total" >> "$OUTFILE"
done

echo "Summary written to $OUTFILE"
