#!/bin/bash

# Merge all gvcf files into one merged gvcf file
# Prerequisites: bcftools, bgzip, tabix

set -euo pipefail

MERGE_THREADS=30

REFERENCE_SEQ='../data/reference/Egrandisvar_TAG0014HAP1_891_v4.1.fa'
FILES_TO_MERGE='../../../../Downloads/ALL_PROV/masked_data/*.gtmask.vcf.gz'
TEMP_NORM_FILE_LIST='temp_norm_file_list.txt'
MERGED_OUTPUT='../results/02_merge/merged.g.vcf.gz'

LOG='../logs/02_merge.log'

now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
log() { echo "[$(now)] $*" | tee -a "$LOG"; }

log "Starting merge process..."
rm -f "$TEMP_NORM_FILE_LIST"

for f in $FILES_TO_MERGE; do
    # Normalize each file and index it
    log "Normalizing: $f"
    out="${f%.vcf.gz}.norm.g.vcf.gz"
    bcftools norm -f "$REFERENCE_SEQ" -m-any -Oz -o "$out" "$f"
    bcftools index -t "$out"
    # save the normalized file name for merging
    echo "$out" >> "$TEMP_NORM_FILE_LIST"
done

# Merge all normalized files
log "Merging normalized files..."
bcftools merge --threads "$MERGE_THREADS" --gvcf "$REFERENCE_SEQ" -m 'both,**' --file-list "$TEMP_NORM_FILE_LIST" -Oz -o "$MERGED_OUTPUT"
if [ $? -ne 0 ]; then
    log "Error: Merge failed. Please check the individual normalized files for issues."
    exit 1
fi

# Index the final merged file
log "Indexing merged file..."
bcftools index -t "$MERGED_OUTPUT"

# Clean up the intermediate file list
rm "$TEMP_NORM_FILE_LIST"