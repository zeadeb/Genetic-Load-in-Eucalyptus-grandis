#!/bin/bash

# Run SnpEff annotation on the final filtered VCF
# Prerequisites: SnpEff database for Egrandis_TAG0014 obtained by 00_setup.sh

INPUT_FILE='../results/03_filter/merged.filtered.norm.vcf.gz'
REFERENCE_SEQ='../data/reference/Egrandisvar_TAG0014HAP1_891_v4.1.fa'
SNPEFF_DIRECTORY='../snpEff'
OUTPUT_DIRECTORY='../results/04_snpeff'

LOG='../logs/04_snpeff.log'

set -euo pipefail

mkdir -p "$OUTPUT_DIRECTORY"

now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
log() { echo "[$(now)] $*" | tee -a "$LOG"; }

# Run SnpEff
log "Running SnpEff..."

java -Xmx8g -jar "$SNPEFF_DIRECTORY/snpEff.jar" -v -stats snpeff_stats.html -csvStats snpeff_stats.csv Egrandis_TAG0014 "$INPUT_FILE" | bgzip > "$OUTPUT_DIRECTORY/merged.anno.vcf.gz"
tabix -p vcf "$OUTPUT_DIRECTORY/merged.anno.vcf.gz"

log "SnpEff annotation completed. Output saved to $OUTPUT_DIRECTORY/merged.anno.vcf.gz"

