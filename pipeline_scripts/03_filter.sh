#!/bin/bash

# Convert masked and merged gvcf to final filtered vcf
# Prerequisites: bcftools, bgzip, tabix

set -euo pipefail

CONVERT_TO_VCF=true
SPLIT_MULTIALLELIC=true
REMOVE_NO_ALT=true
REMOVE_REF_SNPS=true
FILTER_QUAL=true

INPUT_GVCF='../results/02_merge/merged.g.vcf.gz'
OUTPUT_VCF='../results/03_filter/merged.filtered.vcf.gz'
OUTPUT_DIR='../results/03_filter'
REFERENCE_SEQ='../data/reference/Egrandisvar_TAG0014HAP1_891_v4.1.fa'
LOG='../logs/03_filter.log'

now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
log() { echo "[$(now)] $*" | tee -a "$LOG"; }
log "Starting filtering process..."

# ==========================================================
# Convert gvcf to vcf by removing non-variant sites
# ==========================================================
if [ "$CONVERT_TO_VCF" = false ]; then
    log "Skipping gvcf to vcf conversion as per configuration."
else
    log "Converting gvcf to vcf..."
    bcftools view -i 'TYPE!="ref"' -Oz -o "$OUTPUT_DIR/merged.vcf.gz" "$INPUT_GVCF"
fi

# ===========================================================
# Split multiallelic sites into multiple biallelic records and fill tags
# ===========================================================
if [ "$SPLIT_MULTIALLELIC" = false ]; then
    log "Skipping multiallelic splitting as per configuration."
else
    log "Splitting multiallelic sites and filling tags..."
    bcftools norm -f "$REFERENCE_SEQ" "$OUTPUT_DIR/merged.vcf.gz" -Ou | bcftools norm -m-any -Ou | bcftools annotate \
    -x INFO/AC,INFO/AN,INFO/AF,INFO/NS -Ou | bcftools +fill-tags -- -t AC,AN,AF,F_MISSING | \
    bcftools view -Oz -o "$OUTPUT_DIR/merged.norm.vcf.gz"
fi

# ===========================================================
# Remove sites with no ALT genotypes (all 0/0 or ./.)
# ===========================================================
if [ "$REMOVE_NO_ALT" = false ]; then
    log "Skipping removal of sites with no ALT genotypes as per configuration."
else
    log "Removing sites with no ALT genotypes..."
    bcftools view -i 'AC>0' "$OUTPUT_DIR/merged.norm.vcf.gz" -Oz -o "$OUTPUT_DIR/merged.norm.var.only.vcf.gz"
fi

# ===========================================================
# Remove REF SNPs
# ===========================================================
if [ "$REMOVE_REF_SNPS" = true ]; then
    log "Removing REF SNPs..."
    bcftools view -e 'AF=1 && AN>=318' -Oz -o "$OUTPUT_DIR/merged.norm.var.no.ref.vcf.gz" "$OUTPUT_DIR/merged.norm.var.only.vcf.gz"
fi

# ===========================================================
# Filter variants based on quality metrics
# ===========================================================
if [ "$FILTER_QUAL" = false ]; then
    log "Skipping variant filtering based on quality metrics as per configuration."
else
    log "Filtering variants based on quality metrics..."
    bcftools view -i 'QUAL>=30' "$OUTPUT_DIR/merged.norm.var.no.ref.vcf.gz" -Oz -o "$OUTPUT_VCF"
    bcftools index -t "$OUTPUT_VCF"
fi

log "Filtering process completed. Output saved to $OUTPUT_VCF"