#!/bin/bash

# Filter for stop-gained variants from the SnpEff annotated VCF

INPUT_FILE='../results/04_snpeff/merged.anno.vcf.gz'
SNPEFF_DIRECTORY='../snpEff'
OUTPUT_VCF='../results/05_stop_gained/merged.stop_gained.vcf.gz'

LOG='../logs/05_stop_gained.log'

set -euo pipefail

now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
log() { echo "[$(now)] $*" | tee -a "$LOG"; }

mkdir -p "$(dirname "$OUTPUT_VCF")"

# =====================================================================
# Extract stop_gained variants
# =====================================================================
log "Extracting stop_gained variants..."
java -Xmx8g -jar $SNPEFF_DIRECTORY/SnpSift.jar filter "ANN =~ 'stop_gained'" "$INPUT_FILE" | bgzip > "$OUTPUT_VCF"
tabix -p vcf "$OUTPUT_VCF"
log "Stop-gained variants extraction completed. Output saved to $OUTPUT_VCF"

# =====================================================================
# Extract stop-gained SNPs
# =====================================================================
log "Creating stop_gained SNPs only VCF..."
bcftools view -v snps "$OUTPUT_VCF" -Oz -o "${OUTPUT_VCF%.vcf.gz}.snps.vcf.gz"
tabix -p vcf "${OUTPUT_VCF%.vcf.gz}.snps.vcf.gz"
log "SNPs only VCF created at ${OUTPUT_VCF%.vcf.gz}.snps.vcf.gz"

# =====================================================================
# Extract high impact variants
# =====================================================================
log "Extract high impact variants..."
java -Xmx8g -jar "$SNPEFF_DIRECTORY/SnpSift.jar" filter "(ANN[*].IMPACT has 'HIGH')" "$INPUT_FILE" | bgzip > "${OUTPUT_VCF%.vcf.gz}.high_impact.vcf.gz"
tabix -p vcf "${OUTPUT_VCF%.vcf.gz}.high_impact.vcf.gz"
log "High impact variants extraction completed. Output saved to ${OUTPUT_VCF%.vcf.gz}.high_impact.vcf.gz"