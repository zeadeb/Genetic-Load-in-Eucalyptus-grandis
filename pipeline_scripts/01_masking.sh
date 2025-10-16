#!/bin/bash

TEST_MODE=false  # set to false to produce masked output files
INPUT_DIRECTORY='../data/renamed_data/'
OUTPUT_DIRECTORY='../data/masked_data/'
MASKING_TESTING_OUTPUT='../results/01_masking/masking_test.txt'
MASKING_RESULTS='../results/01_masking/masking_results.txt'

INFO_DP_TH=10
FMT_DP_TH=8
FMT_GQ_TH=20

LOG='../logs/01_masking.log'

set -euo pipefail

now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
log() { echo "[$(now)] $*" | tee -a "$LOG"; }

# ========================================================================================
# Test the effect of masking
# ========================================================================================

shopt -s nullglob
if [ "$TEST_MODE" = true ]; then
    log "Test mode enabled..."
    echo -e "Sample\tMasking_Filter\tNum_Records\tNum_Variants\tNum_Non_Var_Masked\tNum_Variants_Masked\tPercentage_Variants_Masked" > "$MASKING_TESTING_OUTPUT"
    for filepath in "$INPUT_DIRECTORY"/*.vcf.gz; do
        base="${filepath%.bam.flt.g.vcf.gz}"
        bn="$(basename "$base")"

        log "Processing file: $bn"
        # Counts before masking
        num_records=$(bcftools view -H "$filepath" | wc -l)
        num_variants_before=$(bcftools view -H -e 'ALT="."' "$filepath" | wc -l)

        # Thresholds
        INFO_DP_TH_PLUS2=$((INFO_DP_TH + 2))
        FMT_DP_TH_MINUS2=$((FMT_DP_TH - 2))
        FMT_DP_TH_MINUS1=$((FMT_DP_TH - 1))

        # Apply INFO/DP filter
        after_info_dp_var=$(bcftools view -i "INFO/DP<$INFO_DP_TH && ALT!='.'" -H "$filepath" | wc -l)
        after_info_dp_nonvar=$(bcftools view -i "INFO/DP<$INFO_DP_TH && ALT='.'" -H "$filepath" | wc -l)

        after_info_dp_2_var=$(bcftools view -i "INFO/DP<$INFO_DP_TH_PLUS2 && ALT!='.'" -H "$filepath" | wc -l)
        after_info_dp_2_nonvar=$(bcftools view -i "INFO/DP<$INFO_DP_TH_PLUS2 + 2 && ALT='.'" -H "$filepath" | wc -l)

        # Apply FMT/DP filter
        after_fmt_dp_var=$(bcftools view -i "FMT/DP<$FMT_DP_TH && ALT!='.'" -H "$filepath" | wc -l)
        after_fmt_dp_nonvar=$(bcftools view -i "FMT/DP<$FMT_DP_TH && ALT='.'" -H "$filepath" | wc -l)

        after_fmt_dp_1_var=$(bcftools view -i "FMT/DP<$FMT_DP_TH_MINUS1 && ALT!='.'" -H "$filepath" | wc -l)
        after_fmt_dp_1_nonvar=$(bcftools view -i "FMT/DP<$FMT_DP_TH_MINUS1 && ALT='.'" -H "$filepath" | wc -l) 

        after_fmt_dp_2_var=$(bcftools view -i "FMT/DP<$FMT_DP_TH_MINUS2 && ALT!='.'" -H "$filepath" | wc -l)
        after_fmt_dp_2_nonvar=$(bcftools view -i "FMT/DP<$FMT_DP_TH_MINUS2 && ALT='.'" -H "$filepath" | wc -l)

        # Apply FMT/GQ filter
        after_fmt_gq_var=$(bcftools view -i "FMT/GQ<$FMT_GQ_TH && ALT!='.'" -H "$filepath" | wc -l)
        after_fmt_gq_nonvar=$(bcftools view -i "FMT/GQ<$FMT_GQ_TH && ALT='.'" -H "$filepath" | wc -l)

        # Apply info and fmt dp
        after_info_fmt_dp_var=$(bcftools view -i "(INFO/DP<$INFO_DP_TH || FMT/DP<$FMT_DP_TH) && ALT!='.'" -H "$filepath" | wc -l)
        after_info_fmt_dp_nonvar=$(bcftools view -i "(INFO/DP<$INFO_DP_TH || FMT/DP<$FMT_DP_TH) && ALT='.'" -H "$filepath" | wc -l)

        # Apply info and fmt gq
        after_info_fmt_gq_var=$(bcftools view -i "(INFO/DP<$INFO_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT!='.'" -H "$filepath" | wc -l)
        after_info_fmt_gq_nonvar=$(bcftools view -i "(INFO/DP<$INFO_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT='.'" -H "$filepath" | wc -l)

        # Apply fmt dp and gq
        after_fmt_fmt_gq_var=$(bcftools view -i "(FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT!='.'" -H "$filepath" | wc -l)
        after_fmt_fmt_gq_nonvar=$(bcftools view -i "(FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT='.'" -H "$filepath" | wc -l)

        # Apply all three filters
        after_all_var=$(bcftools view -i "(INFO/DP<$INFO_DP_TH || FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT!='.'" -H "$filepath" | wc -l)
        after_all_nonvar=$(bcftools view -i "(INFO/DP<$INFO_DP_TH || FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT='.'" -H "$filepath" | wc -l)

        # Add all the lines to the output file
        echo -e "$bn\tINFO/DP<$INFO_DP_TH\t$num_records\t$num_variants_before\t$after_info_dp_nonvar\t$after_info_dp_var\t$(awk -v b="$num_variants_before" -v a="$after_info_dp_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tINFO/DP<$((INFO_DP_TH + 2))\t$num_records\t$num_variants_before\t$after_info_dp_2_nonvar\t$after_info_dp_2_var\t$(awk -v b="$num_variants_before" -v a="$after_info_dp_2_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tFMT/DP<$FMT_DP_TH\t$num_records\t$num_variants_before\t$after_fmt_dp_nonvar\t$after_fmt_dp_var\t$(awk -v b="$num_variants_before" -v a="$after_fmt_dp_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tFMT/DP<$FMT_DP_TH_MINUS1\t$num_records\t$num_variants_before\t$after_fmt_dp_1_nonvar\t$after_fmt_dp_1_var\t$(awk -v b="$num_variants_before" -v a="$after_fmt_dp_1_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT" 

        echo -e "$bn\tFMT/DP<$((FMT_DP_TH - 2))\t$num_records\t$num_variants_before\t$after_fmt_dp_2_nonvar\t$after_fmt_dp_2_var\t$(awk -v b="$num_variants_before" -v a="$after_fmt_dp_2_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tFMT/GQ<$FMT_GQ_TH\t$num_records\t$num_variants_before\t$after_fmt_gq_nonvar\t$after_fmt_gq_var\t$(awk -v b="$num_variants_before" -v a="$after_fmt_gq_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tINFO/DP<$INFO_DP_TH || FMT/DP<$FMT_DP_TH\t$num_records\t$num_variants_before\t$after_info_fmt_dp_nonvar\t$after_info_fmt_dp_var\t$(awk -v b="$num_variants_before" -v a="$after_info_fmt_dp_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tINFO/DP<$INFO_DP_TH || FMT/GQ<$FMT_GQ_TH\t$num_records\t$num_variants_before\t$after_info_fmt_gq_nonvar\t$after_info_fmt_gq_var\t$(awk -v b="$num_variants_before" -v a="$after_info_fmt_gq_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tFMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH\t$num_records\t$num_variants_before\t$after_fmt_fmt_gq_nonvar\t$after_fmt_fmt_gq_var\t$(awk -v b="$num_variants_before" -v a="$after_fmt_fmt_gq_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

        echo -e "$bn\tINFO/DP<$INFO_DP_TH || FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH\t$num_records\t$num_variants_before\t$after_all_nonvar\t$after_all_var\t$(awk -v b="$num_variants_before" -v a="$after_all_var" 'BEGIN {printf "%.2f", (a/b)*100}')" >> "$MASKING_TESTING_OUTPUT"

    done
    log "Masking test completed. Results saved to '$MASKING_TESTING_OUTPUT'."
else
    log "Masking files..."
    mkdir -p "$OUTPUT_DIRECTORY"
    echo -e "Sample\tNum_Records\tNum_Variants\tNum_Non_Var_Masked\tNum_Variants_Masked" > "$MASKING_RESULTS"
    for filepath in "$INPUT_DIRECTORY"/*.vcf.gz; do
        base="${filepath%.bam.flt.g.vcf.gz}"
        bn="$(basename "$base")"
        log "Processing file: $bn"

        # Counts before masking
        num_records=$(bcftools view -H "$filepath" | wc -l)
        num_variants_before=$(bcftools view -H -e 'ALT="."' "$filepath" | wc -l)

        # Counts after masking
        num_variants_masked=$(bcftools view -i "(FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT!='.'" -H "$filepath" | wc -l)
        num_nonvar_masked=$(bcftools view -i "(FMT/DP<$FMT_DP_TH || FMT/GQ<$FMT_GQ_TH) && ALT='.'" -H "$filepath" | wc -l)

        # Log the results
        echo -e "$bn\t$num_records\t$num_variants_before\t$num_nonvar_masked\t$num_variants_masked" >> "$MASKING_RESULTS"

        # Apply masking and save to output directory
        bcftools +setGT "$filepath" -Oz -o "$OUTPUT_DIRECTORY/${bn}.gtmask.vcf.gz" -- -t q -n . -i "FMT/GQ<$FMT_GQ_TH || FMT/DP<$FMT_DP_TH"

    done
    log "Masking completed. Masked files saved to '$OUTPUT_DIRECTORY'."
fi