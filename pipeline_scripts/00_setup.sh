#!/bin/bash

# Rename downloaded files, download snpeff, and create custom snpeff database
# Prerequisites: wget, unzip, java, metadata file with sample ID, family, and provenance codes

RENAME_FILES=true  # set to false to skip renaming
DOWNLOAD_SNPEFF=true  # set to false to skip snpeff download
CREATE_SNPEFF_DB=true  # set to false to skip snpeff db creation

METADATA_FILE='../data/metadata/fam_prov_file.tsv'
RAW_DATA_DIR='../data/raw_data/'
RENAMED_DATA_DIR='../data/renamed_data/'

SNP_EFF_DB_NAME='Egrandis_TAG0014'
SNP_EFF_GENOME_FASTA='../data/reference/Egrandisvar_TAG0014HAP1_891_v4.1.fa'
SNP_EFF_GENOME_GFF='../data/reference/Egrandisvar_TAG0014HAP1_891_v4.1.gene.gff3'

LOG='../logs/00_setup.log'

set -euo pipefail

# ========================================================================================
# Rename downloaded files based on metadata
# ========================================================================================

if [ "$RENAME_FILES" = true ]; then
    mkdir -p "$RENAMED_DATA_DIR"
    now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
    log() { echo "[$(now)] $*" | tee -a "$LOG"; }
    log "Renaming files based on metadata..."

    # for each file in raw data dir
    for filepath in "$RAW_DATA_DIR"/*; do
        filename=$(basename "$filepath")
        # extract sample ID from filename
        sample_id=$(echo "$filename" | cut -d'.' -f1)
        # find corresponding family and provenance from metadata file
        metadata_line=$(grep -w "$sample_id" "$METADATA_FILE" || true)
        family=$(echo "$metadata_line" | awk '{print $1}')
        prov_code=$(echo "$metadata_line" | awk '{print $2}')
        if [[ -z "$prov_code" ]]; then
        log "Warning: No provenance code found for sample ID '$sample_id'. Skipping file '$filename'."
        continue
        fi
        # determine new filename
        extension="${filename#*.}"
        new_filename="${family}_${prov_code}.${extension}"
        new_filepath="$RENAMED_DATA_DIR/$new_filename"

        # rename file if new filename doesn't already exist
        if [[ -f "$new_filepath" ]]; then
            log "File '$new_filename' already exists. Skipping rename for '$filename'."
        else
            cp "$filepath" "$new_filepath"
            log "Renamed '$filename' to '$new_filename'."
        fi
    done
    log "File renaming completed."
else
    now() { date "+%Y-%m-%d %H:%M:%S %Z"; }
    log() { echo "[$(now)] $*" | tee -a "$LOG"; }
    log "Skipping file renaming."
fi

# ========================================================================================
# Download latest snpeff database
# ========================================================================================

if [ "$DOWNLOAD_SNPEFF" = true ]; then
    log "Downloading snpeff database..."
    wget https://snpeff.odsp.astrazeneca.com/versions/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip -d ../
    rm snpEff_latest_core.zip
    log "Snpeff database download completed."
else
    log "Skipping snpeff database download."
fi

# ========================================================================================
# Create custom snpeff database
# ========================================================================================

if [ "$CREATE_SNPEFF_DB" = true ]; then
    log "Creating custom snpeff database..."
    SNPEFF_DIR="../snpEff"
    mkdir -p "$SNPEFF_DIR/data/$SNP_EFF_DB_NAME"
    cp "$SNP_EFF_GENOME_FASTA" "$SNPEFF_DIR/data/$SNP_EFF_DB_NAME/sequences.fa"
    cp "$SNP_EFF_GENOME_GFF" "$SNPEFF_DIR/data/$SNP_EFF_DB_NAME/genes.gff"

    # Add entry to snpEff.config if not already present
    if ! grep -q "^$SNP_EFF_DB_NAME.genome" "$SNPEFF_DIR/snpEff.config"; then
        echo "$SNP_EFF_DB_NAME.genome : Eucalyptus grandis TAG0014" >> "$SNPEFF_DIR/snpEff.config"
    fi

    # Build the snpeff database
    cd "$SNPEFF_DIR" || exit
    java -Xmx8g -jar snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v "$SNP_EFF_DB_NAME"
    cd - || exit
    log "Custom snpeff database creation completed."
else
    log "Skipping custom snpeff database creation."
fi

log "Setup script completed."