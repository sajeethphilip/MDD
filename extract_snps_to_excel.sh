#!/bin/bash

############################################################
# MDD2 SNP Data Extractor for Excel Output
# Extracts SNP information from TSV files and creates Excel report
# Robust version - handles missing data gracefully
############################################################

set -euo pipefail

# Configuration - try to auto-detect paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/mdd2.sh" ]]; then
    # If running from same directory as mdd2.sh
    BASE_DIR="${SCRIPT_DIR}/analysis2"
else
    # Try common locations
    if [[ -d "./analysis2" ]]; then
        BASE_DIR="./analysis2"
    elif [[ -d "../analysis2" ]]; then
        BASE_DIR="../analysis2"
    else
        BASE_DIR="${SCRIPT_DIR}"
        echo "WARNING: Could not auto-detect analysis2 directory. Using: ${BASE_DIR}"
    fi
fi

TSV_DIR="${BASE_DIR}/data/tsv"
OUTPUT_DIR="${BASE_DIR}/reports"
OUTPUT_EXCEL="${OUTPUT_DIR}/mdd_snp_summary.xlsx"
OUTPUT_TSV="${OUTPUT_DIR}/mdd_snp_summary.tsv"

mkdir -p "${OUTPUT_DIR}"

# Logging function
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] INFO: $1"
}

log_warning() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] WARNING: $1" >&2
}

# Function to extract and format SNP data (matrix format)
extract_snp_data_matrix() {
    log "Starting SNP data extraction (matrix format)..."

    # Create temporary working directory
    TEMP_DIR="${OUTPUT_DIR}/temp"
    mkdir -p "${TEMP_DIR}"

    # Step 1: Find all sample directories
    local sample_dirs=()
    if [[ -d "${TSV_DIR}" ]]; then
        sample_dirs=("${TSV_DIR}"/*/)
    fi

    if [[ ${#sample_dirs[@]} -eq 0 ]] || [[ ! -d "${TSV_DIR}" ]]; then
        log_warning "No sample directories found in ${TSV_DIR}"
        log "Looking for TSV files in other locations..."

        # Try to find TSV files anywhere
        local tsv_files=($(find "${BASE_DIR}" -name "*.tsv.gz" -type f 2>/dev/null | head -20))
        if [[ ${#tsv_files[@]} -eq 0 ]]; then
            tsv_files=($(find "${BASE_DIR}" -name "*.tsv" -type f 2>/dev/null | head -20))
        fi

        if [[ ${#tsv_files[@]} -eq 0 ]]; then
            log_warning "No TSV files found anywhere in ${BASE_DIR}"
            echo "Please ensure the pipeline has generated TSV files."
            echo "You can run: ./mdd2.sh run"
            rm -rf "${TEMP_DIR}"
            return 1
        fi

        log "Found ${#tsv_files[@]} TSV files in other locations"

        # Create sample directories structure from found files
        for tsv_file in "${tsv_files[@]}"; do
            local sample_id=$(basename "${tsv_file}" .tsv.gz)
            sample_id=$(basename "${sample_id}" .tsv)
            local sample_temp_dir="${TEMP_DIR}/${sample_id}"
            mkdir -p "${sample_temp_dir}"

            if [[ "${tsv_file}" == *.gz ]]; then
                cp "${tsv_file}" "${sample_temp_dir}/${sample_id}.tsv.gz"
            else
                gzip -c "${tsv_file}" > "${sample_temp_dir}/${sample_id}.tsv.gz" 2>/dev/null || true
            fi
        done

        sample_dirs=("${TEMP_DIR}"/*/)
    fi

    log "Found ${#sample_dirs[@]} samples to process"

    # Step 2: Extract sample IDs
    local sample_ids=()
    for sample_dir in "${sample_dirs[@]}"; do
        local sample_id=$(basename "${sample_dir}")
        sample_ids+=("${sample_id}")
    done

    # Step 3: Initialize data structures
    declare -A snp_data
    declare -A snp_sample_map

    # Step 4: Process each sample
    local processed_samples=0
    for sample_id in "${sample_ids[@]}"; do
        local tsv_file=""

        # Try multiple locations
        if [[ -f "${TSV_DIR}/${sample_id}/${sample_id}.tsv.gz" ]]; then
            tsv_file="${TSV_DIR}/${sample_id}/${sample_id}.tsv.gz"
        elif [[ -f "${TSV_DIR}/${sample_id}.tsv.gz" ]]; then
            tsv_file="${TSV_DIR}/${sample_id}.tsv.gz"
        elif [[ -f "${TEMP_DIR}/${sample_id}/${sample_id}.tsv.gz" ]]; then
            tsv_file="${TEMP_DIR}/${sample_id}/${sample_id}.tsv.gz"
        elif [[ -f "${TEMP_DIR}/${sample_id}.tsv.gz" ]]; then
            tsv_file="${TEMP_DIR}/${sample_id}.tsv.gz"
        fi

        if [[ -z "${tsv_file}" ]] || [[ ! -f "${tsv_file}" ]]; then
            log_warning "TSV file not found for sample ${sample_id}, skipping"
            continue
        fi

        log "Processing sample: ${sample_id}"

        # Decompress and process TSV
        local decompressed_file="${TEMP_DIR}/${sample_id}_decompressed.tsv"
        if command -v pigz > /dev/null 2>&1; then
            pigz -dc "${tsv_file}" 2>/dev/null > "${decompressed_file}" || {
                log_warning "Failed to decompress ${tsv_file} with pigz, trying gzip..."
                gzip -dc "${tsv_file}" 2>/dev/null > "${decompressed_file}" || {
                    log_warning "Failed to decompress ${tsv_file}, skipping sample"
                    continue
                }
            }
        else
            gzip -dc "${tsv_file}" 2>/dev/null > "${decompressed_file}" || {
                log_warning "Failed to decompress ${tsv_file}, skipping sample"
                continue
            }
        fi

        # Parse TSV file
        local line_num=0
        while IFS=$'\t' read -r -a fields; do
            ((line_num++))

            # Skip header and empty lines
            if [[ $line_num -eq 1 ]] || [[ ${#fields[@]} -lt 8 ]]; then
                continue
            fi

            # Extract fields with safe defaults
            local chrom="${fields[0]:-}"
            local pos="${fields[1]:-}"
            local id="${fields[2]:-.}"
            local ref="${fields[3]:-}"
            local alt="${fields[4]:-}"
            local gene="${fields[7]:-UNKNOWN}"

            # Skip if missing essential fields
            if [[ -z "${chrom}" ]] || [[ -z "${pos}" ]] || [[ -z "${ref}" ]] || [[ -z "${alt}" ]]; then
                continue
            fi

            # Create unique SNP key
            local snp_key="${chrom}:${pos}:${ref}:${alt}:${gene}"

            # Store SNP basic data if not already stored
            if [[ -z "${snp_data[${snp_key}]+x}" ]]; then
                # Format chromosome
                local formatted_chrom="${chrom}"
                # Simple chromosome formatting - you can expand this
                if [[ "${chrom}" == chr* ]]; then
                    local chr_num=${chrom#chr}
                    case "${chr_num}" in
                        1) formatted_chrom="NC_000001.11" ;;
                        2) formatted_chrom="NC_000002.12" ;;
                        3) formatted_chrom="NC_000003.12" ;;
                        4) formatted_chrom="NC_000004.12" ;;
                        5) formatted_chrom="NC_000005.10" ;;
                        6) formatted_chrom="NC_000006.12" ;;
                        7) formatted_chrom="NC_000007.14" ;;
                        8) formatted_chrom="NC_000008.11" ;;
                        9) formatted_chrom="NC_000009.12" ;;
                        10) formatted_chrom="NC_000010.11" ;;
                        11) formatted_chrom="NC_000011.10" ;;
                        12) formatted_chrom="NC_000012.12" ;;
                        13) formatted_chrom="NC_000013.11" ;;
                        14) formatted_chrom="NC_000014.9" ;;
                        15) formatted_chrom="NC_000015.10" ;;
                        16) formatted_chrom="NC_000016.10" ;;
                        17) formatted_chrom="NC_000017.11" ;;
                        18) formatted_chrom="NC_000018.10" ;;
                        19) formatted_chrom="NC_000019.10" ;;
                        20) formatted_chrom="NC_000020.11" ;;
                        21) formatted_chrom="NC_000021.9" ;;
                        22) formatted_chrom="NC_000022.11" ;;
                        X) formatted_chrom="NC_000023.11" ;;
                        Y) formatted_chrom="NC_000024.10" ;;
                        M|MT) formatted_chrom="NC_012920.1" ;;
                        *) formatted_chrom="${chrom}" ;;
                    esac
                fi

                # Store SNP data
                snp_data["${snp_key}"]="${formatted_chrom}\t${pos}\t${id}\t${ref}\t${alt}\t${gene}"
            fi

            # Mark this sample as having this SNP (safely)
            snp_sample_map["${snp_key}:${sample_id}"]="1"

        done < "${decompressed_file}"

        # Clean up
        rm -f "${decompressed_file}"
        ((processed_samples++))
    done

    if [[ ${processed_samples} -eq 0 ]]; then
        log_warning "No samples were successfully processed"
        rm -rf "${TEMP_DIR}"
        return 1
    fi

    log "Successfully processed ${processed_samples} samples"
    log "Found ${#snp_data[@]} unique SNPs"

    # Step 5: Generate output matrix
    log "Step 5: Generating SNP matrix..."

    # Write header
    {
        echo -n -e "SAMPLE ID\tGENE\tNo. of SNPs\tCHR\tPOS\tRSID\tREF\tALT"
        for sample_id in "${sample_ids[@]}"; do
            echo -n -e "\t${sample_id}"
        done
        echo ""
    } > "${OUTPUT_TSV}"

    # Count SNPs per gene per sample
    declare -A gene_snp_count

    # First, count SNPs for each gene in each sample
    for snp_key in "${!snp_sample_map[@]}"; do
        # Extract gene from snp_key (format: chrom:pos:ref:alt:gene:sample)
        local gene=$(echo "${snp_key}" | cut -d: -f5)
        local sample_id=$(echo "${snp_key}" | cut -d: -f6)

        if [[ -n "${gene}" ]] && [[ -n "${sample_id}" ]]; then
            local key="${sample_id}:${gene}"
            gene_snp_count["${key}"]=$(( ${gene_snp_count["${key}"]:-0} + 1 ))
        fi
    done

    # Process each unique SNP
    for snp_key in "${!snp_data[@]}"; do
        IFS=$'\t' read -r chrom pos id ref alt gene <<< "${snp_data[${snp_key}]}"

        # For each sample, check if it has this SNP
        local sample_presence=""
        for sample_id in "${sample_ids[@]}"; do
            if [[ -n "${snp_sample_map[${snp_key}:${sample_id}]+x}" ]]; then
                sample_presence+="\t1"
            else
                sample_presence+="\t"
            fi
        done

        # Get SNP count for this gene in the first sample that has it
        local snp_count="0"
        for sample_id in "${sample_ids[@]}"; do
            local key="${sample_id}:${gene}"
            if [[ -n "${gene_snp_count[${key}]+x}" ]]; then
                snp_count="${gene_snp_count[${key}]}"
                # Use this sample ID for the first column
                local display_sample="${sample_id}"
                # Write row
                echo -e "${display_sample}\t${gene}\t${snp_count}\t${chrom}\t${pos}\t${id}\t${ref}\t${alt}${sample_presence}" >> "${OUTPUT_TSV}"
                break
            fi
        done
    done

    # Sort the output file
    sort -t$'\t' -k4,4 -k5,5n "${OUTPUT_TSV}" > "${OUTPUT_TSV}.sorted"
    mv "${OUTPUT_TSV}.sorted" "${OUTPUT_TSV}"

    log "SNP matrix created: ${OUTPUT_TSV}"

    # Step 6: Convert to Excel if possible
    convert_to_excel "${OUTPUT_TSV}" "${OUTPUT_EXCEL}"

    # Cleanup
    rm -rf "${TEMP_DIR}"

    log "Extraction complete!"
    return 0
}

# Function to extract data in detailed format
extract_snp_data_detailed() {
    log "Starting detailed SNP data extraction..."

    local output_file="${OUTPUT_DIR}/snp_data_detailed.tsv"

    # Write header
    echo -e "SAMPLE_ID\tGENE\tSNP_COUNT\tCHROM\tPOS\tRSID\tREF\tALT\tGT\tDP\tAD\tGQ\tQUAL\tFILTER" > "${output_file}"

    # Find and process TSV files
    local processed_files=0

    # Look for TSV files in various locations
    find "${BASE_DIR}" -name "*.tsv.gz" -type f 2>/dev/null | while read -r tsv_file; do
        local sample_id=$(basename "${tsv_file}" .tsv.gz)
        log "Processing: ${sample_id}"

        # Try to decompress and parse
        local temp_file="${OUTPUT_DIR}/temp_${sample_id}.tsv"

        if command -v pigz > /dev/null 2>&1; then
            pigz -dc "${tsv_file}" 2>/dev/null > "${temp_file}" || continue
        else
            gzip -dc "${tsv_file}" 2>/dev/null > "${temp_file}" || continue
        fi

        # Count SNPs per gene
        declare -A gene_counts
        while IFS=$'\t' read -r -a fields; do
            if [[ ${#fields[@]} -ge 8 ]] && [[ "${fields[0]}" != "CHROM" ]]; then
                local gene="${fields[7]:-UNKNOWN}"
                gene_counts["${gene}"]=$(( ${gene_counts["${gene}"]:-0} + 1 ))
            fi
        done < "${temp_file}"

        # Extract and format data
        while IFS=$'\t' read -r -a fields; do
            # Skip header and incomplete lines
            if [[ ${#fields[@]} -lt 8 ]] || [[ "${fields[0]}" == "CHROM" ]]; then
                continue
            fi

            local chrom="${fields[0]}"
            local pos="${fields[1]}"
            local id="${fields[2]:-.}"
            local ref="${fields[3]}"
            local alt="${fields[4]}"
            local qual="${fields[5]:-}"
            local filter="${fields[6]:-}"
            local gene="${fields[7]:-UNKNOWN}"
            local gt="${fields[10]:-}"    # Adjust indices based on your TSV format
            local dp="${fields[11]:-}"
            local ad="${fields[12]:-}"
            local gq="${fields[13]:-}"

            local snp_count="${gene_counts["${gene}"]:-0}"

            # Format chromosome
            local formatted_chrom="${chrom}"
            if [[ "${chrom}" == "chr4" ]]; then
                formatted_chrom="NC_000004.12"
            fi

            echo -e "${sample_id}\t${gene}\t${snp_count}\t${formatted_chrom}\t${pos}\t${id}\t${ref}\t${alt}\t${gt}\t${dp}\t${ad}\t${gq}\t${qual}\t${filter}"
        done < "${temp_file}" >> "${output_file}"

        rm -f "${temp_file}"
        ((processed_files++))
    done

    if [[ ${processed_files} -eq 0 ]]; then
        # Try uncompressed TSV files
        find "${BASE_DIR}" -name "*.tsv" -type f 2>/dev/null | while read -r tsv_file; do
            local sample_id=$(basename "${tsv_file}" .tsv)
            log "Processing uncompressed: ${sample_id}"
            ((processed_files++))
            # Similar processing for uncompressed files...
        done
    fi

    if [[ ${processed_files} -eq 0 ]]; then
        log_warning "No TSV files found for detailed extraction"
        return 1
    fi

    log "Detailed SNP data extracted to: ${output_file}"
    convert_to_excel "${output_file}" "${OUTPUT_DIR}/snp_data_detailed.xlsx"

    return 0
}

# Function to convert TSV to Excel
convert_to_excel() {
    local input_tsv="$1"
    local output_excel="$2"

    if [[ ! -f "${input_tsv}" ]] || [[ ! -s "${input_tsv}" ]]; then
        log_warning "Input TSV file not found or empty: ${input_tsv}"
        return 1
    fi

    # Check if Python is available
    if ! command -v python3 > /dev/null 2>&1; then
        log_warning "Python3 not available. Cannot create Excel file."
        log "TSV file is available at: ${input_tsv}"
        return 1
    fi

    # Check for pandas
    if ! python3 -c "import pandas" 2>/dev/null; then
        log_warning "Python pandas not installed. Installing..."

        # Try to install pandas
        if command -v pip3 > /dev/null 2>&1; then
            pip3 install --user pandas openpyxl 2>/dev/null || {
                log_warning "Failed to install pandas. Creating TSV only."
                return 1
            }
        elif command -v pip > /dev/null 2>&1; then
            pip install --user pandas openpyxl 2>/dev/null || {
                log_warning "Failed to install pandas. Creating TSV only."
                return 1
            }
        else
            log_warning "pip not available. Cannot install pandas."
            return 1
        fi
    fi

    log "Converting to Excel format: ${output_excel}"

    python3 -c "
import pandas as pd
import sys
import os

try:
    input_file = '${input_tsv}'
    output_file = '${output_excel}'

    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    # Replace NaN with empty string
    df = df.fillna('')

    # Create Excel writer
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Main data sheet
        df.to_excel(writer, sheet_name='SNP_Data', index=False)

        # Create summary statistics if we have sample columns
        sample_cols = [col for col in df.columns if col.startswith('SRR')]
        if sample_cols:
            summary_data = []
            for sample in sample_cols:
                # Count non-empty cells for this sample
                sample_snps = df[sample].astype(str).str.strip().ne('').sum()
                # Count unique genes with SNPs
                genes_with_snps = df[df[sample].astype(str).str.strip().ne('')]['GENE'].nunique()

                summary_data.append({
                    'Sample': sample,
                    'Total SNPs': sample_snps,
                    'Genes with SNPs': genes_with_snps,
                    'SNPs per Gene': round(sample_snps / max(genes_with_snps, 1), 2)
                })

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name='Summary', index=False)

        # Auto-adjust column widths
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]

            # Find the maximum length in each column
            for col in worksheet.columns:
                max_length = 0
                column = col[0].column_letter  # Get the column name

                for cell in col:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass

                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column].width = adjusted_width

        print(f'Successfully created Excel file: {output_file}')
        print(f'Sheet 1: SNP_Data - Main SNP matrix')
        if sample_cols:
            print(f'Sheet 2: Summary - Sample statistics')

except Exception as e:
    print(f'Error creating Excel file: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)
" 2>&1 | tee "${OUTPUT_DIR}/excel_conversion.log"

    if [[ $? -eq 0 ]] && [[ -f "${output_excel}" ]]; then
        log "Excel file created successfully: ${output_excel}"
        return 0
    else
        log_warning "Excel conversion had issues. TSV file is available at: ${input_tsv}"
        return 1
    fi
}

# Main function
main() {
    echo ""
    echo "========================================"
    echo "MDD2 SNP Data Extractor (Robust Version)"
    echo "========================================"
    echo ""

    # Display detected paths
    echo "Configuration:"
    echo "  Base directory: ${BASE_DIR}"
    echo "  TSV directory: ${TSV_DIR}"
    echo "  Output directory: ${OUTPUT_DIR}"
    echo ""

    # Check if any TSV files exist
    local tsv_count=$(find "${BASE_DIR}" -name "*.tsv.gz" -type f 2>/dev/null | wc -l)
    local tsv_count2=$(find "${BASE_DIR}" -name "*.tsv" -type f 2>/dev/null | wc -l)
    local total_tsv=$((tsv_count + tsv_count2))

    if [[ ${total_tsv} -eq 0 ]]; then
        echo "WARNING: No TSV files found in ${BASE_DIR}"
        echo "You can:"
        echo "  1. Run the pipeline first: ./mdd2.sh run"
        echo "  2. Place existing TSV files in: ${TSV_DIR}/"
        echo "  3. Continue anyway to create template files"
        echo ""
        read -p "Continue anyway? [y/N]: " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Exiting. Please run the pipeline first."
            exit 0
        fi
    else
        echo "Found ${total_tsv} TSV files"
    fi

    # Ask user which extraction method to use
    echo ""
    echo "Extraction Options:"
    echo "  1. Matrix format (samples as columns - recommended)"
    echo "  2. Detailed format (all data in one table)"
    echo "  3. Both formats"
    echo ""
    read -p "Select option [1/2/3]: " -n 1 -r
    echo

    local success_count=0

    case $REPLY in
        1)
            if extract_snp_data_matrix; then
                ((success_count++))
            fi
            ;;
        2)
            if extract_snp_data_detailed; then
                ((success_count++))
            fi
            ;;
        3)
            if extract_snp_data_matrix; then
                ((success_count++))
            fi
            if extract_snp_data_detailed; then
                ((success_count++))
            fi
            ;;
        *)
            echo "Invalid choice, using matrix format"
            if extract_snp_data_matrix; then
                ((success_count++))
            fi
            ;;
    esac

    # Display summary
    echo ""
    echo "========================================"
    echo "Extraction Summary"
    echo "========================================"
    echo ""

    if [[ ${success_count} -gt 0 ]]; then
        echo "✅ Successfully created ${success_count} output file(s)"
        echo ""
        echo "Output files in ${OUTPUT_DIR}:"
        ls -la "${OUTPUT_DIR}"/*.tsv "${OUTPUT_DIR}"/*.xlsx 2>/dev/null | while read -r line; do
            echo "  ${line}"
        done
        echo ""
        echo "To view the Excel file:"
        echo "  libreoffice ${OUTPUT_EXCEL}  # or open with Excel"
        echo ""
    else
        echo "❌ No output files were created"
        echo ""
        echo "Troubleshooting tips:"
        echo "  1. Ensure TSV files exist in: ${TSV_DIR}/"
        echo "  2. Check file permissions"
        echo "  3. Install Python pandas: pip install pandas openpyxl"
        echo "  4. Run the pipeline first: ./mdd2.sh run"
        echo ""
    fi

    # Create a simple README file
    cat > "${OUTPUT_DIR}/README.txt" << EOF
MDD2 SNP Data Extraction Output
================================

Generated on: $(date)

Output Files:
1. mdd_snp_summary.tsv - Tab-separated matrix format
2. mdd_snp_summary.xlsx - Excel version with multiple sheets
3. snp_data_detailed.tsv - Detailed per-SNP data (if generated)
4. snp_data_detailed.xlsx - Excel version of detailed data (if generated)

File Formats:
- Matrix format: Samples as columns, SNPs as rows
- Detailed format: One row per SNP per sample

To update:
1. Run the pipeline: ./mdd2.sh run
2. Re-run this extractor: ./extract_snps_to_excel.sh

EOF

    echo "README file created: ${OUTPUT_DIR}/README.txt"
    echo ""
}

# Run main function with error trapping
main "$@"
