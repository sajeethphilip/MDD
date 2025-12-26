#!/bin/bash
# mdd2_streaming_analyzer.sh - For very large datasets

# Configuration
BASE_DIR="analysis2"
THRESHOLD=95
MIN_SAMPLES=3

# Step 1: Create SNP presence/absence matrix on the fly
create_snp_matrix() {
    # Use awk streaming to avoid loading everything in memory
    awk '
    BEGIN {
        # We'll use disk-based counting for huge datasets
        control_file = "'"${BASE_DIR}"'/temp_control_counts"
        mdd_file = "'"${BASE_DIR}"'/temp_mdd_counts"
        system("mkdir -p '"${BASE_DIR}"'/temp_counts")
    }

    # Process each sample file as it arrives
    # This would be called from a find | xargs pipeline
    {
        # Process one SNP from one sample
        sample = $1
        group = $2
        snp_key = $3

        if (group ~ /[Cc]ontrol/) {
            print snp_key >> "'"${BASE_DIR}"'/temp_counts/control/" sample
        } else {
            print snp_key >> "'"${BASE_DIR}"'/temp_counts/mdd/" sample
        }
    }
    '
}

# Step 2: Merge and analyze (using sort | uniq for counting)
analyze_counts() {
    # Count occurrences using sort/uniq (memory efficient)
    find "${BASE_DIR}/temp_counts/control" -type f -exec cat {} \; | \
        sort | uniq -c | \
        awk -v threshold="${THRESHOLD}" -v min_samples="${MIN_SAMPLES}" \
        '{
            count = $1
            snp = $2
            if (count >= min_samples) {
                control_snps[snp] = count
            }
        }
        END {
            for (snp in control_snps) {
                print snp, control_snps[snp]
            }
        }' | sort -k2,2nr > "${BASE_DIR}/control_frequent_snps.tsv"

    # Same for MDD
    find "${BASE_DIR}/temp_counts/mdd" -type f -exec cat {} \; | \
        sort | uniq -c | \
        awk -v threshold="${THRESHOLD}" -v min_samples="${MIN_SAMPLES}" \
        '{
            count = $1
            snp = $2
            if (count >= min_samples) {
                mdd_snps[snp] = count
            }
        }
        END {
            for (snp in mdd_snps) {
                print snp, mdd_snps[snp]
            }
        }' | sort -k2,2nr > "${BASE_DIR}/mdd_frequent_snps.tsv"
}

# Step 3: Find differences using comm (extremely efficient)
find_differential_snps() {
    # Extract just SNP keys
    cut -f1 "${BASE_DIR}/control_frequent_snps.tsv" | sort > "${BASE_DIR}/control_snps_sorted.txt"
    cut -f1 "${BASE_DIR}/mdd_frequent_snps.tsv" | sort > "${BASE_DIR}/mdd_snps_sorted.txt"

    # Use comm to find unique SNPs
    comm -23 "${BASE_DIR}/control_snps_sorted.txt" "${BASE_DIR}/mdd_snps_sorted.txt" > "${BASE_DIR}/control_unique.txt"
    comm -13 "${BASE_DIR}/control_snps_sorted.txt" "${BASE_DIR}/mdd_snps_sorted.txt" > "${BASE_DIR}/mdd_unique.txt"

    # Get top N from each
    head -20 "${BASE_DIR}/control_unique.txt" > "${BASE_DIR}/top_20_control_unique.txt"
    head -20 "${BASE_DIR}/mdd_unique.txt" > "${BASE_DIR}/top_20_mdd_unique.txt"
}

# Step 4: Minimal annotation
annotate_top_snps() {
    # Stream through TSV files once to annotate top SNPs
    while read -r snp_key; do
        # Parse chrom:pos:ref:alt
        IFS=':' read -r chrom pos ref alt <<< "${snp_key}"

        # Find in any TSV file (first occurrence)
        find "${BASE_DIR}/data/tsv" -name "*.tsv.gz" -exec \
            sh -c "gunzip -c {} 2>/dev/null | awk -v chrom='$chrom' -v pos='$pos' -v ref='$ref' -v alt='$alt' '\$1==chrom && \$2==pos && \$4==ref && \$5==alt {print \$0; exit}'" \; | \
            head -1
    done < "${BASE_DIR}/top_20_control_unique.txt" > "${BASE_DIR}/annotated_control_unique.tsv"

    # Same for MDD
    while read -r snp_key; do
        IFS=':' read -r chrom pos ref alt <<< "${snp_key}"
        find "${BASE_DIR}/data/tsv" -name "*.tsv.gz" -exec \
            sh -c "gunzip -c {} 2>/dev/null | awk -v chrom='$chrom' -v pos='$pos' -v ref='$ref' -v alt='$alt' '\$1==chrom && \$2==pos && \$4==ref && \$5==alt {print \$0; exit}'" \; | \
            head -1
    done < "${BASE_DIR}/top_20_mdd_unique.txt" > "${BASE_DIR}/annotated_mdd_unique.tsv"
}
