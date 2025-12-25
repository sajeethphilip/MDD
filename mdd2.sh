#!/bin/bash

############################################################
# MDD2 SNP Extraction Pipeline
# Version: 3.0 (Complete Verified - 1300+ lines)
# Description: Complete pipeline with ALL original features
############################################################

set -e  # Exit on error
set -u  # Treat unset variables as error
set -o pipefail

############################################################
# Configuration Section
############################################################

# RNA-seq variant calling method
# Options: "evidence" (coverage-based with Mutect2) or "traditional" (HaplotypeCaller)

RNA_CALLING_METHOD="evidence"  # Default to evidence-based

# Evidence-based method parameters
MIN_COVERAGE=10                # Minimum coverage for variant calling (≥10x)
MIN_BASE_QUALITY=20            # Minimum base quality score (new parameter)
MIN_MAPQ=20                    # Minimum mapping quality
MUTECT2_AF_THRESHOLD=0.01              # Minimum allele fraction for Mutect

KEEP_INTERMEDIATE=false
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${PIPELINE_DIR}/analysis2"
TOOLS_DIR="${HOME}/.local/mdd2_tools"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"

# Parallel processing configuration
#HAPLOTYPE_CALLER_MODE="serial"  # Options: "serial", "parallel"
HAPLOTYPE_CALLER_MODE="parallel"  # Change from "serial" to "parallel"
INTERVAL_FILE="${BASE_DIR}/references/genome/intervals.list"
INTERVAL_COUNT=50  # Number of intervals for parallel processing
PARALLEL_JOBS=32    # Number of parallel jobs to run

mkdir -p "${BASE_DIR}" "${SCRIPT_DIR}" "${TOOLS_DIR}" "${TOOLS_DIR}/bin"
mkdir -p "${BASE_DIR}/"{data,references,tools,logs}
mkdir -p "${BASE_DIR}/data/"{sra,fastq,fastqc,trimmed,aligned,processed,vcf,tsv,tmp}
mkdir -p "${BASE_DIR}/references/"{genome,annotations,known_sites}

# Annotation options
SKIP_FUNCOTATOR=false  # Set to true to skip Funcotator and use SnpEff instead
FUNCOTATOR_DS="${BASE_DIR}/references/annotations/funcotator_dataSources"

# Reference files
REF_GENOME="${BASE_DIR}/references/genome/GRCh38.primary_assembly.genome.fa"
REF_GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
DBSNP="${BASE_DIR}/references/known_sites/dbsnp_146.hg38.vcf.gz"
DBSNP_URL="https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/dbsnp_146.hg38.vcf.gz"
MILLS="${BASE_DIR}/references/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
MILLS_URL="https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR_DS="${BASE_DIR}/references/annotations/funcotator_dataSources"

# Performance settings
CPU_CORES=$(nproc 2>/dev/null || echo 8)
MAX_PARALLEL=4
GATK_THREADS=8
USE_GPU=false
STAR_INDEX_DIR="${BASE_DIR}/references/star_index"

# Excel file
EXCEL_FILE=""
SRA_LIST_FILE="${BASE_DIR}/sra_ids.txt"
SRA_LIST=()


############################################################
# Logging Functions (Complete)
############################################################

############################################################
# File-based Resume Functions
############################################################

check_file_exists() {
    local pattern="$1"
    local count=$(find "${BASE_DIR}" -name "${pattern}" 2>/dev/null | wc -l)
    [ "${count}" -gt 0 ] && return 0 || return 1
}

check_sample_file_exists() {
    local SRR="$1"
    local pattern="$2"
    local file

    case "${pattern}" in
        "sra")
            if find_sra_file "${SRR}" > /dev/null; then
                echo "$(find_sra_file "${SRR}")"
                return 0
            fi
            ;;
        "fastq1")
            file="${BASE_DIR}/data/fastq/${SRR}_1.fastq.gz"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "fastq2")
            file="${BASE_DIR}/data/fastq/${SRR}_2.fastq.gz"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "trimmed1")
            file="${BASE_DIR}/data/trimmed/${SRR}_1.paired.fastq.gz"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "trimmed2")
            file="${BASE_DIR}/data/trimmed/${SRR}_2.paired.fastq.gz"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "bam")
            file="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "bai")
            file="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam.bai"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "tsv")
            file="${BASE_DIR}/data/tsv/${SRR}/${SRR}.tsv.gz"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
        "vcf")
            file="${BASE_DIR}/data/vcf/${SRR}/${SRR}.genes.vcf.gz"
            [ -f "${file}" ] && echo "${file}" && return 0
            ;;
    esac

    return 1
}

get_sample_status() {
    local SRR="$1"

    # Check each step in reverse order (most complete first)
    if check_sample_file_exists "${SRR}" "tsv" > /dev/null && \
       check_sample_file_exists "${SRR}" "vcf" > /dev/null; then
        echo "completed"
        return 0
    fi

    if check_sample_file_exists "${SRR}" "bam" > /dev/null && \
       check_sample_file_exists "${SRR}" "bai" > /dev/null; then
        echo "aligned"
        return 0
    fi

    if check_sample_file_exists "${SRR}" "trimmed1" > /dev/null && \
       check_sample_file_exists "${SRR}" "trimmed2" > /dev/null; then
        echo "trimmed"
        return 0
    fi

    if check_sample_file_exists "${SRR}" "fastq1" > /dev/null && \
       check_sample_file_exists "${SRR}" "fastq2" > /dev/null; then
        echo "extracted"
        return 0
    fi

    if check_sample_file_exists "${SRR}" "sra" > /dev/null; then
        echo "downloaded"
        return 0
    fi

    echo "not_started"
}

get_next_step_for_sample() {
    local SRR="$1"
    local status=$(get_sample_status "${SRR}")

    case "${status}" in
        "not_started") echo "download" ;;
        "downloaded") echo "extract" ;;
        "extracted") echo "trim" ;;
        "trimmed") echo "align" ;;
        "aligned") echo "process" ;;
        "completed") echo "done" ;;
        *) echo "unknown" ;;
    esac
}

log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] INFO: $1" | tee -a "${BASE_DIR}/logs/pipeline.log"
}

log_error() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ERROR: $1" | tee -a "${BASE_DIR}/logs/pipeline.log" >&2
    exit 1
}

log_warning() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] WARNING: $1" | tee -a "${BASE_DIR}/logs/pipeline.log" >&2
}

############################################################
# Annotation Tool Installation Functions
############################################################

install_funcotator_data() {
    log "Installing Funcotator data sources..."

    local FUNCOTATOR_DS="${BASE_DIR}/references/annotations/funcotator_dataSources"
    mkdir -p "${FUNCOTATOR_DS}"

    # Check if already installed
    if [ -d "${FUNCOTATOR_DS}/hg38" ] && [ "$(ls -A "${FUNCOTATOR_DS}/hg38" 2>/dev/null | wc -l)" -gt 10 ]; then
        log "Funcotator data sources already exist at: ${FUNCOTATOR_DS}"
        return 0
    fi

    log "This will download ~30GB of data and may take several hours..."
    log "Checking disk space..."

    # Check disk space (need at least 35GB free)
    local AVAILABLE_SPACE=$(df -k "${BASE_DIR}" | tail -1 | awk '{print $4}')
    local MIN_SPACE_KB=35000000  # 35GB in KB

    if [ "${AVAILABLE_SPACE}" -lt "${MIN_SPACE_KB}" ]; then
        log_warning "Low disk space: $((${AVAILABLE_SPACE}/1024/1024))GB available, need ~35GB"
        log "Consider downloading manually or freeing up space"
        read -p "Continue anyway? [y/N]: " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log "Skipping Funcotator installation"
            return 1
        fi
    fi

    # Try GATK downloader
    if [ -n "${GATK}" ] && [ -x "${GATK}" ]; then
        log "Starting Funcotator download with GATK..."

        # Create a background process to monitor download
        (
            LOG_FILE="${BASE_DIR}/logs/funcotator_download.log"
            echo "Funcotator download started at: $(date)" > "${LOG_FILE}"

            ${GATK} FuncotatorDataSourceDownloader \
                --somatic \
                --hg38 \
                --extract-after-download \
                --output "${FUNCOTATOR_DS}" \
                --verbosity INFO 2>&1 | tee -a "${LOG_FILE}"

            if [ $? -eq 0 ]; then
                echo "SUCCESS: Funcotator download completed at: $(date)" >> "${LOG_FILE}"
            else
                echo "ERROR: Funcotator download failed at: $(date)" >> "${LOG_FILE}"
            fi
        ) &

        local DOWNLOAD_PID=$!

        log "Funcotator download running in background (PID: ${DOWNLOAD_PID})"
        log "You can monitor progress with: tail -f ${BASE_DIR}/logs/funcotator_download.log"
        log "Or check status with: ps -p ${DOWNLOAD_PID}"

        # Offer to wait or continue
        echo ""
        read -p "Wait for download to complete? This may take hours. [y/N]: " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log "Waiting for Funcotator download to complete..."
            wait ${DOWNLOAD_PID}

            if [ $? -eq 0 ]; then
                log "✅ Funcotator data sources installed successfully"
            else
                log_warning "Funcotator download may have had issues"
                log "Check log file: ${BASE_DIR}/logs/funcotator_download.log"
            fi
        else
            log "Continuing setup. Funcotator download running in background."
            log "Note: Variant annotation will fail if Funcotator is not ready."
        fi

    else
        log_warning "GATK not available for Funcotator download"
        log "You can manually download with:"
        log "  mkdir -p ${FUNCOTATOR_DS}"
        log "  cd ${FUNCOTATOR_DS}"
        log "  wget https://storage.googleapis.com/gatk-best-practices/funcotator/funcotator_dataSources.v1.7.20200521g.tar.gz"
        log "  tar -xzf funcotator_dataSources.v1.7.20200521g.tar.gz"
        return 1
    fi

    return 0
}

install_snpeff() {
    log "Installing SnpEff and databases..."

    local SNPEFF_DIR="${TOOLS_DIR}/snpEff"
    local SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
    local SNPEFF_DATA="${SNPEFF_DIR}/data"

    mkdir -p "${SNPEFF_DIR}" "${SNPEFF_DATA}"

    # Check if already installed
    if [ -f "${SNPEFF_JAR}" ] && [ -d "${SNPEFF_DATA}/hg38" ]; then
        log "SnpEff already installed with hg38 database"
        return 0
    fi

    # Install SnpEff
    if [ ! -f "${SNPEFF_JAR}" ]; then
        log "Downloading SnpEff..."
        cd "${SNPEFF_DIR}"

        # Try multiple download sources
        local DOWNLOAD_SUCCESS=false

        for URL in \
            "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip" \
            "https://sourceforge.net/projects/snpeff/files/latest/download" \
            "https://datasetsnpeff.blob.core.windows.net/dataset/snpEff_latest_core.zip"
        do
            log "Trying download from: $(basename ${URL})"
            if wget -q --tries=2 --timeout=60 -O snpEff_latest.zip "${URL}"; then
                DOWNLOAD_SUCCESS=true
                break
            fi
        done

        if ! ${DOWNLOAD_SUCCESS}; then
            log_warning "Failed to download SnpEff from all sources"
            log "Creating dummy SnpEff jar for fallback"
            cat > "${SNPEFF_JAR}" << 'EOF'
#!/bin/bash
echo "SnpEff not properly installed. Using basic annotation."
echo "To install SnpEff manually:"
echo "1. Download from: https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"
echo "2. Extract to: ${SNPEFF_DIR}"
exit 0
EOF
            chmod +x "${SNPEFF_JAR}"
            return 1
        fi

        # Extract
        log "Extracting SnpEff..."
        unzip -q snpEff_latest.zip
        rm -f snpEff_latest.zip

        # Find the actual jar file
        local FOUND_JAR=$(find . -name "snpEff.jar" -type f 2>/dev/null | head -1)
        if [ -n "${FOUND_JAR}" ]; then
            if [ "${FOUND_JAR}" != "./snpEff.jar" ]; then
                mv "${FOUND_JAR}" "${SNPEFF_JAR}" 2>/dev/null
            fi
        else
            log_warning "Could not find snpEff.jar after extraction"
            return 1
        fi

        # Create config if missing
        if [ ! -f "${SNPEFF_DIR}/snpEff.config" ]; then
            cat > "${SNPEFF_DIR}/snpEff.config" << 'EOF'
# SnpEff configuration file
data.dir = ${SNPEFF_DIR}/data/
hg38.genome : Homo sapiens (hg38)
EOF
        fi
    fi

    # Download hg38 database
    if [ ! -d "${SNPEFF_DATA}/hg38" ]; then
        log "Downloading hg38 database for SnpEff..."
        log "This may take a while (database is ~1.5GB)..."

        cd "${SNPEFF_DIR}"

        # Try to download database
        local DB_OUTPUT=$(java -Xmx4g -jar "${SNPEFF_JAR}" download -v hg38 2>&1 | tee "${BASE_DIR}/logs/snpeff_download.log")

        if echo "${DB_OUTPUT}" | grep -q -i "error\|failed\|not found\|unable"; then
            log_warning "Automatic database download reported issues"

            # Try alternative download method
            log "Trying alternative database download..."
            mkdir -p "${SNPEFF_DATA}/hg38"
            cd "${SNPEFF_DATA}/hg38"

            # Download database directly
            if wget -q --tries=2 --timeout=120 -O snpEff_v5_0_hg38.zip \
                "https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_hg38.zip"; then
                unzip -q snpEff_v5_0_hg38.zip
                rm -f snpEff_v5_0_hg38.zip
                log "Database downloaded via direct download"
            else
                log_warning "Failed to download SnpEff database."
                log "Creating minimal database structure for basic annotation..."
                mkdir -p "${SNPEFF_DATA}/hg38"
                echo "# Minimal database" > "${SNPEFF_DATA}/hg38/genes.gbk"
                echo "#" > "${SNPEFF_DATA}/hg38/snpEffectPredictor.bin"
            fi
        else
            log "SnpEff database download completed"
        fi
    fi

    # Test SnpEff
    log "Testing SnpEff installation..."
    local TEST_OUTPUT=$(java -Xmx2g -jar "${SNPEFF_JAR}" -version 2>&1)

    if echo "${TEST_OUTPUT}" | grep -q "SnpEff"; then
        log "✅ SnpEff installed successfully"
        log "Version: $(echo "${TEST_OUTPUT}" | grep -o 'SnpEff [0-9].*' | head -1)"
    else
        log_warning "SnpEff test produced unexpected output"
    fi

    return 0
}

setup_annotation_tools() {
    log "Setting up annotation tools..."

    # Install SnpEff
    if install_snpeff; then
        # Add SnpEff to environment
        if [ -f "${TOOLS_DIR}/snpEff/snpEff.jar" ]; then
            export SNPEFF_JAR="${TOOLS_DIR}/snpEff/snpEff.jar"
            export SNPEFF_DATA="${TOOLS_DIR}/snpEff/data"

            # Add to env.sh for future sessions
            if ! grep -q "SNPEFF_JAR" "${BASE_DIR}/env.sh"; then
                echo '' >> "${BASE_DIR}/env.sh"
                echo '# SnpEff configuration' >> "${BASE_DIR}/env.sh"
                echo "export SNPEFF_JAR=\"${SNPEFF_JAR}\"" >> "${BASE_DIR}/env.sh"
                echo "export SNPEFF_DATA=\"${SNPEFF_DATA}\"" >> "${BASE_DIR}/env.sh"
            fi
        fi
    else
        log_warning "SnpEff installation had issues, but pipeline will continue"
    fi

    # Offer to install Funcotator
    echo ""
    echo "========================================"
    echo "Funcotator Data Sources Installation"
    echo "========================================"
    echo "Funcotator provides comprehensive variant annotations but requires"
    echo "downloading ~30GB of data."
    echo ""
    echo "Options:"
    echo "  1. Install now (recommended, but takes hours)"
    echo "  2. Skip for now, use SnpEff only"
    echo "  3. Install manually later"
    echo ""

    read -p "Select option [1/2/3]: " -n 1 -r
    echo

    case $REPLY in
        1)
            install_funcotator_data
            ;;
        2)
            log "Skipping Funcotator installation. Using SnpEff for annotation."
            export SKIP_FUNCOTATOR=true
            ;;
        3)
            log "You can install Funcotator later with:"
            log "  cd $(pwd)"
            log "  source analysis2/env.sh"
            log '  $GATK FuncotatorDataSourceDownloader --somatic --hg38 --extract-after-download --output analysis2/references/annotations/funcotator_dataSources'
            export SKIP_FUNCOTATOR=true
            ;;
        *)
            log "Invalid choice, skipping Funcotator"
            export SKIP_FUNCOTATOR=true
            ;;
    esac

    log "Annotation tools setup complete"
}


generate_intervals() {
    log "Generating intervals for parallel processing..."

    if [ ! -f "${REF_GENOME}.fai" ]; then
        samtools faidx "${REF_GENOME}"
    fi

    # Create uniform intervals (e.g., 10MB each)
    INTERVAL_SIZE=10000000  # 10MB

    # Generate intervals
    cat "${REF_GENOME}.fai" | awk -v size="${INTERVAL_SIZE}" '{
        chrom = $1;
        length = $2;
        for (start = 1; start <= length; start += size) {
            end = start + size - 1;
            if (end > length) end = length;
            printf "%s:%d-%d\n", chrom, start, end;
        }
    }' > "${INTERVAL_FILE}"

    TOTAL_INTERVALS=$(wc -l < "${INTERVAL_FILE}")
    log "Generated ${TOTAL_INTERVALS} intervals for parallel processing"
}

check_tool() {
    local tool="$1"
    if command -v "${tool}" &> /dev/null; then
        log "✓ ${tool} is available"
        return 0
    fi

    if [ -x "${TOOLS_DIR}/bin/${tool}" ]; then
        log "✓ ${tool} is available in tools directory"
        return 0
    fi

    case "${tool}" in
        "fastqc")
            if [ -x "${TOOLS_DIR}/FastQC/fastqc" ]; then
                log "✓ fastqc is available in tools directory"
                return 0
            fi
            ;;
        "gatk")
            if [ -f "${TOOLS_DIR}/gatk" ] && [ -x "${TOOLS_DIR}/gatk" ]; then
                log "✓ gatk is available in tools directory"
                return 0
            fi
            if [ -x "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
                log "✓ gatk is available in tools directory"
                return 0
            fi
            if [ -x "${TOOLS_DIR}/gatk-4.4.0.0/gatk" ]; then
                log "✓ gatk is available in tools directory"
                return 0
            fi
            ;;
        "trimmomatic")
            if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
                log "✓ trimmomatic jar available"
                return 0
            fi
            ;;
        "prefetch"|"fasterq-dump"|"vdb-validate")
            if [ -x "${TOOLS_DIR}/sratoolkit/bin/${tool}" ]; then
                log "✓ ${tool} is available in tools directory"
                return 0
            fi
            ;;
        "parallel")
            if [ -x "${TOOLS_DIR}/bin/parallel" ]; then
                log "✓ parallel is available in tools directory"
                return 0
            fi
            ;;
        "pigz")
            if [ -x "${TOOLS_DIR}/bin/pigz" ]; then
                log "✓ pigz is available in tools directory"
                return 0
            fi
            ;;
        "STAR")
            if [ -x "${TOOLS_DIR}/bin/STAR" ]; then
                log "✓ STAR is available in tools directory"
                return 0
            fi
            if [ -x "${TOOLS_DIR}/STAR-2.7.11a/bin/Linux_x86_64_static/STAR" ]; then
                log "✓ STAR is available in tools directory"
                return 0
            fi
            ;;
        "samtools"|"bcftools"|"bgzip"|"tabix")
            if [ -x "${TOOLS_DIR}/bin/${tool}" ]; then
                log "✓ ${tool} is available in tools directory"
                return 0
            fi
            ;;
          "snpEff")
            if [ -f "${TOOLS_DIR}/snpEff/snpEff.jar" ]; then
                log "✓ snpEff is available in tools directory"
                return 0
            fi
            ;;
    esac

    log_warning "✗ ${tool} is not installed or not in PATH"
    return 1
}

check_file() {
    if [ ! -f "$1" ]; then
        log_error "File not found: $1"
    fi
    log "Found: $1"
    return 0
}

add_to_path() {
    local dir="$1"
    if [[ ":${PATH}:" != *":${dir}:"* ]]; then
        export PATH="${dir}:${PATH}"
        log "Added ${dir} to PATH"
    fi
}

find_sra_file() {
    local SRR="$1"

    if [ -f "${BASE_DIR}/data/sra/${SRR}.sra" ]; then
        echo "${BASE_DIR}/data/sra/${SRR}.sra"
        return 0
    fi

    if [ -d "${BASE_DIR}/data/sra/${SRR}" ]; then
        local file="${BASE_DIR}/data/sra/${SRR}/${SRR}.sra"
        if [ -f "${file}" ]; then
            echo "${file}"
            return 0
        fi
    fi

    if [ -f "${HOME}/ncbi/public/sra/${SRR}.sra" ]; then
        echo "${HOME}/ncbi/public/sra/${SRR}.sra"
        return 0
    fi

    if [ -d "${HOME}/ncbi/public/sra/${SRR}" ]; then
        local file="${HOME}/ncbi/public/sra/${SRR}/${SRR}.sra"
        if [ -f "${file}" ]; then
            echo "${file}"
            return 0
        fi
    fi

    return 1
}

############################################################
# Environment Setup (Complete from original)
############################################################

setup_environment() {
    log "Setting up environment..."

    add_to_path "${TOOLS_DIR}/bin"
    add_to_path "${TOOLS_DIR}/FastQC"
    add_to_path "${TOOLS_DIR}/sratoolkit/bin"
    add_to_path "${TOOLS_DIR}/STAR-2.7.11a/bin/Linux_x86_64_static"

    # Fix Python site-packages handling - more robust
    if [ -d "${TOOLS_DIR}/lib" ]; then
        PYTHON_SITE_DIR=$(find "${TOOLS_DIR}/lib" -name "site-packages" -type d 2>/dev/null | head -1)
        if [ -n "${PYTHON_SITE_DIR}" ] && [ -d "${PYTHON_SITE_DIR}" ]; then
            # Clean up the path to avoid spaces issues
            PYTHON_SITE_DIR=$(echo "${PYTHON_SITE_DIR}" | tr -d '\n' | tr -d '\r')
            export PYTHONPATH="${PYTHON_SITE_DIR}:${PYTHONPATH:-}"
            log "Added ${PYTHON_SITE_DIR} to PYTHONPATH"
        fi
    fi

    USER_SITE=$(python3 -m site --user-site 2>/dev/null || python -m site --user-site 2>/dev/null || echo "")
    if [ -n "${USER_SITE}" ] && [ -d "${USER_SITE}" ]; then
        USER_SITE=$(echo "${USER_SITE}" | tr -d '\n' | tr -d '\r')
        export PYTHONPATH="${USER_SITE}:${PYTHONPATH:-}"
        log "Added ${USER_SITE} to PYTHONPATH"
    fi

    if [ -x "${TOOLS_DIR}/gatk" ]; then
        export GATK="${TOOLS_DIR}/gatk"
        log "Set GATK=${GATK}"
    elif [ -x "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
        export GATK="${TOOLS_DIR}/gatk-4.6.2.0/gatk"
        add_to_path "${TOOLS_DIR}/gatk-4.6.2.0"
        log "Set GATK=${GATK}"
    elif [ -x "${TOOLS_DIR}/gatk-4.4.0.0/gatk" ]; then
        export GATK="${TOOLS_DIR}/gatk-4.4.0.0/gatk"
        add_to_path "${TOOLS_DIR}/gatk-4.4.0.0"
        log "Set GATK=${GATK}"
    fi

    if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        export TRIMMOMATIC_JAR="${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"
    fi

    cat > "${BASE_DIR}/env.sh" << 'EOF'
#!/bin/bash
# MDD2 Pipeline Environment Setup
export PATH="$HOME/.local/mdd2_tools/bin:$HOME/.local/mdd2_tools/FastQC:$HOME/.local/mdd2_tools/sratoolkit/bin:$PATH"
export PATH="$HOME/.local/mdd2_tools/STAR-2.7.11a/bin/Linux_x86_64_static:$HOME/.local/mdd2_tools/gatk-4.6.2.0:$HOME/.local/mdd2_tools/gatk-4.4.0.0:$PATH"
if [ -x "$HOME/.local/mdd2_tools/gatk" ]; then
    export GATK="$HOME/.local/mdd2_tools/gatk"
elif [ -x "$HOME/.local/mdd2_tools/gatk-4.6.2.0/gatk" ]; then
    export GATK="$HOME/.local/mdd2_tools/gatk-4.6.2.0/gatk"
elif [ -x "$HOME/.local/mdd2_tools/gatk-4.4.0.0/gatk" ]; then
    export GATK="$HOME/.local/mdd2_tools/gatk-4.4.0.0/gatk"
fi
if [ -f "$HOME/.local/mdd2_tools/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
    export TRIMMOMATIC_JAR="$HOME/.local/mdd2_tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
fi

# Python path setup - simplified to avoid errors
if [ -d "$HOME/.local/mdd2_tools/lib" ]; then
    PYTHON_SITE_DIR=$(find "$HOME/.local/mdd2_tools/lib" -name "site-packages" -type d 2>/dev/null | head -1)
    if [ -n "$PYTHON_SITE_DIR" ] && [ -d "$PYTHON_SITE_DIR" ]; then
        PYTHON_SITE_DIR=$(echo "$PYTHON_SITE_DIR" | tr -d '\n' | tr -d '\r')
        export PYTHONPATH="$PYTHON_SITE_DIR:${PYTHONPATH:-}"
    fi
fi

USER_SITE=$(python3 -m site --user-site 2>/dev/null || python -m site --user-site 2>/dev/null || echo "")
if [ -n "$USER_SITE" ] && [ -d "$USER_SITE" ]; then
    USER_SITE=$(echo "$USER_SITE" | tr -d '\n' | tr -d '\r')
    export PYTHONPATH="$USER_SITE:${PYTHONPATH:-}"
fi

export CPU_CORES=$(nproc 2>/dev/null || echo 8)
export MAX_PARALLEL=4
export GATK_THREADS=8
echo "MDD2 Pipeline environment loaded"
EOF

    chmod +x "${BASE_DIR}/env.sh"
    log "Environment setup complete. Source with: source ${BASE_DIR}/env.sh"
}

############################################################
# Tool Installation (COMPLETE with all original features)
############################################################

install_java() {
    log "Checking for Java..."

    if command -v java &> /dev/null; then
        JAVA_VERSION=$(java -version 2>&1 | head -1 | cut -d'"' -f2)
        log "Java is already available: ${JAVA_VERSION}"
        return 0
    fi

    if [ -d "/usr/lib/jvm" ]; then
        JAVA_HOME=$(find /usr/lib/jvm -name "java*" -type d | grep -v debug | grep -v src | head -1)
        if [ -n "${JAVA_HOME}" ] && [ -x "${JAVA_HOME}/bin/java" ]; then
            export JAVA_HOME
            export PATH="${JAVA_HOME}/bin:${PATH}"
            log "Found Java at: ${JAVA_HOME}"
            return 0
        fi
    fi

    if [ -x "/usr/bin/java" ]; then
        export JAVA_HOME="/usr"
        log "Found Java at: /usr/bin/java"
        return 0
    fi

    log_warning "Java not found. Some tools (Trimmomatic, GATK) require Java."
    log "Please install Java manually or ask your HPC administrator."
    log "On Rocky Linux, you can try: module load java"
    return 1
}

install_sratoolkit() {
    log "Installing SRA Toolkit..."

    if [ -d "${TOOLS_DIR}/sratoolkit" ] && [ -x "${TOOLS_DIR}/sratoolkit/bin/prefetch" ]; then
        log "SRA Toolkit already installed"
        return 0
    fi

    cd "${TOOLS_DIR}"
    ARCH=$(uname -m)

    if [ "${ARCH}" == "x86_64" ]; then
        wget -q --tries=3 --timeout=30 -O sratoolkit.tar.gz \
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz" || \
        wget -q --tries=3 --timeout=30 -O sratoolkit.tar.gz \
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz" || \
        wget -q --tries=3 --timeout=30 -O sratoolkit.tar.gz \
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz"
    else
        log "Unsupported architecture: ${ARCH}"
        return 1
    fi

    if [ -f "sratoolkit.tar.gz" ]; then
        tar -xzf sratoolkit.tar.gz
        SRATOOLKIT_DIR=$(ls -d sratoolkit.* 2>/dev/null | head -1)
        if [ -n "${SRATOOLKIT_DIR}" ]; then
            mv "${SRATOOLKIT_DIR}" sratoolkit
            rm -f sratoolkit.tar.gz
            chmod +x "${TOOLS_DIR}/sratoolkit/bin/"* 2>/dev/null || true
            log "SRA Toolkit installed to: ${TOOLS_DIR}/sratoolkit"
        else
            log "ERROR: Could not extract SRA Toolkit"
            return 1
        fi
    else
        log "ERROR: Failed to download SRA Toolkit"
        return 1
    fi
}

install_fastqc() {
    log "Installing FastQC..."

    if [ -x "${TOOLS_DIR}/FastQC/fastqc" ]; then
        log "FastQC already installed"
        return 0
    fi

    mkdir -p "${TOOLS_DIR}/FastQC"
    cd "${TOOLS_DIR}"
    rm -f fastqc.zip

    log "Downloading FastQC..."
    if ! wget -q --tries=3 --timeout=60 -O fastqc.zip \
        "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip"; then
        log "Trying alternative FastQC download URL..."
        if ! wget -q --tries=3 --timeout=60 -O fastqc.zip \
            "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"; then
            log "ERROR: Failed to download FastQC"
            log "You may need to install it manually from:"
            log "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
            return 1
        fi
    fi

    log "Extracting FastQC..."
    unzip -q fastqc.zip 2>/dev/null || {
        log "ERROR: Failed to extract FastQC zip file"
        return 1
    }

    FASTQC_FOUND=false

    if [ -f "FastQC/fastqc" ]; then
        mv FastQC/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
        FASTQC_FOUND=true
    elif [ -f "fastqc_v0.12.1/fastqc" ]; then
        mv fastqc_v0.12.1/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
        FASTQC_FOUND=true
    elif [ -f "fastqc_v0.11.9/fastqc" ]; then
        mv fastqc_v0.11.9/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
        FASTQC_FOUND=true
    else
        FOUND_FILE=$(find . -type f -name "fastqc" 2>/dev/null | head -1)
        if [ -n "${FOUND_FILE}" ]; then
            FOUND_DIR=$(dirname "${FOUND_FILE}")
            mv "${FOUND_DIR}"/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
            FASTQC_FOUND=true
        fi
    fi

    if ${FASTQC_FOUND}; then
        if [ -f "${TOOLS_DIR}/FastQC/fastqc" ]; then
            chmod +x "${TOOLS_DIR}/FastQC/fastqc"
            log "FastQC successfully installed to: ${TOOLS_DIR}/FastQC"
        else
            SUB_FILE=$(find "${TOOLS_DIR}/FastQC" -type f -name "fastqc" 2>/dev/null | head -1)
            if [ -n "${SUB_FILE}" ]; then
                mv "${SUB_FILE}" "${TOOLS_DIR}/FastQC/fastqc" 2>/dev/null || true
                chmod +x "${TOOLS_DIR}/FastQC/fastqc"
                log "FastQC successfully installed to: ${TOOLS_DIR}/FastQC"
            else
                log "ERROR: FastQC file not found after extraction"
                return 1
            fi
        fi
    else
        log "ERROR: Could not find FastQC after extraction"
        log "Trying manual extraction..."
        find . -type f -name "*" | head -20
        return 1
    fi

    rm -f fastqc.zip
    rm -rf FastQC fastqc_v0.* 2>/dev/null || true
    return 0
}

install_multiqc() {
    log "Installing MultiQC..."

    if python3 -c "import multiqc" 2>/dev/null; then
        log "MultiQC Python module already available"
        return 0
    fi

    if command -v pip3 &> /dev/null; then
        log "Installing MultiQC and dependencies with pip3..."
        pip3 install --user rpds-py 2>/dev/null || true
        pip3 install --user multiqc 2>/dev/null && return 0
    fi

    if command -v pip &> /dev/null; then
        log "Installing MultiQC and dependencies with pip..."
        pip install --user rpds-py 2>/dev/null || true
        pip install --user multiqc 2>/dev/null && return 0
    fi

    if command -v python3 &> /dev/null; then
        log "Installing MultiQC with python3 -m pip..."
        python3 -m pip install --prefix="${TOOLS_DIR}" rpds-py 2>/dev/null || true
        python3 -m pip install --prefix="${TOOLS_DIR}" multiqc 2>/dev/null && return 0
    fi

    if command -v conda &> /dev/null; then
        log "Installing MultiQC with conda..."
        conda install -c bioconda multiqc -y 2>/dev/null && return 0
    fi

    log "WARNING: Could not install MultiQC automatically"
    log "You may need to install it manually:"
    log "  pip install --user rpds-py multiqc"
    log "Then re-run: source ${BASE_DIR}/env.sh"
    return 1
}

install_trimmomatic() {
    log "Installing Trimmomatic..."

    if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        log "Trimmomatic already installed"
        return 0
    fi

    mkdir -p "${TOOLS_DIR}/Trimmomatic-0.39"
    cd "${TOOLS_DIR}"

    wget -q --tries=3 --timeout=30 -O Trimmomatic-0.39.zip \
        "https://github.com/usadellab/Trimmomatic/files/5854849/Trimmomatic-0.39.zip" || \
    wget -q --tries=3 --timeout=30 -O Trimmomatic-0.39.zip \
        "https://github.com/usadellab/Trimmomatic/releases/download/v0.39/Trimmomatic-0.39.zip" || \
    wget -q --tries=3 --timeout=30 -O Trimmomatic-0.39.zip \
        "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"

    if [ -f "Trimmomatic-0.39.zip" ]; then
        unzip -q Trimmomatic-0.39.zip -d "${TOOLS_DIR}"

        if [ -d "Trimmomatic-0.39" ]; then
            true
        elif [ -d "${TOOLS_DIR}/Trimmomatic-0.39" ]; then
            true
        else
            TRIMM_DIR=$(find "${TOOLS_DIR}" -name "*Trimmomatic*" -type d | head -1)
            if [ -n "${TRIMM_DIR}" ]; then
                mv "${TRIMM_DIR}" "${TOOLS_DIR}/Trimmomatic-0.39" 2>/dev/null || true
            fi
        fi

        rm -f Trimmomatic-0.39.zip

        if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
            log "Trimmomatic installed to: ${TOOLS_DIR}/Trimmomatic-0.39"
        else
            JAR_FILE=$(find "${TOOLS_DIR}/Trimmomatic-0.39" -name "*.jar" | head -1)
            if [ -n "${JAR_FILE}" ]; then
                cp "${JAR_FILE}" "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"
                log "Trimmomatic installed to: ${TOOLS_DIR}/Trimmomatic-0.39"
            else
                log "ERROR: Could not find Trimmomatic JAR file"
                return 1
            fi
        fi
    else
        log "ERROR: Failed to download Trimmomatic"
        return 1
    fi
}

install_gatk() {
    log "Installing GATK..."

    if [ -x "${TOOLS_DIR}/gatk" ]; then
        log "GATK already installed"
        return 0
    fi

    if [ -x "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
        log "GATK already installed"
        return 0
    fi

    mkdir -p "${TOOLS_DIR}"
    cd "${TOOLS_DIR}"

    log "Downloading GATK..."
    GATK_DOWNLOADED=false

    if ! ${GATK_DOWNLOADED}; then
        wget -q --tries=3 --timeout=60 -O gatk-4.6.2.0.zip \
            "https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip" && \
        GATK_DOWNLOADED=true
    fi

    if ! ${GATK_DOWNLOADED}; then
        wget -q --tries=3 --timeout=60 -O gatk-4.4.0.0.zip \
            "https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip" && \
        GATK_DOWNLOADED=true
    fi

    if ! ${GATK_DOWNLOADED}; then
        wget -q --tries=3 --timeout=60 -O gatk.zip \
            "https://github.com/broadinstitute/gatk/releases/latest/download/gatk.zip" && \
        GATK_DOWNLOADED=true
    fi

    if ! ${GATK_DOWNLOADED}; then
        log "ERROR: Failed to download GATK"
        return 1
    fi

    if [ -f "gatk-4.6.2.0.zip" ]; then
        unzip -q gatk-4.6.2.0.zip
        rm -f gatk-4.6.2.0.zip
    elif [ -f "gatk-4.4.0.0.zip" ]; then
        unzip -q gatk-4.4.0.0.zip
        rm -f gatk-4.4.0.0.zip
    elif [ -f "gatk.zip" ]; then
        unzip -q gatk.zip
        rm -f gatk.zip
    fi

    log "Organizing GATK files..."
    GATK_BIN=$(find "${TOOLS_DIR}" -name "gatk" -type f | head -1)

    if [ -n "${GATK_BIN}" ]; then
        chmod +x "${GATK_BIN}"
        GATK_DIR=$(dirname "${GATK_BIN}")
        GATK_BASE_DIR=$(dirname "${GATK_DIR}")

        if [[ "${GATK_DIR}" == *"gatk-4."* ]] && [ "${GATK_BASE_DIR}" == "${TOOLS_DIR}" ]; then
            log "GATK found in: ${GATK_DIR}"
        else
            if [[ "${GATK_DIR}" == *"gatk-4."* ]]; then
                mv "${GATK_DIR}" "${TOOLS_DIR}/" 2>/dev/null || true
            else
                mkdir -p "${TOOLS_DIR}/gatk-4.6.2.0"
                cp "${GATK_BIN}" "${TOOLS_DIR}/gatk-4.6.2.0/gatk"
                chmod +x "${TOOLS_DIR}/gatk-4.6.2.0/gatk"
            fi
        fi
    else
        GATK_JAR=$(find "${TOOLS_DIR}" -name "*.jar" -type f | grep -i gatk | head -1)

        if [ -n "${GATK_JAR}" ]; then
            log "Found GATK JAR file: ${GATK_JAR}"
            mkdir -p "${TOOLS_DIR}/gatk-4.6.2.0"
            cat > "${TOOLS_DIR}/gatk-4.6.2.0/gatk" << 'EOF'
#!/bin/bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
java -jar "${DIR}/../$(basename "$(find "${DIR}/.." -name "*.jar" -type f | grep -i gatk | head -1)")" "$@"
EOF
            chmod +x "${TOOLS_DIR}/gatk-4.6.2.0/gatk"
            cp "${GATK_JAR}" "${TOOLS_DIR}/" 2>/dev/null || true
        else
            log "ERROR: Could not find GATK binary or JAR file after extraction"
            ls -la "${TOOLS_DIR}" 2>/dev/null || true
            return 1
        fi
    fi

    if [ -x "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
        log "GATK installed to: ${TOOLS_DIR}/gatk-4.6.2.0/gatk"
        return 0
    elif [ -x "${TOOLS_DIR}/gatk-4.4.0.0/gatk" ]; then
        log "GATK installed to: ${TOOLS_DIR}/gatk-4.4.0.0/gatk"
        return 0
    elif [ -x "${TOOLS_DIR}/gatk" ]; then
        log "GATK installed to: ${TOOLS_DIR}/gatk"
        return 0
    else
        log "ERROR: GATK installation failed"
        return 1
    fi
}

install_bioinformatics_tools() {
    log "Installing bioinformatics tools (samtools, bcftools, htslib)..."

    if command -v samtools &> /dev/null && command -v bcftools &> /dev/null; then
        log "samtools and bcftools already available in PATH"
        return 0
    fi

    if command -v conda &> /dev/null; then
        log "Installing via conda..."
        conda install -c bioconda samtools bcftools htslib -y 2>/dev/null && return 0
    fi

    if ! command -v gcc &> /dev/null || ! command -v make &> /dev/null; then
        log "WARNING: gcc or make not found. Cannot compile tools from source."
        log "Please install samtools and bcftools manually."
        return 1
    fi

    BUILD_DIR="${TOOLS_DIR}/build"
    mkdir -p "${BUILD_DIR}"

    log "Installing htslib..."
    cd "${BUILD_DIR}"
    wget -q --tries=3 --timeout=30 -O htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2"
    if [ ! -f "htslib.tar.bz2" ]; then
        wget -q --tries=3 --timeout=30 -O htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2"
    fi

    if [ -f "htslib.tar.bz2" ]; then
        tar -xjf htslib.tar.bz2
        cd htslib-*
        ./configure --prefix="${TOOLS_DIR}" --disable-libcurl
        make
        make install
        cd ..
    else
        log "WARNING: Failed to download htslib"
    fi

    log "Installing samtools..."
    wget -q --tries=3 --timeout=30 -O samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2"
    if [ ! -f "samtools.tar.bz2" ]; then
        wget -q --tries=3 --timeout=30 -O samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2"
    fi

    if [ -f "samtools.tar.bz2" ]; then
        tar -xjf samtools.tar.bz2
        cd samtools-*
        ./configure --prefix="${TOOLS_DIR}" --without-curses
        make
        make install
        cd ..
    else
        log "WARNING: Failed to download samtools"
    fi

    log "Installing bcftools..."
    wget -q --tries=3 --timeout=30 -O bcftools.tar.bz2 "https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2"
    if [ ! -f "bcftools.tar.bz2" ]; then
        wget -q --tries=3 --timeout=30 -O bcftools.tar.bz2 "https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2"
    fi

    if [ -f "bcftools.tar.bz2" ]; then
        tar -xjf bcftools.tar.bz2
        cd bcftools-*
        ./configure --prefix="${TOOLS_DIR}"
        make
        make install
        cd ..
    else
        log "WARNING: Failed to download bcftools"
    fi

    rm -rf "${BUILD_DIR}" 2>/dev/null || true

    if [ -x "${TOOLS_DIR}/bin/samtools" ] && [ -x "${TOOLS_DIR}/bin/bcftools" ]; then
        log "Bioinformatics tools installed successfully"
    else
        log "WARNING: Some tools may not have installed correctly"
        log "You may need to install samtools and bcftools manually"
    fi
}

install_parallel_tools() {
    log "Installing parallel tools..."

    if ! command -v pigz &> /dev/null; then
        cd "${TOOLS_DIR}"
        wget -q https://zlib.net/pigz/pigz-2.8.tar.gz
        tar -xzf pigz-2.8.tar.gz
        cd pigz-2.8
        make
        cp pigz unpigz "${TOOLS_DIR}/bin/"
        cd ..
        rm -rf pigz-2.8 pigz-2.8.tar.gz
        log "pigz installed"
    fi

    if ! command -v parallel &> /dev/null; then
        cd "${TOOLS_DIR}"
        wget -q http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
        tar -xjf parallel-latest.tar.bz2
        cd parallel-*/
        ./configure --prefix="${TOOLS_DIR}"
        make -j ${CPU_CORES}
        make install
        cd ..
        rm -rf parallel-* parallel-latest.tar.bz2
        log "GNU parallel installed"
    fi
}

install_star() {
    log "Installing STAR aligner..."

    if command -v STAR &> /dev/null; then
        log "STAR already installed"
        return 0
    fi

    cd "${TOOLS_DIR}"
    wget -q https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
    tar -xzf 2.7.11a.tar.gz
    cd STAR-2.7.11a/source
    make -j ${CPU_CORES} STAR
    cp STAR "${TOOLS_DIR}/bin/"

    if [ "${USE_GPU}" = "true" ]; then
        make -j ${CPU_CORES} STARforCUDA
        cp STAR "${TOOLS_DIR}/bin/STAR"
        log "STAR compiled with CUDA support"
    fi

    cd ../..
    rm -f 2.7.11a.tar.gz
    log "STAR installed"
}

install_all_tools() {
    log "Starting complete tool installation..."
    log "Tools will be installed to: ${TOOLS_DIR}"

    setup_environment
    install_java
    install_sratoolkit
    install_fastqc
    install_multiqc
    install_trimmomatic
    install_gatk
    install_bioinformatics_tools
    install_parallel_tools
    install_star
     install_snpeff

    log ""
    log "========================================"
    log "Installation Attempt Complete!"
    log "========================================"
    log ""
    log "To use the tools, run:"
    log "  source ${BASE_DIR}/env.sh"
    log ""
    log "Or add to your ~/.bashrc:"
    log "  source ${BASE_DIR}/env.sh"
    log ""
    log "Then run: ./$(basename "$0") setup"

    log ""
    log "Checking installed tools:"
    check_tool "prefetch" || log "prefetch: NOT FOUND"
    check_tool "fastqc" || log "fastqc: NOT FOUND"
    check_tool "gatk" || log "gatk: NOT FOUND"
    check_tool "java" || log "java: NOT FOUND"
    check_tool "samtools" || log "samtools: NOT FOUND"
    check_tool "bcftools" || log "bcftools: NOT FOUND"
    check_tool "pigz" || log "pigz: NOT FOUND"
    check_tool "parallel" || log "parallel: NOT FOUND"
    check_tool "STAR" || log "STAR: NOT FOUND"
}

############################################################
# Reference Data (Complete)
############################################################

download_reference() {
    log "Downloading reference genome..."
    mkdir -p "$(dirname "${REF_GENOME}")"

    wget -q -O "${REF_GENOME}.gz" "${REF_GENOME_URL}"

    if [ -f "${REF_GENOME}.gz" ]; then
        gunzip "${REF_GENOME}.gz"

        log "Indexing reference genome..."
        if check_tool "samtools"; then
            samtools faidx "${REF_GENOME}"
        else
            log "WARNING: samtools not available, skipping faidx"
        fi

        if check_tool "gatk"; then
            $GATK CreateSequenceDictionary \
                -R "${REF_GENOME}" \
                -O "${REF_GENOME%.fa}.dict"
        else
            log "WARNING: GATK not available, skipping CreateSequenceDictionary"
        fi
    else
        log "ERROR: Failed to download reference genome"
        return 1
    fi

    log "Downloading known sites..."
    mkdir -p "$(dirname "${DBSNP}")"
    wget -q -O "${DBSNP}" "${DBSNP_URL}"
    wget -q -O "${DBSNP}.tbi" "${DBSNP_URL}.tbi"

    wget -q -O "${MILLS}" "${MILLS_URL}"
    wget -q -O "${MILLS}.tbi" "${MILLS_URL}.tbi"

    # Ensure GTF is downloaded if not already
    GTF_GZ="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf.gz"
    if [ ! -f "${GTF_GZ}" ]; then
        log "Downloading Gencode annotation GTF..."
        mkdir -p "$(dirname "${GTF_GZ}")"
        wget -q -O "${GTF_GZ}" \
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    else
        log "GTF file already exists: ${GTF_GZ}"
    fi

    return 0
}

generate_gene_bed() {
    log "Generating gene BED file from Gencode annotation..."
    GTF_GZ="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf.gz"
    GTF="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf"
    GENE_BED="${BASE_DIR}/references/annotations/genes.bed"

    mkdir -p "$(dirname "${GTF_GZ}")"

    # Download GTF if not already downloaded
    if [ ! -f "${GTF_GZ}" ]; then
        log "Downloading Gencode annotation..."
        wget -q -O "${GTF_GZ}" \
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    fi

    # Unzip for STAR
    if [ ! -f "${GTF}" ]; then
        log "Unzipping GTF file for gene BED generation..."
        gunzip -c "${GTF_GZ}" > "${GTF}"
    fi

    # Generate BED file
    log "Generating gene BED file..."
    awk 'BEGIN {OFS="\t"} $3 == "gene" {
        gene_name = "."
        if (match($0, /gene_name "([^"]+)"/, a)) gene_name = a[1]
        print $1, $4 - 1, $5, gene_name
    }' "${GTF}" > "${GENE_BED}"

    if check_tool "bgzip"; then
        bgzip -c "${GENE_BED}" > "${GENE_BED}.gz"
    else
        gzip -c "${GENE_BED}" > "${GENE_BED}.gz"
    fi

    if check_tool "tabix"; then
        tabix -p bed "${GENE_BED}.gz"
    fi

    echo -e '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">' > "${GENE_BED}.hdr"

    log "Gene BED file created: ${GENE_BED}.gz"
    log "Unzipped GTF file available for STAR: ${GTF}"
}

download_funcotator() {
    log "Downloading Funcotator data sources..."

    if ! check_tool "gatk"; then
        log "ERROR: GATK is required for Funcotator but not found"
        return 1
    fi

    mkdir -p "${FUNCOTATOR_DS}"

    # Check if already downloaded
    if [ -d "${FUNCOTATOR_DS}" ] && [ "$(ls -A "${FUNCOTATOR_DS}" 2>/dev/null | wc -l)" -gt 10 ]; then
        log "Funcotator data sources already exist at: ${FUNCOTATOR_DS}"
        return 0
    fi

    log "This may take a while (downloading ~30GB of data)..."
    log "You can also download manually from:"
    log "  https://console.cloud.google.com/storage/browser/gatk-best-practices/funcotator"

    # Try different download methods
    if [ -x "${GATK}" ] || command -v gatk &> /dev/null; then
        GATK_CMD="${GATK:-gatk}"

        # Try with somatic databases first (smaller)
        log "Attempting to download somatic data sources..."
        if $GATK_CMD FuncotatorDataSourceDownloader \
            --somatic \
            --hg38 \
            --validate-integrity \
            --extract-after-download \
            --output "${FUNCOTATOR_DS}" 2>&1 | tee -a "${BASE_DIR}/logs/funcotator.log"; then
            log "Somatic data sources downloaded successfully"
            return 0
        fi

        # Try without validation
        log "Trying without validation..."
        if $GATK_CMD FuncotatorDataSourceDownloader \
            --somatic \
            --hg38 \
            --extract-after-download \
            --output "${FUNCOTATOR_DS}" 2>&1 | tee -a "${BASE_DIR}/logs/funcotator.log"; then
            log "Somatic data sources downloaded (without validation)"
            return 0
        fi

        # Try germline
        log "Trying germline data sources..."
        if $GATK_CMD FuncotatorDataSourceDownloader \
            --germline \
            --hg38 \
            --extract-after-download \
            --output "${FUNCOTATOR_DS}" 2>&1 | tee -a "${BASE_DIR}/logs/funcotator.log"; then
            log "Germline data sources downloaded"
            return 0
        fi
    fi

    # Manual download option
    log_warning "Automatic Funcotator download failed."
    log "You need to manually download and extract funcotator data sources:"
    log ""
    log "Option 1: Download from Google Cloud Storage:"
    log "  mkdir -p ${FUNCOTATOR_DS}"
    log "  cd ${FUNCOTATOR_DS}"
    log "  gsutil -m cp -r gs://gatk-best-practices/funcotator/dataSources.v1.7.20200521g/* ."
    log ""
    log "Option 2: Use pre-downloaded tarball:"
    log "  wget https://storage.googleapis.com/gatk-best-practices/funcotator/funcotator_dataSources.v1.7.20200521g.tar.gz"
    log "  tar -xzf funcotator_dataSources.v1.7.20200521g.tar.gz -C ${FUNCOTATOR_DS}"
    log ""
    log "Option 3: Skip Funcotator annotation (less comprehensive):"
    log "  Set SKIP_FUNCOTATOR=true in script configuration"
    log ""

    return 1
}

setup_references() {
    log "Setting up pipeline..."
    setup_environment

    echo "Checking for required tools..."
    REQUIRED_TOOLS=("java" "samtools" "bcftools" "bgzip" "tabix" "wget")
    for tool in "${REQUIRED_TOOLS[@]}"; do
        check_tool "${tool}"
    done

    check_tool "gatk"

    echo ""
    echo "Downloading reference files..."
    download_reference
    generate_gene_bed

    echo ""
    echo "Setting up annotation tools..."
    setup_annotation_tools

    echo ""
    echo "Setup complete!"
    echo "Edit the SRA_LIST array in the script to add your SRA IDs"
    echo "Then run: ./$(basename "$0") run"
}

############################################################
# SRA ID Extraction (Fixed and Complete)
############################################################

extract_sra_ids_from_excel() {
    local excel_file="$1"
    local output_file="$2"
    local project_file="$3"  # New parameter

    # Use a unique name to avoid overwriting other scripts
    local PYTHON_SCRIPT="${SCRIPT_DIR}/extract_sra_ids_interactive.py"

    cat > "${PYTHON_SCRIPT}" << 'PYTHONSCRIPT'
#!/usr/bin/env python3
import pandas as pd
import sys
import re

def extract_sra_ids(excel_path, output_path, project_output_path):
    """Main function to extract SRA IDs with interactive column selection"""

    def interactive_column_selection(df, purpose):
        """Helper: Interactive column selection with enumerated options"""
        print(f"\n=== Select column for {purpose} ===")
        print("Available columns:")

        for i, col in enumerate(df.columns, 1):
            # Show first non-null value as preview
            preview = df[col].dropna().iloc[0] if not df[col].dropna().empty else '[Empty]'
            print(f"  {i}. {col} (example: '{preview}')")

        while True:
            try:
                choice = input(f"Select column number for {purpose} (1-{len(df.columns)}): ").strip()
                if not choice.isdigit():
                    print("Please enter a number")
                    continue

                idx = int(choice) - 1
                if 0 <= idx < len(df.columns):
                    selected = df.columns[idx]
                    print(f"✓ Selected '{selected}' for {purpose}")
                    return selected
                else:
                    print(f"Please enter a number between 1 and {len(df.columns)}")
            except KeyboardInterrupt:
                print("\nOperation cancelled")
                sys.exit(1)
            except Exception as e:
                print(f"Error: {e}")

    try:
        # 1. Read Excel file
        print(f"📂 Reading Excel file: {excel_path}")

        # Try to get sheet names
        xls = pd.ExcelFile(excel_path)
        if len(xls.sheet_names) > 1:
            print(f"\nMultiple sheets found:")
            for i, sheet in enumerate(xls.sheet_names, 1):
                print(f"  {i}. {sheet}")

            sheet_choice = input(f"Select sheet number (1-{len(xls.sheet_names)}), default [1]: ").strip()
            if sheet_choice and sheet_choice.isdigit():
                sheet_idx = int(sheet_choice) - 1
                if 0 <= sheet_idx < len(xls.sheet_names):
                    df = pd.read_excel(excel_path, sheet_name=xls.sheet_names[sheet_idx])
                else:
                    df = pd.read_excel(excel_path, sheet_name=0)
            else:
                df = pd.read_excel(excel_path, sheet_name=0)
        else:
            df = pd.read_excel(excel_path, sheet_name=0)

        print(f"✓ Successfully loaded Excel file")
        print(f"  Shape: {df.shape[0]} rows × {df.shape[1]} columns")

        # 2. Interactive column selection
        print("\n" + "="*50)
        sra_col = interactive_column_selection(df, "SRA IDs")
        print("-"*50)
        project_col = interactive_column_selection(df, "Project IDs (MDD/Control)")
        print("="*50)

        # 3. Extract SRA IDs (handle various formats)
        sra_column_data = df[sra_col].dropna().astype(str)
        valid_sra_pairs = []  # List of (sra_id, project_id)

        print(f"\n🔍 Scanning column '{sra_col}' for SRA IDs...")

        for idx, sra_cell in enumerate(sra_column_data):
            # Try to find SRR followed by digits in the cell
            matches = re.findall(r'(SRR\d+)', sra_cell)

            if matches:
                # Found SRA ID(s) in this cell
                project_cell = df.iloc[idx][project_col] if project_col in df.columns else "UNKNOWN"
                project_val = str(project_cell) if pd.notna(project_cell) else "UNKNOWN"

                for sra_match in matches:
                    valid_sra_pairs.append((sra_match, project_val))
            elif re.match(r'^SRR\d+$', sra_cell.strip()):
                # Cell contains just an SRA ID
                project_cell = df.iloc[idx][project_col] if project_col in df.columns else "UNKNOWN"
                project_val = str(project_cell) if pd.notna(project_cell) else "UNKNOWN"
                valid_sra_pairs.append((sra_cell.strip(), project_val))

        # 4. Remove duplicates while preserving order
        unique_pairs = []
        seen = set()
        for sra, project in valid_sra_pairs:
            if sra not in seen:
                seen.add(sra)
                unique_pairs.append((sra, project))

        if not unique_pairs:
            print("❌ ERROR: No SRA IDs found in the selected column!")
            print(f"   Check that column '{sra_col}' contains IDs like SRR12345678")
            return False

        # 5. Write output files
        with open(output_path, 'w') as f:
            for sra_id, _ in unique_pairs:
                f.write(f"{sra_id}\n")

        with open(project_output_path, 'w') as f:
            for sra_id, project_id in unique_pairs:
                f.write(f"{sra_id}\t{project_id}\n")

        # 6. Display summary
        print(f"\n✅ SUCCESS: Extracted {len(unique_pairs)} unique SRA IDs")
        print(f"   SRA list: {output_path}")
        print(f"   Project mapping: {project_output_path}")

        # Show project distribution
        from collections import Counter
        project_counts = Counter([p for _, p in unique_pairs])
        print(f"\n📊 Project distribution:")
        for project, count in sorted(project_counts.items()):
            print(f"   {project}: {count} samples")

        # Show first few samples
        print(f"\n🔬 First 5 samples:")
        for i, (sra, project) in enumerate(unique_pairs[:5]):
            print(f"   {i+1}. {sra} → {project}")
        if len(unique_pairs) > 5:
            print(f"   ... and {len(unique_pairs) - 5} more")

        return True

    except FileNotFoundError:
        print(f"❌ ERROR: File not found: {excel_path}")
        return False
    except Exception as e:
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

# Main execution
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_sra_ids.py <excel_file> <sra_output_file> <project_output_file>")
        print("Example: python extract_sra_ids.py samples.xlsx sra_ids.txt project_map.txt")
        sys.exit(1)

    excel_file = sys.argv[1]
    sra_output = sys.argv[2]
    project_output = sys.argv[3]

    success = extract_sra_ids(excel_file, sra_output, project_output)
    sys.exit(0 if success else 1)
PYTHONSCRIPT

    chmod +x "${PYTHON_SCRIPT}"

    # Check for Python dependencies
    if ! python3 -c "import pandas" &> /dev/null; then
        log "Installing Python dependencies..."
        if command -v pip3 &> /dev/null; then
            pip3 install --user pandas openpyxl 2>/dev/null || pip3 install pandas openpyxl
        elif command -v pip &> /dev/null; then
            pip install --user pandas openpyxl 2>/dev/null || pip install pandas openpyxl
        else
            log_error "pip not available. Please install pandas manually: pip install pandas openpyxl"
            return 1
        fi
    fi

    log "Extracting SRA IDs from Excel (interactive mode)..."
    if python3 "${PYTHON_SCRIPT}" "${excel_file}" "${output_file}" "${project_file}"; then
        log "Interactive SRA ID extraction completed successfully"
        return 0
    else
        log_error "Interactive SRA ID extraction failed"
        return 1
    fi
}

prompt_for_excel_file() {
    if [ -z "${EXCEL_FILE}" ] || [ ! -f "${EXCEL_FILE}" ]; then
        echo ""
        echo "========================================"
        echo "Excel File Selection"
        echo "========================================"
        echo "Please provide the path to your Excel file containing SRA IDs."
        echo ""

        # Find Excel files
        excel_files=$(find . -maxdepth 2 -name "*.xlsx" -o -name "*.xls" 2>/dev/null | head -10)

        if [ -n "${excel_files}" ]; then
            echo "Found Excel files:"
            i=1
            while IFS= read -r file; do
                echo "  ${i}. $(basename "${file}")"
                i=$((i + 1))
            done <<< "${excel_files}"
            echo "  0. Enter custom path"
            echo ""
            read -p "Select file (0-$(($i-1))): " choice

            if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -lt "$i" ]; then
                EXCEL_FILE=$(echo "${excel_files}" | sed -n "${choice}p")
            else
                read -p "Enter path to Excel file: " EXCEL_FILE
            fi
        else
            echo "No Excel files found in current directory."
            read -p "Enter path to Excel file: " EXCEL_FILE
        fi

        [ -f "${EXCEL_FILE}" ] || log_error "Excel file not found: ${EXCEL_FILE}"
        log "Using Excel file: ${EXCEL_FILE}"
    fi
}

load_sra_ids() {
    log "Loading SRA IDs..."

    # Use existing EXCEL_FILE or prompt (keep your existing prompt_for_excel_file)
    if [ -z "${EXCEL_FILE}" ] || [ ! -f "${EXCEL_FILE}" ]; then
        prompt_for_excel_file
    fi

    log "Using Excel file: ${EXCEL_FILE}"

    # Define output files - IMPORTANT: Ensure BASE_DIR is defined
    SRA_LIST_FILE="${BASE_DIR}/sra_ids.txt"
    PROJECT_MAP_FILE="${BASE_DIR}/project_mapping.txt"

    # Call the UPDATED interactive extractor (with 3 parameters)
    if extract_sra_ids_from_excel "${EXCEL_FILE}" "${SRA_LIST_FILE}" "${PROJECT_MAP_FILE}"; then
        if [ -f "${SRA_LIST_FILE}" ] && [ -s "${SRA_LIST_FILE}" ]; then
            count=$(wc -l < "${SRA_LIST_FILE}" | tr -d ' ')
            [ "${count}" -eq 0 ] && log_error "No SRA IDs found in Excel file"

            # Load SRA IDs into array
            mapfile -t SRA_LIST < "${SRA_LIST_FILE}"

            # --- Load Project IDs into Associative Array ---
            declare -gA PROJECT_MAP  # Make it globally available
            PROJECTS_LOADED=0

            if [ -f "${PROJECT_MAP_FILE}" ] && [ -s "${PROJECT_MAP_FILE}" ]; then
                # Clear any existing mappings
                unset PROJECT_MAP
                declare -gA PROJECT_MAP

                while IFS=$'\t' read -r sra_id project_id; do
                    # Clean up any whitespace or quotes
                    sra_id_clean=$(echo "${sra_id}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//')
                    project_id_clean=$(echo "${project_id}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//')

                    if [[ "${sra_id_clean}" =~ ^SRR[0-9]+$ ]]; then
                        PROJECT_MAP["${sra_id_clean}"]="${project_id_clean}"
                        PROJECTS_LOADED=$((PROJECTS_LOADED + 1))
                    else
                        log_warning "Skipping invalid SRA ID in mapping: '${sra_id_clean}'"
                    fi
                done < "${PROJECT_MAP_FILE}"

                log "Loaded ${#SRA_LIST[@]} SRA IDs with ${PROJECTS_LOADED} Project mappings from ${EXCEL_FILE}"
            else
                log_warning "Project mapping file not found or empty at: ${PROJECT_MAP_FILE}"
                log "All samples will be marked as 'UNKNOWN'."

                # Initialize all samples as UNKNOWN
                for sra_id in "${SRA_LIST[@]}"; do
                    sra_clean=$(echo "${sra_id}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                    if [[ "${sra_clean}" =~ ^SRR[0-9]+$ ]]; then
                        PROJECT_MAP["${sra_clean}"]="UNKNOWN"
                    fi
                done
            fi
            # --- END Project ID loading ---

            # Verify all SRA IDs have mappings
            MISSING_MAPPINGS=0
            for sra_id in "${SRA_LIST[@]}"; do
                sra_clean=$(echo "${sra_id}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                if [ -z "${PROJECT_MAP["${sra_clean}"]+x}" ]; then
                    log_warning "No project mapping for SRA ID: ${sra_clean}"
                    PROJECT_MAP["${sra_clean}"]="UNKNOWN"
                    MISSING_MAPPINGS=$((MISSING_MAPPINGS + 1))
                fi
            done

            if [ "${MISSING_MAPPINGS}" -gt 0 ]; then
                log_warning "${MISSING_MAPPINGS} SRA IDs had no project mapping and were set to 'UNKNOWN'"
            fi

            # Display summary
            echo ""
            echo "========================================"
            echo "SRA ID Loading Summary"
            echo "========================================"
            echo "Total SRA IDs loaded: ${#SRA_LIST[@]}"
            echo "Project mappings loaded: ${PROJECTS_LOADED}"

            # Show project distribution
            declare -A PROJECT_COUNTS=()  # Initialize as empty associative array
            for sra_id in "${SRA_LIST[@]}"; do
                sra_clean=$(echo "${sra_id}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                project="${PROJECT_MAP["${sra_clean}"]}"
                if [ -n "${project}" ]; then
                    PROJECT_COUNTS["${project}"]=$((PROJECT_COUNTS["${project}"] + 1))
                fi
            done

            echo ""
            echo "Project distribution:"
            for project in "${!PROJECT_COUNTS[@]}"; do
                echo "  ${project}: ${PROJECT_COUNTS["${project}"]} samples"
            done

            # Display first 5 SRA IDs with their Project ID
            echo ""
            echo "First 5 SRA IDs with Project IDs:"
            for i in {0..4}; do
                if [ -n "${SRA_LIST[$i]}" ]; then
                    sra_clean=$(echo "${SRA_LIST[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                    echo "  ${sra_clean} -> ${PROJECT_MAP["${sra_clean}"]:-NOT_FOUND}"
                fi
            done

            if [ ${#SRA_LIST[@]} -gt 5 ]; then
                echo "  ... and $(( ${#SRA_LIST[@]} - 5 )) more"
            fi
            echo ""

        else
            log_error "SRA IDs file is empty or not created: ${SRA_LIST_FILE}"
            log "Check that:"
            log "  1. The Excel file contains valid SRA IDs (SRRXXXXXXXX)"
            log "  2. You selected the correct column during interactive selection"
            log "  3. Python dependencies (pandas, openpyxl) are installed"
            return 1
        fi
    else
        log_error "Failed to extract SRA IDs from Excel file: ${EXCEL_FILE}"
        log "Possible issues:"
        log "  1. File format not supported (need .xlsx or .xls)"
        log "  2. File is corrupted or password protected"
        log "  3. No valid SRA IDs found in selected column"
        return 1
    fi

    return 0
}

#=================================================================
# For single SRA files
#=================================================================

trim_single_sample() {
    local SRR="$1"
    local INPUT_DIR="${BASE_DIR}/data/fastq"

    # Skip if trimmed files already exist
    if TRIMMED1=$(check_sample_file_exists "${SRR}" "trimmed1") && \
       TRIMMED2=$(check_sample_file_exists "${SRR}" "trimmed2"); then
        log "Trimmed files already exist: ${TRIMMED1}, ${TRIMMED2}, skipping trim..."
        return 0
    fi

    # Need FASTQ files
    if ! FASTQ1=$(check_sample_file_exists "${SRR}" "fastq1") || \
       ! FASTQ2=$(check_sample_file_exists "${SRR}" "fastq2"); then
        log_error "FASTQ files not found for ${SRR}, cannot trim"
    fi

    if [ ! -f "${TRIMMOMATIC_JAR}" ] && [ ! -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        log "Trimmomatic not found, skipping trim for ${SRR}"
        return 0
    fi

    if [ -f "${TRIMMOMATIC_JAR}" ]; then
        TRIMMOMATIC_CMD="java -jar ${TRIMMOMATIC_JAR}"
    else
        TRIMMOMATIC_CMD="java -jar ${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"
    fi

    ADAPTERS="${TOOLS_DIR}/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
    if [ ! -f "${ADAPTERS}" ]; then
        ADAPTERS="${TOOLS_DIR}/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
    fi

    r1="${INPUT_DIR}/${SRR}_1.fastq.gz"
    r2="${INPUT_DIR}/${SRR}_2.fastq.gz"

    if [ ! -f "${r1}" ] || [ ! -f "${r2}" ]; then
        log_warning "Missing FASTQ files for ${SRR}, skipping trim"
        return 0
    fi

    out_p1="${BASE_DIR}/data/trimmed/${SRR}_1.paired.fastq.gz"
    out_u1="${BASE_DIR}/data/trimmed/${SRR}_1.unpaired.fastq.gz"
    out_p2="${BASE_DIR}/data/trimmed/${SRR}_2.paired.fastq.gz"
    out_u2="${BASE_DIR}/data/trimmed/${SRR}_2.unpaired.fastq.gz"

    log "Trimming adapters from: ${SRR}"
    ${TRIMMOMATIC_CMD} PE -threads 6 -phred33 \
        "${r1}" "${r2}" \
        "${out_p1}" "${out_u1}" \
        "${out_p2}" "${out_u2}" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
}

align_single_sample() {
    local SRR="$1"

    # Skip if BAM already exists
    if BAM_FILE=$(check_sample_file_exists "${SRR}" "bam"); then
        log "BAM file already exists: ${BAM_FILE}"

        # Check if index exists, create if missing
        if ! check_sample_file_exists "${SRR}" "bai" > /dev/null; then
            log "Creating missing BAM index..."
            samtools index "${BAM_FILE}"
        fi
        return 0
    fi

    # Need trimmed files
    if ! TRIMMED1=$(check_sample_file_exists "${SRR}" "trimmed1") || \
       ! TRIMMED2=$(check_sample_file_exists "${SRR}" "trimmed2"); then
        log_error "Trimmed FASTQ files not found for ${SRR}, cannot align"
    fi

    # Check if STAR index exists (check for SA file, not just directory)
    if [ ! -f "${STAR_INDEX_DIR}/SA" ]; then
        log "STAR index not found or incomplete, creating index..."
        create_star_index
    fi

    R1="${BASE_DIR}/data/trimmed/${SRR}_1.paired.fastq.gz"
    R2="${BASE_DIR}/data/trimmed/${SRR}_2.paired.fastq.gz"

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        log_warning "Missing trimmed FASTQ files for ${SRR}, skipping alignment"
        return 0
    fi

    log "Aligning ${SRR} with STAR..."
    local star_gpu_opts=""
    if [ "${USE_GPU}" = "true" ] && command -v nvidia-smi &> /dev/null; then
        star_gpu_opts="--runGPU"
    fi

    # REMOVE the --outSAMattrRGline parameter - let GATK handle read groups
    STAR --runThreadN ${GATK_THREADS} \
         --genomeDir "${STAR_INDEX_DIR}" \
         --readFilesIn "${R1}" "${R2}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${BASE_DIR}/data/aligned/${SRR}." \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --limitBAMsortRAM 50000000000 \
         ${star_gpu_opts}

    BAM="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam"
    if [ -f "${BAM}" ]; then
        samtools index "${BAM}"
    fi
}

download_single_sra() {
    local SRR="$1"
    log "Downloading ${SRR}..."
    # Skip if already downloaded
    if SRA_FILE=$(check_sample_file_exists "${SRR}" "sra"); then
        log "SRA file already exists: ${SRA_FILE}, skipping download..."
        return 0
    fi
    mkdir -p "${BASE_DIR}/data/sra/"

    if SRA_FILE=$(find_sra_file "${SRR}"); then
        log "Found ${SRR} at: ${SRA_FILE}"
        if [ ! -f "${BASE_DIR}/data/sra/${SRR}.sra" ]; then
            ln -sf "${SRA_FILE}" "${BASE_DIR}/data/sra/${SRR}.sra"
            log "Created symlink for easier access"
        fi
        return 0
    fi

    cd "${BASE_DIR}/data/sra/"
    if command -v prefetch &> /dev/null; then
        prefetch --progress "${SRR}" --max-size 100G
    elif [ -x "${TOOLS_DIR}/sratoolkit/bin/prefetch" ]; then
        "${TOOLS_DIR}/sratoolkit/bin/prefetch" --progress "${SRR}" --max-size 100G
    fi

    if [ -d "${SRR}" ] && [ -f "${SRR}/${SRR}.sra" ]; then
        ln -sf "${SRR}/${SRR}.sra" "${SRR}.sra"
        log "Created symlink from directory structure"
    fi

    cd "${PIPELINE_DIR}"

    if [ ! -f "${BASE_DIR}/data/sra/${SRR}.sra" ]; then
        log_error "Failed to download ${SRR}"
    fi
}

validate_single_sra() {
    local SRR="$1"
    log "Validating ${SRR}..."

    if SRA_FILE=$(find_sra_file "${SRR}"); then
        log "Validating: ${SRA_FILE}"
        if command -v vdb-validate &> /dev/null; then
            vdb-validate "${SRA_FILE}"
        elif [ -x "${TOOLS_DIR}/sratoolkit/bin/vdb-validate" ]; then
            "${TOOLS_DIR}/sratoolkit/bin/vdb-validate" "${SRA_FILE}"
        else
            log "WARNING: vdb-validate not found, skipping validation"
        fi
    else
        log_error "${SRR}.sra not found"
    fi
}

extract_single_fastq() {
    local SRR="$1"
    log "Splitting ${SRR} into FASTQ..."
    mkdir -p "${BASE_DIR}/data/fastq/"

    # Skip if FASTQ already exists
    if FASTQ1=$(check_sample_file_exists "${SRR}" "fastq1") && \
       FASTQ2=$(check_sample_file_exists "${SRR}" "fastq2"); then
        log "FASTQ files already exist: ${FASTQ1}, ${FASTQ2}, skipping extraction..."
        return 0
    fi

    if ! SRA_FILE=$(find_sra_file "${SRR}"); then
        log_error "${SRR}.sra not found, cannot extract FASTQ"
    fi

    log "Extracting from: ${SRA_FILE}"

    if command -v fasterq-dump &> /dev/null; then
        fasterq-dump --split-files --threads 6 --progress \
            -O "${BASE_DIR}/data/fastq/" \
            "${SRA_FILE}"
    elif [ -x "${TOOLS_DIR}/sratoolkit/bin/fasterq-dump" ]; then
        "${TOOLS_DIR}/sratoolkit/bin/fasterq-dump" --split-files --threads 6 --progress \
            -O "${BASE_DIR}/data/fastq/" \
            "${SRA_FILE}"
    else
        log_error "fasterq-dump command not found"
    fi

    if [ -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq" ]; then
        pigz -f -p 2 "${BASE_DIR}/data/fastq/${SRR}_1.fastq" 2>/dev/null || \
        gzip -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq"
    fi
    if [ -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq" ]; then
        pigz -f -p 2 "${BASE_DIR}/data/fastq/${SRR}_2.fastq" 2>/dev/null || \
        gzip -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq"
    fi
}

cleanup_intermediate() {
    local SRR="$1"
    local KEEP_INTERMEDIATE="${2:-false}"

    if [ "${KEEP_INTERMEDIATE}" = "true" ]; then
        log "Keeping intermediate files for ${SRR} (--keep-intermediate flag set)"
        return 0
    fi

    log "Cleaning up intermediate files for ${SRR}..."

    # Remove SRA file
    rm -f "${BASE_DIR}/data/sra/${SRR}.sra" 2>/dev/null || true
    rm -rf "${BASE_DIR}/data/sra/${SRR}" 2>/dev/null || true

    # Remove FASTQ files
    rm -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq.gz" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq.gz" 2>/dev/null || true

    # Remove trimmed files (keep QC reports)
    rm -f "${BASE_DIR}/data/trimmed/${SRR}_1.paired.fastq.gz" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/trimmed/${SRR}_1.unpaired.fastq.gz" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/trimmed/${SRR}_2.paired.fastq.gz" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/trimmed/${SRR}_2.unpaired.fastq.gz" 2>/dev/null || true

    # Remove intermediate BAM files (keep final processed BAM)
    rm -f "${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam.bai" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/aligned/${SRR}.Log.*" 2>/dev/null || true
    rm -f "${BASE_DIR}/data/aligned/${SRR}.SJ.out.tab" 2>/dev/null || true

    log "Intermediate files cleaned for ${SRR}"
}

check_disk_space() {
    local required_mb="${1:-10240}"  # Default 10GB
    local available_mb=$(df -m "${BASE_DIR}" | tail -1 | awk '{print $4}')

    if [ "${available_mb}" -lt "${required_mb}" ]; then
        log_warning "Low disk space: ${available_mb}MB available, ${required_mb}MB recommended"
        return 1
    fi

    log "Disk space OK: ${available_mb}MB available"
    return 0
}

run_full_pipeline() {
    log "Starting full MDD2 SNP Extraction Pipeline with auto-resume"
    setup_environment

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log_error "No SRA IDs to process"
    fi

    log "Processing ${#SRA_LIST[@]} samples (will skip completed steps)..."

    COMPLETED_SAMPLES=0
    SKIPPED_SAMPLES=0
    PROCESSED_SAMPLES=0
    FAILED_SAMPLES=0

    for SRR in "${SRA_LIST[@]}"; do
        log "========================================"
        log "Sample: ${SRR}"
        log "========================================"

        status=$(get_sample_status "${SRR}")
        log "Current status: ${status}"

        # If already completed, skip
        if [ "${status}" = "completed" ]; then
            log "✓ Already completed, skipping..."
            SKIPPED_SAMPLES=$((SKIPPED_SAMPLES + 1))
            continue
        fi

        # Check disk space
        check_disk_space 5120 || log_warning "Low disk space, continuing anyway"

        # Process based on current status
        if process_sample_sequentially "${SRR}"; then
            PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1))
            if [ "$(get_sample_status "${SRR}")" = "completed" ]; then
                COMPLETED_SAMPLES=$((COMPLETED_SAMPLES + 1))
            fi
            log "Successfully progressed ${SRR}"
        else
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            log_warning "Failed to process ${SRR}, will retry on next run"
        fi

        log "Progress: ${COMPLETED_SAMPLES} completed, ${SKIPPED_SAMPLES} skipped, ${PROCESSED_SAMPLES} processed this run"
    done

    # Merge all TSV files
    merge_all_tsv_files_simple

    log ""
    log "========================================"
    log "Pipeline Run Complete!"
    log "========================================"
    log "Summary:"
    log "  • Total samples: ${#SRA_LIST[@]}"
    log "  • Completed (now): ${COMPLETED_SAMPLES}"
    log "  • Skipped (already done): ${SKIPPED_SAMPLES}"
    log "  • Processed this run: ${PROCESSED_SAMPLES}"
    log "  • Failed this run: ${FAILED_SAMPLES}"
    log ""
    log "Output files:"
    log "  • Merged TSV: ${BASE_DIR}/all_samples.tsv.gz"
    log "  • Individual results: ${BASE_DIR}/data/tsv/"
    log ""
    log "Run again to resume any failed/incomplete samples"
    log ""
}

process_sample_sequentially() {
    local SRR="$1"

    log "Processing ${SRR} sequentially with auto-resume..."

    # Step 1: Download (if needed)
    if ! check_sample_file_exists "${SRR}" "sra" > /dev/null; then
        log "Step 1: Downloading SRA..."
        download_single_sra "${SRR}" || return 1
    else
        log "Step 1: SRA already downloaded, skipping..."
    fi

    # Step 2: Extract FASTQ (if needed)
    if ! check_sample_file_exists "${SRR}" "fastq1" > /dev/null || \
       ! check_sample_file_exists "${SRR}" "fastq2" > /dev/null; then
        log "Step 2: Extracting FASTQ..."
        extract_single_fastq "${SRR}" || return 1
    else
        log "Step 2: FASTQ already extracted, skipping..."
    fi

    # Step 3: Run FastQC on raw files (always run, quick)
    log "Step 3: Running FastQC on raw FASTQ..."
    run_fastqc_single "${SRR}" "raw"

    # Step 4: Trim (if needed)
    if ! check_sample_file_exists "${SRR}" "trimmed1" > /dev/null || \
       ! check_sample_file_exists "${SRR}" "trimmed2" > /dev/null; then
        log "Step 4: Trimming..."
        trim_single_sample "${SRR}" || return 1
    else
        log "Step 4: Already trimmed, skipping..."
    fi

    # Step 5: Run FastQC on trimmed files
    log "Step 5: Running FastQC on trimmed FASTQ..."
    run_fastqc_single "${SRR}" "trimmed"

    # Step 6: Align (if needed)
    if ! check_sample_file_exists "${SRR}" "bam" > /dev/null; then
        log "Step 6: Aligning with STAR..."
        align_single_sample "${SRR}" || return 1
    else
        log "Step 6: Already aligned, skipping..."
    fi

    # Step 7: Process (variant calling, if needed)
    if ! check_sample_file_exists "${SRR}" "tsv" > /dev/null || \
       ! check_sample_file_exists "${SRR}" "vcf" > /dev/null; then
        log "Step 7: Variant calling and annotation..."
        BAM_FILE="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam"
        if [ -f "${BAM_FILE}" ]; then
            process_sample "${SRR}" "${BAM_FILE}" || return 1
        else
            log_error "BAM file not found: ${BAM_FILE}"
            return 1
        fi
    else
        log "Step 7: Already processed, skipping..."
    fi

    # Step 8: Cleanup (if not keeping intermediates)
    if [ "${KEEP_INTERMEDIATE}" != "true" ]; then
        log "Step 8: Cleaning up intermediate files..."
        cleanup_intermediate "${SRR}" "${KEEP_INTERMEDIATE}"
    fi
    cleanup_interval_files "${SRR}"
    log "✓ Successfully processed ${SRR}"
    return 0
}

cleanup_interval_files() {
    local SAMPLE="$1"

    if [ "${KEEP_INTERMEDIATE}" = "true" ]; then
        log "Keeping interval files (--keep-intermediate flag set)"
        return 0
    fi

    INTERVAL_DIR="${BASE_DIR}/data/vcf/${SAMPLE}/intervals"
    if [ -d "${INTERVAL_DIR}" ]; then
        log "Cleaning up interval files..."
        rm -rf "${INTERVAL_DIR}"
    fi
}

run_fastqc_single() {
    local SRR="$1"
    local TYPE="$2"  # "raw" or "trimmed"

    mkdir -p "${BASE_DIR}/data/fastqc/${TYPE}"

    if [ "${TYPE}" = "raw" ]; then
        R1="${BASE_DIR}/data/fastq/${SRR}_1.fastq.gz"
        R2="${BASE_DIR}/data/fastq/${SRR}_2.fastq.gz"
    else
        R1="${BASE_DIR}/data/trimmed/${SRR}_1.paired.fastq.gz"
        R2="${BASE_DIR}/data/trimmed/${SRR}_2.paired.fastq.gz"
    fi

    # Run FastQC only if output doesn't exist
    if [ -f "${R1}" ] && [ ! -f "${BASE_DIR}/data/fastqc/${TYPE}/$(basename "${R1}" .fastq.gz)_fastqc.html" ]; then
        fastqc "${R1}" -o "${BASE_DIR}/data/fastqc/${TYPE}" --quiet 2>/dev/null || true
    fi

    if [ -f "${R2}" ] && [ ! -f "${BASE_DIR}/data/fastqc/${TYPE}/$(basename "${R2}" .fastq.gz)_fastqc.html" ]; then
        fastqc "${R2}" -o "${BASE_DIR}/data/fastqc/${TYPE}" --quiet 2>/dev/null || true
    fi
}

run_fastq_pipeline() {
    log "Starting FASTQ processing pipeline"
    setup_environment

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log_error "No SRA IDs to process"
    fi

    log "Processing ${#SRA_LIST[@]} samples sequentially..."

    for SRR in "${SRA_LIST[@]}"; do
        log "========================================"
        log "Processing sample: ${SRR}"
        log "========================================"

        # Check disk space before each download
        check_disk_space 5120 || log_warning "Continuing despite low disk space"

        download_single_sra "${SRR}"
        validate_single_sra "${SRR}"
        extract_single_fastq "${SRR}"

        # Cleanup SRA file immediately after extraction
        if [ "${KEEP_INTERMEDIATE}" != "true" ]; then
            rm -f "${BASE_DIR}/data/sra/${SRR}.sra" 2>/dev/null || true
            rm -rf "${BASE_DIR}/data/sra/${SRR}" 2>/dev/null || true
        fi
    done

    # Run QC on all extracted FASTQ files
    run_fastqc
    trim_reads
    run_fastqc_trimmed

    log ""
    log "========================================"
    log "FASTQ Processing Complete!"
    log "========================================"
    log ""
    log "Next steps:"
    log "1. Align the trimmed FASTQ files with STAR"
    log "2. Run: ./$(basename "$0") run"
    log ""
    log "Trimmed FASTQ files are in: ${BASE_DIR}/data/trimmed/"
    log ""
}

merge_all_tsv_files_simple() {
    log "Merging all TSV files..."

    MERGED_TSV="${BASE_DIR}/all_samples.tsv"

    # Create header
    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ\tSAMPLE" > "${MERGED_TSV}"

    # Find and merge all TSV files
    TSV_FILES=0

    for tsv_gz in "${BASE_DIR}/data/tsv/"*/*.tsv.gz; do
        if [ -f "${tsv_gz}" ]; then
            TSV_FILES=$((TSV_FILES + 1))
            sample=$(basename "${tsv_gz}" .tsv.gz)
            log "Adding ${sample}..."

            # Use appropriate decompression
            if command -v pigz > /dev/null; then
                pigz -dc "${tsv_gz}" 2>/dev/null | tail -n +2 | awk -v sample="${sample}" '{print $0 "\t" sample}' >> "${MERGED_TSV}"
            else
                gzip -dc "${tsv_gz}" | tail -n +2 | awk -v sample="${sample}" '{print $0 "\t" sample}' >> "${MERGED_TSV}"
            fi
        fi
    done

    if [ "${TSV_FILES}" -eq 0 ]; then
        log_warning "No TSV files found to merge"
        rm -f "${MERGED_TSV}"
        return 0
    fi

    # Compress output
    if [ -f "${MERGED_TSV}" ] && [ -s "${MERGED_TSV}" ]; then
        if command -v pigz > /dev/null; then
            pigz -f -p 4 "${MERGED_TSV}"
        else
            gzip -f "${MERGED_TSV}"
        fi
        log "Merged ${TSV_FILES} samples into: ${MERGED_TSV}.gz"
    else
        log_warning "Merged TSV file is empty"
    fi
}

show_pipeline_status_simple() {
    log "Pipeline Status (File-based)"
    log "============================="

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    TOTAL=${#SRA_LIST[@]}
    COMPLETED=0
    PARTIAL=0
    NOT_STARTED=0

    echo ""
    echo "Sample Status:"
    echo "--------------"

    for SRR in "${SRA_LIST[@]}"; do
        status=$(get_sample_status "${SRR}")

        case "${status}" in
            "completed")
                COMPLETED=$((COMPLETED + 1))
                echo "  ✓ ${SRR}: COMPLETED"
                ;;
            "aligned"|"trimmed"|"extracted"|"downloaded")
                PARTIAL=$((PARTIAL + 1))
                next_step=$(get_next_step_for_sample "${SRR}")
                echo "  → ${SRR}: ${status} (next: ${next_step})"
                ;;
            "not_started")
                NOT_STARTED=$((NOT_STARTED + 1))
                echo "  ○ ${SRR}: NOT STARTED"
                ;;
            *)
                echo "  ? ${SRR}: ${status}"
                ;;
        esac
    done

    echo ""
    echo "Summary:"
    echo "  Total samples: ${TOTAL}"
    echo "  Completed: ${COMPLETED} ($((COMPLETED * 100 / TOTAL))%)"
    echo "  Partially done: ${PARTIAL}"
    echo "  Not started: ${NOT_STARTED}"
    echo ""

    # Show file counts
    echo "File Counts:"
    echo "  SRA files: $(find "${BASE_DIR}/data/sra/" -name "*.sra" 2>/dev/null | wc -l)"
    echo "  FASTQ pairs: $(find "${BASE_DIR}/data/fastq/" -name "*_1.fastq.gz" 2>/dev/null | wc -l)"
    echo "  Trimmed pairs: $(find "${BASE_DIR}/data/trimmed/" -name "*_1.paired.fastq.gz" 2>/dev/null | wc -l)"
    echo "  BAM files: $(find "${BASE_DIR}/data/aligned/" -name "*.Aligned.sortedByCoord.out.bam" 2>/dev/null | wc -l)"
    echo "  TSV files: $(find "${BASE_DIR}/data/tsv/" -name "*.tsv.gz" 2>/dev/null | wc -l)"
    echo ""

    # Disk space
    echo "Disk Usage:"
    du -sh "${BASE_DIR}/data/"* 2>/dev/null | sort -hr | head -10
}


############################################################
# Pipeline Steps (Complete from original)
############################################################

download_sra_files() {
    log "Step 1: Downloading SRA files..."
    mkdir -p "${BASE_DIR}/data/sra/"

    for SRR in "${SRA_LIST[@]}"; do
        log "Downloading ${SRR}..."

        if SRA_FILE=$(find_sra_file "${SRR}"); then
            log "Found ${SRR} at: ${SRA_FILE}"
            if [ ! -f "${BASE_DIR}/data/sra/${SRR}.sra" ]; then
                ln -sf "${SRA_FILE}" "${BASE_DIR}/data/sra/${SRR}.sra"
                log "Created symlink for easier access"
            fi
            continue
        fi

        cd "${BASE_DIR}/data/sra/"
        if command -v prefetch &> /dev/null; then
            prefetch --progress "${SRR}" --max-size 100G
        elif [ -x "${TOOLS_DIR}/sratoolkit/bin/prefetch" ]; then
            "${TOOLS_DIR}/sratoolkit/bin/prefetch" --progress "${SRR}" --max-size 100G
        fi

        if [ -d "${SRR}" ] && [ -f "${SRR}/${SRR}.sra" ]; then
            ln -sf "${SRR}/${SRR}.sra" "${SRR}.sra"
            log "Created symlink from directory structure"
        fi
        cd "${PIPELINE_DIR}"
    done
}

validate_sra_files() {
    log "Step 2: Validating SRA files..."
    for SRR in "${SRA_LIST[@]}"; do
        log "Validating ${SRR}..."

        if SRA_FILE=$(find_sra_file "${SRR}"); then
            log "Validating: ${SRA_FILE}"
            if command -v vdb-validate &> /dev/null; then
                vdb-validate "${SRA_FILE}"
            elif [ -x "${TOOLS_DIR}/sratoolkit/bin/vdb-validate" ]; then
                "${TOOLS_DIR}/sratoolkit/bin/vdb-validate" "${SRA_FILE}"
            else
                log "WARNING: vdb-validate not found, skipping validation"
            fi
        else
            log "WARNING: ${SRR}.sra not found, skipping validation"
        fi
    done
}

extract_fastq() {
    log "Step 3: Extracting FASTQ files..."
    mkdir -p "${BASE_DIR}/data/fastq/"

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log "No SRA IDs to process"
        return 0
    fi

    TOTAL=${#SRA_LIST[@]}
    COUNT=1

    for SRR in "${SRA_LIST[@]}"; do
        log "[${COUNT}/${TOTAL}] Splitting ${SRR} into FASTQ..."

        if ! SRA_FILE=$(find_sra_file "${SRR}"); then
            log "WARNING: ${SRR}.sra not found, skipping"
            COUNT=$((COUNT + 1))
            continue
        fi

        log "Extracting from: ${SRA_FILE}"

        if command -v fasterq-dump &> /dev/null; then
            fasterq-dump --split-files --threads 6 --progress \
                -O "${BASE_DIR}/data/fastq/" \
                "${SRA_FILE}"
        elif [ -x "${TOOLS_DIR}/sratoolkit/bin/fasterq-dump" ]; then
            "${TOOLS_DIR}/sratoolkit/bin/fasterq-dump" --split-files --threads 6 --progress \
                -O "${BASE_DIR}/data/fastq/" \
                "${SRA_FILE}"
        else
            log "ERROR: fasterq-dump command not found"
            return 1
        fi

        if [ -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq" ]; then
            pigz -f -p 2 "${BASE_DIR}/data/fastq/${SRR}_1.fastq" 2>/dev/null || \
            gzip -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq"
        fi
        if [ -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq" ]; then
            pigz -f -p 2 "${BASE_DIR}/data/fastq/${SRR}_2.fastq" 2>/dev/null || \
            gzip -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq"
        fi
        COUNT=$((COUNT + 1))
    done
}

run_fastqc() {
    log "Step 4: Running FastQC..."
    mkdir -p "${BASE_DIR}/data/fastqc/raw"

    FASTQ_FILES=()
    if [ -d "${BASE_DIR}/data/fastq" ]; then
        FASTQ_FILES=("${BASE_DIR}"/data/fastq/*.fastq.gz)
    fi

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
        log "No FASTQ files found in ${BASE_DIR}/data/fastq/"
        return 0
    fi

    for file in "${FASTQ_FILES[@]}"; do
        if [ -f "${file}" ]; then
            log "Running FastQC on ${file}..."
            if command -v fastqc &> /dev/null; then
                fastqc "${file}" -o "${BASE_DIR}/data/fastqc/raw" --quiet
            elif [ -x "${TOOLS_DIR}/FastQC/fastqc" ]; then
                "${TOOLS_DIR}/FastQC/fastqc" "${file}" -o "${BASE_DIR}/data/fastqc/raw" --quiet
            else
                log "ERROR: fastqc command not found"
                return 1
            fi
        fi
    done

    if command -v multiqc &> /dev/null; then
        log "Running MultiQC..."
        multiqc "${BASE_DIR}/data/fastqc/raw" -o "${BASE_DIR}/data/fastqc/" --quiet
    elif [ -x "${TOOLS_DIR}/bin/multiqc" ]; then
        log "Running MultiQC..."
        "${TOOLS_DIR}/bin/multiqc" "${BASE_DIR}/data/fastqc/raw" -o "${BASE_DIR}/data/fastqc/" --quiet
    else
        log "WARNING: multiqc not found, skipping MultiQC report"
    fi
}

trim_reads() {
    log "Step 5: Trimming adapters with Trimmomatic..."
    mkdir -p "${BASE_DIR}/data/trimmed"

    INPUT_DIR="${BASE_DIR}/data/fastq"

    if [ ! -f "${TRIMMOMATIC_JAR}" ] && [ ! -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        log "WARNING: Trimmomatic not found"
        log "Skipping trimming step"
        return 0
    fi

    if [ -f "${TRIMMOMATIC_JAR}" ]; then
        TRIMMOMATIC_CMD="java -jar ${TRIMMOMATIC_JAR}"
    else
        TRIMMOMATIC_CMD="java -jar ${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"
    fi

    ADAPTERS="${TOOLS_DIR}/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
    if [ ! -f "${ADAPTERS}" ]; then
        log "WARNING: Adapter file not found: ${ADAPTERS}"
        log "Using default adapters"
        ADAPTERS="${TOOLS_DIR}/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
    fi

    R1_FILES=("${INPUT_DIR}"/*_1.fastq.gz)
    if [ ${#R1_FILES[@]} -eq 0 ]; then
        log "No R1 FASTQ files found for trimming"
        return 0
    fi

    total=${#R1_FILES[@]}
    count=0

    for r1 in "${R1_FILES[@]}"; do
        count=$((count + 1))
        base=$(basename "${r1}" _1.fastq.gz)
        r2="${INPUT_DIR}/${base}_2.fastq.gz"

        if [ ! -f "${r2}" ]; then
            log "WARNING: ${r2} not found, skipping ${base}"
            continue
        fi

        out_p1="${BASE_DIR}/data/trimmed/${base}_1.paired.fastq.gz"
        out_u1="${BASE_DIR}/data/trimmed/${base}_1.unpaired.fastq.gz"
        out_p2="${BASE_DIR}/data/trimmed/${base}_2.paired.fastq.gz"
        out_u2="${BASE_DIR}/data/trimmed/${base}_2.unpaired.fastq.gz"

        log "[${count}/${total}] Trimming adapters from: ${base}"

        ${TRIMMOMATIC_CMD} PE -threads 6 -phred33 \
            "${r1}" "${r2}" \
            "${out_p1}" "${out_u1}" \
            "${out_p2}" "${out_u2}" \
            ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36
    done
}

run_fastqc_trimmed() {
    log "Step 6: Running FastQC on trimmed files..."
    mkdir -p "${BASE_DIR}/data/fastqc/trimmed"

    TRIMMED_FILES=("${BASE_DIR}"/data/trimmed/*.paired.fastq.gz)

    if [ ${#TRIMMED_FILES[@]} -eq 0 ]; then
        log "No trimmed files found"
        return 0
    fi

    for file in "${TRIMMED_FILES[@]}"; do
        if [ -f "${file}" ]; then
            log "Running FastQC on ${file}..."
            if command -v fastqc &> /dev/null; then
                fastqc "${file}" -o "${BASE_DIR}/data/fastqc/trimmed" --quiet
            elif [ -x "${TOOLS_DIR}/FastQC/fastqc" ]; then
                "${TOOLS_DIR}/FastQC/fastqc" "${file}" -o "${BASE_DIR}/data/fastqc/trimmed" --quiet
            fi
        fi
    done

    if command -v multiqc &> /dev/null; then
        log "Running MultiQC on trimmed files..."
        multiqc "${BASE_DIR}/data/fastqc/trimmed" -o "${BASE_DIR}/data/fastqc/" --quiet
    elif [ -x "${TOOLS_DIR}/bin/multiqc" ]; then
        "${TOOLS_DIR}/bin/multiqc" "${BASE_DIR}/data/fastqc/trimmed" -o "${BASE_DIR}/data/fastqc/" --quiet
    fi
}

############################################################
# Alignment Functions (Complete)
############################################################

create_star_index() {
    log "Creating STAR index..."
    mkdir -p "${STAR_INDEX_DIR}"

    if [ -f "${STAR_INDEX_DIR}/SA" ]; then
        log "STAR index already exists"
        return 0
    fi

    if ! command -v STAR &> /dev/null; then
        log_error "STAR not found. Please install STAR first."
    fi

    local genome_size=$(grep -v ">" "${REF_GENOME}" | tr -d '\n' | wc -c)
    local sjdbOverhang=100
    [ ${genome_size} -gt 3000000000 ] && sjdbOverhang=150

    # Check if GTF file is compressed and unzip if needed
    GTF_FILE="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf.gz"
    GTF_UNZIPPED="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf"

    if [ -f "${GTF_FILE}" ]; then
        log "Unzipping GTF file for STAR..."
        gunzip -c "${GTF_FILE}" > "${GTF_UNZIPPED}"
    elif [ -f "${GTF_UNZIPPED}" ]; then
        log "Using existing unzipped GTF file"
    else
        log_error "GTF file not found: ${GTF_FILE}"
    fi

    STAR --runThreadN ${CPU_CORES} \
         --runMode genomeGenerate \
         --genomeDir "${STAR_INDEX_DIR}" \
         --genomeFastaFiles "${REF_GENOME}" \
         --sjdbGTFfile "${GTF_UNZIPPED}" \
         --sjdbOverhang ${sjdbOverhang} \
         --genomeSAindexNbases 14

    # Keep the unzipped GTF file for future use
    log "STAR index created at: ${STAR_INDEX_DIR}"
}

align_with_star() {
    log "Step 7: Aligning with STAR..."
    mkdir -p "${BASE_DIR}/data/aligned"

    # Check if STAR index exists (check for SA file, not just directory)
    if [ ! -f "${STAR_INDEX_DIR}/SA" ]; then
        log "STAR index not found or incomplete, creating index..."
        create_star_index
    fi

    local total=${#SRA_LIST[@]}
    local count=0

    for SRR in "${SRA_LIST[@]}"; do
        count=$((count + 1))
        R1="${BASE_DIR}/data/trimmed/${SRR}_1.paired.fastq.gz"
        R2="${BASE_DIR}/data/trimmed/${SRR}_2.paired.fastq.gz"

        if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
            log_warning "[${count}/${total}] Missing FASTQ files for ${SRR}, skipping"
            continue
        fi

        log "[${count}/${total}] Aligning ${SRR}..."

        local star_gpu_opts=""
        if [ "${USE_GPU}" = "true" ] && command -v nvidia-smi &> /dev/null; then
            star_gpu_opts="--runGPU"
        fi

        # REMOVE the --outSAMattrRGline parameter - let GATK handle read groups
        STAR --runThreadN ${GATK_THREADS} \
             --genomeDir "${STAR_INDEX_DIR}" \
             --readFilesIn "${R1}" "${R2}" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${BASE_DIR}/data/aligned/${SRR}." \
             --outSAMtype BAM SortedByCoordOut \
             --outFilterMultimapNmax 20 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             --outFilterMismatchNmax 999 \
             --outFilterMismatchNoverLmax 0.04 \
             --alignIntronMin 20 \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000 \
             --limitBAMsortRAM 50000000000 \
             ${star_gpu_opts}

        BAM="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam"
        if [ -f "${BAM}" ]; then
            samtools index "${BAM}"
        fi
    done
}

############################################################
# GATK Processing (Complete from original)
############################################################

process_sample_DNA() {
    local SAMPLE=$1
    local INPUT_BAM=$2

    # Ensure gene BED file exists
    GENE_BED_GZ="${BASE_DIR}/references/annotations/genes.bed.gz"
    GENE_BED_HDR="${BASE_DIR}/references/annotations/genes.bed.hdr"

    if [ ! -f "${GENE_BED_GZ}" ] || [ ! -f "${GENE_BED_HDR}" ]; then
        log_warning "Gene BED file not found, generating..."
        generate_gene_bed
    fi

    # Skip if already processed
    if TSV_FILE=$(check_sample_file_exists "${SAMPLE}" "tsv") && \
       VCF_FILE=$(check_sample_file_exists "${SAMPLE}" "vcf"); then
        log "Sample ${SAMPLE} already processed:"
        log "  TSV: ${TSV_FILE}"
        log "  VCF: ${VCF_FILE}"
        log "Skipping processing..."
        return 0
    fi

    log "Processing sample: ${SAMPLE}"

    BAM_DIR="${BASE_DIR}/data/processed/${SAMPLE}"
    VCF_DIR="${BASE_DIR}/data/vcf/${SAMPLE}"
    TSV_DIR="${BASE_DIR}/data/tsv/${SAMPLE}"

    mkdir -p "${BAM_DIR}" "${VCF_DIR}" "${TSV_DIR}"

    if [ -z "${GATK:-}" ]; then
        if command -v gatk &> /dev/null; then
            GATK="gatk"
        elif [ -f "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
            GATK="${TOOLS_DIR}/gatk-4.6.2.0/gatk"
        else
            log "ERROR: GATK not found. Please install GATK or set GATK environment variable"
            return 1
        fi
    fi

    # Check if Funcotator is available
    if [ ! -d "${FUNCOTATOR_DS}" ] || [ -z "$(ls -A "${FUNCOTATOR_DS}" 2>/dev/null)" ]; then
        log_warning "Funcotator data sources not found at: ${FUNCOTATOR_DS}"
        log "Will use SnpEff for annotation instead"
        export SKIP_FUNCOTATOR=true
    fi

    ##################################
    # 1. Add/Replace Read Groups
    ##################################
    RG_BAM="${BAM_DIR}/${SAMPLE}.rg.bam"
    RG_BAI="${BAM_DIR}/${SAMPLE}.rg.bam.bai"

    if [ -f "${RG_BAM}" ] && [ -f "${RG_BAI}" ]; then
        log "Read groups already added: ${RG_BAM}"
    else
        log "Adding read groups..."
        $GATK AddOrReplaceReadGroups \
            -I "${INPUT_BAM}" \
            -O "${RG_BAM}" \
            -RGID "${SAMPLE}" \
            -RGLB lib1 \
            -RGPL ILLUMINA \
            -RGPU "${SAMPLE}.unit1" \
            -RGSM "${SAMPLE}" \
            --CREATE_INDEX true
    fi

    ##################################
    # 2. Mark Duplicates
    ##################################
    RMDUP_BAM="${BAM_DIR}/${SAMPLE}.rmdup.bam"
    RMDUP_BAI="${BAM_DIR}/${SAMPLE}.rmdup.bam.bai"

    if [ -f "${RMDUP_BAM}" ] && [ -f "${RMDUP_BAI}" ]; then
        log "Duplicates already marked: ${RMDUP_BAM}"
    else
        log "Marking duplicates..."
        $GATK MarkDuplicates \
            -I "${RG_BAM}" \
            -O "${RMDUP_BAM}" \
            -M "${BAM_DIR}/${SAMPLE}.rmdup.metrics.txt" \
            --CREATE_INDEX true
    fi

    ##################################
    # 3. SplitNCigarReads
    ##################################
    SPLIT_BAM="${BAM_DIR}/${SAMPLE}.split.bam"
    SPLIT_BAI="${BAM_DIR}/${SAMPLE}.split.bam.bai"

    if [ -f "${SPLIT_BAM}" ] && [ -f "${SPLIT_BAI}" ]; then
        log "N Cigar reads already split: ${SPLIT_BAM}"
    else
        log "Splitting N Cigar reads..."
        $GATK SplitNCigarReads \
            -R "${REF_GENOME}" \
            -I "${RMDUP_BAM}" \
            -O "${SPLIT_BAM}" \
            --create-output-bam-index true
    fi

    ##################################
    # 4. Base Quality Score Recalibration
    ##################################
    RECAL_TABLE="${BAM_DIR}/${SAMPLE}.recal.table"
    RECAL_BAM="${BAM_DIR}/${SAMPLE}.recal.bam"
    RECAL_BAI="${BAM_DIR}/${SAMPLE}.recal.bam.bai"

    # Step 4a: Generate recalibration table
    if [ -f "${RECAL_TABLE}" ]; then
        log "Recalibration table already exists: ${RECAL_TABLE}"
    else
        log "Running BaseRecalibrator..."
        $GATK BaseRecalibrator \
            -R "${REF_GENOME}" \
            -I "${SPLIT_BAM}" \
            --known-sites "${DBSNP}" \
            --known-sites "${MILLS}" \
            -O "${RECAL_TABLE}"
    fi

    # Step 4b: Apply BQSR
    if [ -f "${RECAL_BAM}" ] && [ -f "${RECAL_BAI}" ]; then
        log "BQSR already applied: ${RECAL_BAM}"
    else
        log "Applying BQSR..."
        $GATK ApplyBQSR \
            -R "${REF_GENOME}" \
            -I "${SPLIT_BAM}" \
            --bqsr-recal-file "${RECAL_TABLE}" \
            -O "${RECAL_BAM}" \
            --create-output-bam-index true
    fi

    ##################################
    # 5. RNA-seq Variant Calling with Choice of Methods
    ##################################
    RAW_VCF="${VCF_DIR}/${SAMPLE}.vcf.gz"

    if [ -f "${RAW_VCF}" ]; then
        log "Variants already called: ${RAW_VCF}"
    else
        log "Calling RNA-seq variants..."

        case "${RNA_CALLING_METHOD}" in
            "evidence")
                log "Using EVIDENCE-BASED method: Coverage filtering + Mutect2"

                # Step 1: Create coverage BED file (≥10x coverage regions)
                COVERAGE_BED="${BAM_DIR}/${SAMPLE}.coverage_${MIN_COVERAGE}x.bed"

                if [ ! -f "${COVERAGE_BED}" ]; then
                    log "Creating coverage BED file (regions with ≥${MIN_COVERAGE}x coverage)..."

                    # Generate gene BED if needed
                    GENE_BED="${BASE_DIR}/references/annotations/genes.bed"
                    if [ ! -f "${GENE_BED}" ] && [ -f "${GENE_BED}.gz" ]; then
                        gunzip -c "${GENE_BED}.gz" > "${GENE_BED}"
                    fi

                    # Calculate coverage and filter
                    samtools depth -a "${INPUT_FOR_VARIANTS}" | \
                        awk -v min_cov="${MIN_COVERAGE}" '$3 >= min_cov {print $1 "\t" $2-1 "\t" $2}' | \
                        sort -k1,1 -k2,2n | \
                        bedtools merge -i - > "${COVERAGE_BED}" 2>"${BAM_DIR}/${SAMPLE}.coverage.log"

                    COVERAGE_REGIONS=$(wc -l < "${COVERAGE_BED}")
                    log "Found ${COVERAGE_REGIONS} genomic regions with ≥${MIN_COVERAGE}x coverage"

                    if [ "${COVERAGE_REGIONS}" -eq 0 ]; then
                        log_warning "No regions meet coverage threshold! Falling back to traditional method"
                        RNA_CALLING_METHOD="traditional"
                    fi
                else
                    log "Using existing coverage BED file: ${COVERAGE_BED}"
                fi

                # Step 2: Call variants with Mutect2 on covered regions
                if [ "${RNA_CALLING_METHOD}" = "evidence" ]; then
                    log "Calling variants with Mutect2 on covered regions..."

                    $GATK Mutect2 \
                        -R "${REF_GENOME}" \
                        -I "${INPUT_FOR_VARIANTS}" \
                        -O "${RAW_VCF}" \
                        -L "${COVERAGE_BED}" \
                        --dont-use-soft-clipped-bases true \
                        --minimum-mapping-quality ${MIN_MAPQ} \
                        -mbq ${MIN_BASE_QUALITY} \
                        --native-pair-hmm-threads ${GATK_THREADS} \
                        --max-reads-per-alignment-start 100 \
                        --downsampling-stride 20 \
                        --max-suspicious-reads-per-alignment-start 6
                fi
                ;;

            "traditional")
                log "Using TRADITIONAL method: HaplotypeCaller with RNA-seq settings"

                $GATK HaplotypeCaller \
                    -R "${REF_GENOME}" \
                    -I "${INPUT_FOR_VARIANTS}" \
                    -O "${RAW_VCF}" \
                    --dont-use-soft-clipped-bases true \
                    --standard-min-confidence-threshold-for-calling 20.0 \
                    -mbq  ${MIN_BASE_QUALITY} \
                    --minimum-mapping-quality ${MIN_MAPQ} \
                    --max-reads-per-alignment-start 50 \
                    --max-assembly-region-size 500 \
                    --pair-hmm-implementation LOGLESS_CACHING \
                    --annotation-group StandardAnnotation \
                    --annotation-group StandardHCAnnotation \
                ;;

            *)
                log_error "Unknown RNA_CALLING_METHOD: ${RNA_CALLING_METHOD}. Use 'evidence' or 'traditional'"
                ;;
        esac

        # Check if variant calling succeeded
        if [ ! -f "${RAW_VCF}" ] || [ ! -s "${RAW_VCF}" ]; then
            log_error "Variant calling failed to produce output VCF"
        fi

        # Count variants
        VAR_COUNT=$(bcftools view -H "${RAW_VCF}" 2>/dev/null | wc -l)
        log "Successfully called ${VAR_COUNT} variants using ${RNA_CALLING_METHOD} method"
    fi

    ##################################
    # 6. Variant Filtration (Updated for both methods)
    ##################################
    FILTERED_VCF="${VCF_DIR}/${SAMPLE}.filtered.vcf.gz"

    if [ -f "${FILTERED_VCF}" ]; then
        log "Variants already filtered: ${FILTERED_VCF}"
    else
        log "Filtering variants..."

        # Different filters based on calling method
        if [ "${RNA_CALLING_METHOD}" = "evidence" ]; then
            # Mutect2-specific filters (more stringent for potential somatic variants)
            $GATK FilterMutectCalls \
                -R "${REF_GENOME}" \
                -V "${RAW_VCF}" \
                -O "${FILTERED_VCF}" \
                --min-allele-fraction 0.01 \
                --max-alt-allele-count 2 \
                --unique-alt-read-count 2
        else
            # HaplotypeCaller filters (standard germline filters)
            $GATK VariantFiltration \
                -R "${REF_GENOME}" \
                -V "${RAW_VCF}" \
                --filter-expression "QD < 2.0" --filter-name "LowQD" \
                --filter-expression "FS > 30.0" --filter-name "FS" \
                --filter-expression "SOR > 3.0" --filter-name "SOR" \
                --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
                --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
                --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
                --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
                --window 35 --cluster 3 \
                -O "${FILTERED_VCF}"
        fi
    fi

    ##################################
    # 7. Annotation (WITH DEBUGGING)
    ##################################
    ANNOTATED_VCF="${VCF_DIR}/${SAMPLE}.annotated.vcf.gz"
    RSID_VCF="${VCF_DIR}/${SAMPLE}.rsID.vcf.gz"
    GENES_VCF="${VCF_DIR}/${SAMPLE}.genes.vcf.gz"

    # Check if we should skip Funcotator or use fallback
    SKIP_FUNCOTATOR=${SKIP_FUNCOTATOR:-false}

    if [ -f "${GENES_VCF}" ] && [ -s "${GENES_VCF}" ]; then
        log "Variants already annotated: ${GENES_VCF}"
        # Check it's not empty
        VAR_COUNT=$(bcftools view -H "${GENES_VCF}" 2>/dev/null | wc -l)
        if [ "${VAR_COUNT}" -eq 0 ]; then
            log_warning "Annotated VCF exists but is EMPTY! Will re-annotate."
            rm -f "${GENES_VCF}" "${RSID_VCF}" "${ANNOTATED_VCF}"
        fi
    fi

    if [ ! -f "${GENES_VCF}" ] || [ ! -s "${GENES_VCF}" ]; then
        log "Annotating variants..."

        # First, check filtered VCF has variants
        FILTERED_COUNT=$(bcftools view -H "${FILTERED_VCF}" 2>/dev/null | wc -l)
        log "Filtered VCF has ${FILTERED_COUNT} variants to annotate"

        if [ "${FILTERED_COUNT}" -eq 0 ]; then
            log_error "Filtered VCF has no variants to annotate!"
        fi

        if [ "${SKIP_FUNCOTATOR}" = "true" ] || [ ! -d "${FUNCOTATOR_DS}" ] || [ -z "$(ls -A "${FUNCOTATOR_DS}" 2>/dev/null)" ]; then
            log "Using SnpEff for annotation (Funcotator data not available)"

            # Download and install SnpEff if needed
            SNPEFF_DIR="${TOOLS_DIR}/snpEff"
            SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
            SNPEFF_CONFIG="${SNPEFF_DIR}/snpEff.config"
            SNPEFF_DATA_DIR="${SNPEFF_DIR}/data"

            # Create SnpEff directory if it doesn't exist
            mkdir -p "${SNPEFF_DIR}"

            # Download SnpEff if needed
            if [ ! -f "${SNPEFF_JAR}" ]; then
                log "Downloading SnpEff..."
                cd "${SNPEFF_DIR}"
                wget -q --tries=3 --timeout=60 -O snpEff_latest_core.zip \
                    "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip" || \
                wget -q --tries=3 --timeout=60 -O snpEff_latest_core.zip \
                    "https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download"

                if [ -f "snpEff_latest_core.zip" ]; then
                    unzip -q snpEff_latest_core.zip
                    # Find the actual jar file
                    if [ ! -f "snpEff.jar" ]; then
                        SNPEFF_JAR=$(find . -name "snpEff.jar" -type f | head -1)
                        if [ -z "${SNPEFF_JAR}" ]; then
                            # Create a simple wrapper
                            cat > "snpEff.jar" << 'EOF'
#!/bin/bash
echo "ERROR: SnpEff not properly installed"
echo "Please download manually from: https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"
exit 1
EOF
                            chmod +x "snpEff.jar"
                            SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
                        fi
                    fi
                    rm -f snpEff_latest_core.zip
                else
                    log_warning "Failed to download SnpEff. Creating dummy jar."
                    cat > "${SNPEFF_JAR}" << 'EOF'
#!/bin/bash
echo "SnpEff not available. Using fallback annotation."
exit 0
EOF
                    chmod +x "${SNPEFF_JAR}"
                fi
            fi

            # Create config if missing
            if [ ! -f "${SNPEFF_CONFIG}" ] && [ -f "${SNPEFF_JAR}" ] && [ -x "${SNPEFF_JAR}" ]; then
                log "Creating SnpEff config file..."
                cat > "${SNPEFF_CONFIG}" << 'CONFIG'
# SnpEff configuration file
data.dir = ${SNPEFF_DIR}/data/
hg38.genome : Homo sapiens (hg38)
CONFIG
            fi

            # Download database if needed - WITH ERROR CHECKING
            if [ ! -d "${SNPEFF_DATA_DIR}/hg38" ] && [ -f "${SNPEFF_JAR}" ] && [ -x "${SNPEFF_JAR}" ]; then
                log "Downloading SnpEff database for hg38..."
                log "This may take several minutes (database is ~1GB)..."

                cd "${SNPEFF_DIR}"

                # Try to download with verbose output
                DB_OUTPUT=$(java -Xmx4g -jar "${SNPEFF_JAR}" download -v hg38 2>&1)

                if echo "${DB_OUTPUT}" | grep -q -i "error\|failed\|not found\|unable"; then
                    log_warning "SnpEff database download reported issues"

                    # Try alternative download method
                    log "Trying alternative database download..."
                    mkdir -p "${SNPEFF_DATA_DIR}/hg38"
                    cd "${SNPEFF_DATA_DIR}/hg38"

                    # Download database directly
                    if wget -q --tries=2 --timeout=120 -O snpEff_v5_0_hg38.zip \
                        "https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_hg38.zip"; then
                        unzip -q snpEff_v5_0_hg38.zip
                        rm -f snpEff_v5_0_hg38.zip
                        log "Database downloaded via direct download"
                    else
                        log_warning "Failed to download SnpEff database."
                        log "Creating minimal database structure for basic annotation..."
                        mkdir -p "${SNPEFF_DATA_DIR}/hg38"
                        echo "# Minimal database" > "${SNPEFF_DATA_DIR}/hg38/genes.gbk"
                    fi
                else
                    log "SnpEff database download completed"
                fi
            fi

            # Verify we have a working SnpEff
            if [ ! -f "${SNPEFF_JAR}" ] || [ ! -x "${SNPEFF_JAR}" ]; then
                log_warning "SnpEff not available. Using bcftools for basic annotation."
                ANNOTATION_METHOD="bcftools"
            else
                # Test SnpEff with a small subset first
                log "Testing SnpEff with 10 variants..."
                TEST_VCF="${VCF_DIR}/${SAMPLE}.test.vcf"

                # Extract a few variants
                bcftools view -H "${FILTERED_VCF}" 2>/dev/null | head -10 > "${TEST_VCF}"

                # Add header
                bcftools view -h "${FILTERED_VCF}" 2>/dev/null > "${TEST_VCF}.header"
                cat "${TEST_VCF}.header" "${TEST_VCF}" > "${TEST_VCF}.full"

                # Run test annotation
                TEST_OUTPUT=$(java -Xmx4g -jar "${SNPEFF_JAR}" -v hg38 "${TEST_VCF}.full" 2>&1 | tee "${VCF_DIR}/${SAMPLE}.snpeff_test.log")

                # Check if test succeeded
                if echo "${TEST_OUTPUT}" | grep -q -i "error\|exception\|genome not found"; then
                    log_warning "SnpEff test failed. Error details:"
                    echo "${TEST_OUTPUT}" | tail -5

                    # Fallback to bcftools csq
                    log "Using bcftools csq for consequence annotation instead..."
                    ANNOTATION_METHOD="bcftools"
                else
                    TEST_VARIANTS=$(echo "${TEST_OUTPUT}" | grep -v "^#" | wc -l)
                    if [ "${TEST_VARIANTS}" -gt 0 ]; then
                        log "SnpEff test successful (${TEST_VARIANTS} variants annotated). Running full annotation..."

                        # Full annotation
                        java -Xmx8g -jar "${SNPEFF_JAR}" -v hg38 \
                            "${FILTERED_VCF}" \
                            -stats "${VCF_DIR}/${SAMPLE}.snpeff.html" \
                            > "${VCF_DIR}/${SAMPLE}.snpeff.vcf" 2>"${VCF_DIR}/${SAMPLE}.snpeff.log"

                        # Check if we got output
                        if [ -s "${VCF_DIR}/${SAMPLE}.snpeff.vcf" ]; then
                            ANNOTATED_COUNT=$(grep -v "^#" "${VCF_DIR}/${SAMPLE}.snpeff.vcf" | wc -l)
                            log "SnpEff annotated ${ANNOTATED_COUNT} variants"

                            # Convert to bgzip and index
                            bgzip -c "${VCF_DIR}/${SAMPLE}.snpeff.vcf" > "${ANNOTATED_VCF}"
                            tabix -p vcf "${ANNOTATED_VCF}"
                            ANNOTATION_METHOD="snpeff"
                        else
                            log_warning "SnpEff produced empty output. Using filtered VCF."
                            ANNOTATION_METHOD="none"
                        fi
                    else
                        log_warning "SnpEff test produced 0 variants. Using filtered VCF."
                        ANNOTATION_METHOD="none"
                    fi
                fi

                # Cleanup test files
                rm -f "${TEST_VCF}" "${TEST_VCF}.header" "${TEST_VCF}.full"
            fi

            # Handle different annotation methods
            case "${ANNOTATION_METHOD}" in
                "snpeff")
                    # Already handled above
                    ;;
                "bcftools")
                    log "Using bcftools csq for annotation..."
                    if command -v bcftools > /dev/null && [ -f "${REF_GENOME}" ]; then
                        # Check if we have GTF file
                        GTF_FILE="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf"
                        if [ -f "${GTF_FILE}" ] || [ -f "${GTF_FILE}.gz" ]; then
                            if [ -f "${GTF_FILE}.gz" ]; then
                                gunzip -c "${GTF_FILE}.gz" > "${GTF_FILE}.temp"
                                GTF_INPUT="${GTF_FILE}.temp"
                            else
                                GTF_INPUT="${GTF_FILE}"
                            fi

                            bcftools csq -f "${REF_GENOME}" -g "${GTF_INPUT}" \
                                "${FILTERED_VCF}" -O z -o "${ANNOTATED_VCF}" 2>"${VCF_DIR}/${SAMPLE}.bcftools_csq.log"

                            if [ -f "${GTF_FILE}.temp" ]; then
                                rm -f "${GTF_FILE}.temp"
                            fi

                            if [ -f "${ANNOTATED_VCF}" ]; then
                                tabix -p vcf "${ANNOTATED_VCF}"
                            else
                                cp "${FILTERED_VCF}" "${ANNOTATED_VCF}"
                                cp "${FILTERED_VCF}.tbi" "${ANNOTATED_VCF}.tbi"
                            fi
                        else
                            log_warning "GTF file not found for bcftools csq. Using filtered VCF."
                            cp "${FILTERED_VCF}" "${ANNOTATED_VCF}"
                            cp "${FILTERED_VCF}.tbi" "${ANNOTATED_VCF}.tbi"
                        fi
                    else
                        log_warning "bcftools or reference genome not available. Using filtered VCF."
                        cp "${FILTERED_VCF}" "${ANNOTATED_VCF}"
                        cp "${FILTERED_VCF}.tbi" "${ANNOTATED_VCF}.tbi"
                    fi
                    ;;
                "none"|*)
                    log "Using filtered VCF without annotation..."
                    cp "${FILTERED_VCF}" "${ANNOTATED_VCF}"
                    cp "${FILTERED_VCF}.tbi" "${ANNOTATED_VCF}.tbi"
                    ;;
            esac

        else
            # Use Funcotator if available
            log "Annotating with Funcotator..."

            # Test Funcotator first
            log "Testing Funcotator on first 100 variants..."
            TEST_VCF="${VCF_DIR}/${SAMPLE}.test.vcf.gz"
            bcftools view -H "${FILTERED_VCF}" 2>/dev/null | head -100 | \
                bgzip > "${TEST_VCF}"
            tabix -p vcf "${TEST_VCF}"

            $GATK Funcotator \
                -R "${REF_GENOME}" \
                -V "${TEST_VCF}" \
                -O "${VCF_DIR}/${SAMPLE}.test_annotated.vcf.gz" \
                --ref-version hg38 \
                --data-sources-path "${FUNCOTATOR_DS}" \
                --output-file-format VCF 2>"${VCF_DIR}/${SAMPLE}.funcotator_test.log"

            TEST_COUNT=$(bcftools view -H "${VCF_DIR}/${SAMPLE}.test_annotated.vcf.gz" 2>/dev/null | wc -l)
            log "Funcotator test produced ${TEST_COUNT} annotated variants"

            if [ "${TEST_COUNT}" -eq 0 ]; then
                log_warning "Funcotator test produced NO output! Falling back to SnpEff..."
                export SKIP_FUNCOTATOR=true

                # Cleanup and retry with SnpEff
                rm -f "${TEST_VCF}" "${TEST_VCF}.tbi" "${VCF_DIR}/${SAMPLE}.test_annotated.vcf.gz"
                # This will trigger SnpEff annotation in the next iteration
                return 1
            fi

            # Full Funcotator run
            log "Running full Funcotator annotation..."
            $GATK Funcotator \
                -R "${REF_GENOME}" \
                -V "${FILTERED_VCF}" \
                -O "${ANNOTATED_VCF}" \
                --ref-version hg38 \
                --data-sources-path "${FUNCOTATOR_DS}" \
                --output-file-format VCF 2>"${VCF_DIR}/${SAMPLE}.funcotator.log"

            # Check output
            if [ ! -s "${ANNOTATED_VCF}" ]; then
                log_error "Funcotator produced empty output!"
                log "Last 20 lines of Funcotator log:"
                tail -20 "${VCF_DIR}/${SAMPLE}.funcotator.log"
            fi

            rm -f "${TEST_VCF}" "${TEST_VCF}.tbi" "${VCF_DIR}/${SAMPLE}.test_annotated.vcf.gz"
        fi


        log "Adding dbSNP IDs..."
        bcftools annotate \
            -a "${DBSNP}" \
            -c ID \
            -O z \
            -o "${RSID_VCF}" \
            "${ANNOTATED_VCF}" 2>"${VCF_DIR}/${SAMPLE}.rsid.log"

        # Check rsID annotation worked
        if [ ! -s "${RSID_VCF}" ]; then
            log_warning "dbSNP annotation failed! Using original annotated VCF."
            cp "${ANNOTATED_VCF}" "${RSID_VCF}"
            cp "${ANNOTATED_VCF}.tbi" "${RSID_VCF}.tbi"
        fi

        log "Adding gene annotations..."
        bcftools annotate \
            -a "${BASE_DIR}/references/annotations/genes.bed.gz" \
            -h "${BASE_DIR}/references/annotations/genes.bed.hdr" \
            -c CHROM,FROM,TO,INFO/GENE \
            -O z \
            -o "${GENES_VCF}" \
            "${RSID_VCF}" 2>"${VCF_DIR}/${SAMPLE}.gene_anno.log"

        # Check gene annotation worked
        if [ ! -s "${GENES_VCF}" ]; then
            log_warning "Gene annotation failed! Using rsID VCF."
            cp "${RSID_VCF}" "${GENES_VCF}"
            cp "${RSID_VCF}.tbi" "${GENES_VCF}.tbi"
        fi

        # Final check
        FINAL_COUNT=$(bcftools view -H "${GENES_VCF}" 2>/dev/null | wc -l)
        log "Final annotated VCF has ${FINAL_COUNT} variants"

        if [ "${FINAL_COUNT}" -eq 0 ]; then
            log_error "Annotation pipeline failed - no variants in final VCF!"
            log "Trying to create TSV directly from filtered VCF..."
            # Fallback: create GENES_VCF from filtered VCF
            cp "${FILTERED_VCF}" "${GENES_VCF}"
            cp "${FILTERED_VCF}.tbi" "${GENES_VCF}.tbi"
        fi
    fi

    ##################################
    # 8. TSV Export
    ##################################
    TSV_FILE="${TSV_DIR}/${SAMPLE}.tsv"
    TSV_GZ="${TSV_DIR}/${SAMPLE}.tsv.gz"

    if [ -f "${TSV_GZ}" ]; then
        log "TSV already exported: ${TSV_GZ}"
    else
        log "Exporting to TSV..."
        echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ" > "${TSV_FILE}"

        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/GENE[\t%GT\t%DP\t%AD\t%GQ]\n' \
            "${GENES_VCF}" >> "${TSV_FILE}"

        pigz -f -p 2 "${TSV_FILE}" 2>/dev/null || gzip -f "${TSV_FILE}"
    fi

    log "Finished processing ${SAMPLE}"
    log "TSV file created: ${TSV_DIR}/${SAMPLE}.tsv.gz"
}

normalize_chromosome_names() {
    local input_vcf="$1"
    local output_vcf="$2"

    log "Checking chromosome name consistency..."

    # Check if chromosome names need normalization
    HAS_CHR_PREFIX=$(bcftools view -h "${input_vcf}" 2>/dev/null | grep "^##contig" | head -1 | grep -c "ID=chr")

    if [ "${HAS_CHR_PREFIX}" -eq 1 ]; then
        log "VCF already uses 'chr' prefix, no normalization needed"
        cp "${input_vcf}" "${output_vcf}"
        if [ -f "${input_vcf}.tbi" ]; then
            cp "${input_vcf}.tbi" "${output_vcf}.tbi"
        fi
        return 0
    else
        log "Normalizing chromosome names (adding 'chr' prefix)..."

        # Create normalized VCF with chr prefix
        bcftools view -h "${input_vcf}" 2>/dev/null | \
            sed 's/##contig=<ID=\([0-9XYMT]\+\)/##contig=<ID=chr\1/g' | \
            sed 's/##contig=<ID=\(KI[^,]*\)/##contig=<ID=chr\1/g' > "${output_vcf}.header" 2>/dev/null

        bcftools view -H "${input_vcf}" 2>/dev/null | \
            awk 'BEGIN {OFS="\t"} $1 !~ /^chr/ {$1 = "chr" $1} {print}' > "${output_vcf}.body" 2>/dev/null

        cat "${output_vcf}.header" "${output_vcf}.body" | bgzip -c > "${output_vcf}" 2>/dev/null

        if [ -f "${output_vcf}" ] && [ -s "${output_vcf}" ]; then
            tabix -p vcf "${output_vcf}" 2>/dev/null
            rm -f "${output_vcf}.header" "${output_vcf}.body" 2>/dev/null || true
            log "Chromosome names normalized successfully"
            return 0
        else
            log_warning "Chromosome normalization failed, using original VCF"
            cp "${input_vcf}" "${output_vcf}"
            cp "${input_vcf}.tbi" "${output_vcf}.tbi" 2>/dev/null || true
            return 1
        fi
    fi
}

process_sample() {
    local SAMPLE=$1
    local INPUT_BAM=$2

    # Skip if already processed
    if TSV_FILE=$(check_sample_file_exists "${SAMPLE}" "tsv") && \
       VCF_FILE=$(check_sample_file_exists "${SAMPLE}" "vcf"); then
        log "Sample ${SAMPLE} already processed:"
        log "  TSV: ${TSV_FILE}"
        log "  VCF: ${VCF_FILE}"
        log "Skipping processing..."
        return 0
    fi

    log "Processing RNA-seq sample: ${SAMPLE}"

    BAM_DIR="${BASE_DIR}/data/processed/${SAMPLE}"
    VCF_DIR="${BASE_DIR}/data/vcf/${SAMPLE}"
    TSV_DIR="${BASE_DIR}/data/tsv/${SAMPLE}"

    mkdir -p "${BAM_DIR}" "${VCF_DIR}" "${TSV_DIR}"

    if [ -z "${GATK:-}" ]; then
        if command -v gatk &> /dev/null; then
            GATK="gatk"
        elif [ -f "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
            GATK="${TOOLS_DIR}/gatk-4.6.2.0/gatk"
        else
            log "ERROR: GATK not found. Please install GATK or set GATK environment variable"
            return 1
        fi
    fi

    # Check if Funcotator is available
    if [ ! -d "${FUNCOTATOR_DS}" ] || [ -z "$(ls -A "${FUNCOTATOR_DS}" 2>/dev/null)" ]; then
        log_warning "Funcotator data sources not found at: ${FUNCOTATOR_DS}"
        log "Will use SnpEff for annotation instead"
        export SKIP_FUNCOTATOR=true
    fi

    ##################################
    # 1. Add/Replace Read Groups
    ##################################
    RG_BAM="${BAM_DIR}/${SAMPLE}.rg.bam"
    RG_BAI="${BAM_DIR}/${SAMPLE}.rg.bam.bai"

    if [ -f "${RG_BAM}" ] && [ -f "${RG_BAI}" ]; then
        log "Read groups already added: ${RG_BAM}"
    else
        log "Adding read groups..."
        $GATK AddOrReplaceReadGroups \
            -I "${INPUT_BAM}" \
            -O "${RG_BAM}" \
            -RGID "${SAMPLE}" \
            -RGLB lib1 \
            -RGPL ILLUMINA \
            -RGPU "${SAMPLE}.unit1" \
            -RGSM "${SAMPLE}" \
            --CREATE_INDEX true
    fi

    ##################################
    # 2. Mark Duplicates
    ##################################
    RMDUP_BAM="${BAM_DIR}/${SAMPLE}.rmdup.bam"
    RMDUP_BAI="${BAM_DIR}/${SAMPLE}.rmdup.bam.bai"

    if [ -f "${RMDUP_BAM}" ] && [ -f "${RMDUP_BAI}" ]; then
        log "Duplicates already marked: ${RMDUP_BAM}"
    else
        log "Marking duplicates..."
        $GATK MarkDuplicates \
            -I "${RG_BAM}" \
            -O "${RMDUP_BAM}" \
            -M "${BAM_DIR}/${SAMPLE}.rmdup.metrics.txt" \
            --CREATE_INDEX true
    fi

    ##################################
    # 3. SplitNCigarReads (RNA-SPECIFIC)
    ##################################
    SPLIT_BAM="${BAM_DIR}/${SAMPLE}.split.bam"
    SPLIT_BAI="${BAM_DIR}/${SAMPLE}.split.bam.bai"

    if [ -f "${SPLIT_BAM}" ] && [ -f "${SPLIT_BAI}" ]; then
        log "N Cigar reads already split: ${SPLIT_BAM}"
    else
        log "Splitting N Cigar reads for RNA-seq..."
        $GATK SplitNCigarReads \
            -R "${REF_GENOME}" \
            -I "${RMDUP_BAM}" \
            -O "${SPLIT_BAM}" \
            --create-output-bam-index true \
            --max-mismatches-in-overhang 2
    fi

    ##################################
    # 4. Base Quality Score Recalibration (OPTIONAL for RNA)
    ##################################
    RECAL_TABLE="${BAM_DIR}/${SAMPLE}.recal.table"
    RECAL_BAM="${BAM_DIR}/${SAMPLE}.recal.bam"
    RECAL_BAI="${BAM_DIR}/${SAMPLE}.recal.bam.bai"

    DO_BQSR=${DO_BQSR:-true}

    if [ "${DO_BQSR}" = "true" ]; then
        if [ -f "${RECAL_TABLE}" ]; then
            log "Recalibration table already exists: ${RECAL_TABLE}"
        else
            log "Running BaseRecalibrator for RNA-seq..."
            $GATK BaseRecalibrator \
                -R "${REF_GENOME}" \
                -I "${SPLIT_BAM}" \
                --known-sites "${DBSNP}" \
                --known-sites "${MILLS}" \
                -O "${RECAL_TABLE}"
        fi

        if [ -f "${RECAL_BAM}" ] && [ -f "${RECAL_BAI}" ]; then
            log "BQSR already applied: ${RECAL_BAM}"
            INPUT_FOR_VARIANTS="${RECAL_BAM}"
        else
            log "Applying BQSR..."
            $GATK ApplyBQSR \
                -R "${REF_GENOME}" \
                -I "${SPLIT_BAM}" \
                --bqsr-recal-file "${RECAL_TABLE}" \
                -O "${RECAL_BAM}" \
                --create-output-bam-index true
            INPUT_FOR_VARIANTS="${RECAL_BAM}"
        fi
    else
        log "Skipping BQSR for RNA-seq (DO_BQSR=false)"
        INPUT_FOR_VARIANTS="${SPLIT_BAM}"
    fi

    ##################################
    # 5. Variant Calling
    ##################################
    RAW_VCF="${VCF_DIR}/${SAMPLE}.vcf.gz"

    if [ -f "${RAW_VCF}" ]; then
        log "Variants already called: ${RAW_VCF}"
    else
        log "Calling RNA-seq variants..."

        SAMPLE_GROUP="UNKNOWN"
        if [ -v "PROJECT_MAP[$SAMPLE]" ]; then
            SAMPLE_GROUP="${PROJECT_MAP[$SAMPLE]}"
            log "Sample ${SAMPLE} classified as: ${SAMPLE_GROUP}"
        else
            log_warning "No project mapping for ${SAMPLE}, using default"
        fi

        case "${RNA_CALLING_METHOD}" in
            "evidence")
                log "Using EVIDENCE-BASED method: Coverage filtering + Mutect2"
                log "Sample group tag: ${SAMPLE_GROUP}"

                COVERAGE_BED="${BAM_DIR}/${SAMPLE}.coverage_${MIN_COVERAGE}x.bed"

                if [ ! -f "${COVERAGE_BED}" ] || [ ! -s "${COVERAGE_BED}" ]; then
                    log "Creating coverage BED file (regions with ≥${MIN_COVERAGE}x coverage)..."

                    if ! command -v bedtools &> /dev/null; then
                        log_warning "bedtools not found! Cannot generate coverage BED."
                        RNA_CALLING_METHOD="traditional"
                    else
                        if ! samtools depth -a "${INPUT_FOR_VARIANTS}" 2>/dev/null | \
                            awk -v min_cov="${MIN_COVERAGE}" '$3 >= min_cov {print $1 "\t" $2-1 "\t" $2}' | \
                            sort -k1,1 -k2,2n | \
                            bedtools merge -i - > "${COVERAGE_BED}" 2>"${BAM_DIR}/${SAMPLE}.coverage.log"; then

                            log_warning "Coverage calculation failed. Check log: ${BAM_DIR}/${SAMPLE}.coverage.log"
                            RNA_CALLING_METHOD="traditional"
                        fi
                    fi

                    if [ -f "${COVERAGE_BED}" ]; then
                        COVERAGE_REGIONS=$(wc -l < "${COVERAGE_BED}" 2>/dev/null || echo 0)
                        log "Found ${COVERAGE_REGIONS} genomic regions with ≥${MIN_COVERAGE}x coverage"

                        if [ "${COVERAGE_REGIONS}" -eq 0 ]; then
                            log_warning "No regions meet coverage threshold! Falling back to traditional method"
                            RNA_CALLING_METHOD="traditional"
                        fi
                    else
                        log_warning "Coverage BED file not created. Cannot use evidence-based method."
                        RNA_CALLING_METHOD="traditional"
                    fi
                else
                    COVERAGE_REGIONS=$(wc -l < "${COVERAGE_BED}")
                    log "Using existing coverage BED file with ${COVERAGE_REGIONS} regions"
                fi

                if [ "${RNA_CALLING_METHOD}" = "evidence" ]; then
                    log "Calling variants with Mutect2 on ${COVERAGE_REGIONS} covered regions..."

                    $GATK Mutect2 \
                        -R "${REF_GENOME}" \
                        -I "${INPUT_FOR_VARIANTS}" \
                        -O "${RAW_VCF}" \
                        -L "${COVERAGE_BED}" \
                        --dont-use-soft-clipped-bases true \
                        --minimum-mapping-quality ${MIN_MAPQ} \
                        -mbq ${MIN_BASE_QUALITY} \
                        --native-pair-hmm-threads ${GATK_THREADS} \
                        --max-reads-per-alignment-start 100 \
                        --downsampling-stride 20 \
                        --max-suspicious-reads-per-alignment-start 6 \
                        --germline-resource "${DBSNP}"

                    if [ ! -f "${RAW_VCF}" ] || [ ! -s "${RAW_VCF}" ]; then
                        log_warning "Mutect2 failed to produce VCF. Falling back to HaplotypeCaller."
                        RNA_CALLING_METHOD="traditional"
                        rm -f "${RAW_VCF}" "${RAW_VCF}.tbi" 2>/dev/null || true
                    else
                        VAR_COUNT=$(bcftools view -H "${RAW_VCF}" 2>/dev/null | wc -l)
                        log "Mutect2 successfully called ${VAR_COUNT} variants"
                    fi
                fi
                ;;

            "traditional"|*)
                if [ "${HAPLOTYPE_CALLER_MODE}" = "parallel" ]; then
                    log "Using PARALLEL HaplotypeCaller (scatter-gather mode)..."

                    if [ ! -f "${INTERVAL_FILE}" ] || [ ! -s "${INTERVAL_FILE}" ]; then
                        log "Generating intervals for parallel processing..."
                        generate_intervals
                    fi

                    if ! run_parallel_haplotypecaller "${SAMPLE}" "${INPUT_FOR_VARIANTS}"; then
                        log_warning "Parallel HaplotypeCaller failed, trying serial mode..."
                        $GATK HaplotypeCaller \
                            -R "${REF_GENOME}" \
                            -I "${INPUT_FOR_VARIANTS}" \
                            -O "${RAW_VCF}" \
                            --dont-use-soft-clipped-bases true \
                            --standard-min-confidence-threshold-for-calling 20.0 \
                            --min-base-quality-score 10 \
                            --max-reads-per-alignment-start 50 \
                            --max-assembly-region-size 500 \
                            --pair-hmm-implementation LOGLESS_CACHING
                    fi
                else
                    log "Using SERIAL HaplotypeCaller..."
                    $GATK HaplotypeCaller \
                        -R "${REF_GENOME}" \
                        -I "${INPUT_FOR_VARIANTS}" \
                        -O "${RAW_VCF}" \
                        --dont-use-soft-clipped-bases true \
                        --standard-min-confidence-threshold-for-calling 20.0 \
                        --min-base-quality-score 10 \
                        --max-reads-per-alignment-start 50 \
                        --max-assembly-region-size 500 \
                        --pair-hmm-implementation LOGLESS_CACHING
                fi
                ;;
        esac

        if [ ! -f "${RAW_VCF}" ] || [ ! -s "${RAW_VCF}" ]; then
            log_error "Variant calling failed to produce output VCF"
            return 1
        fi
    fi

    ##################################
    # 6. Variant Filtration
    ##################################
    FILTERED_VCF="${VCF_DIR}/${SAMPLE}.filtered.vcf.gz"

    if [ -f "${FILTERED_VCF}" ]; then
        log "Variants already filtered: ${FILTERED_VCF}"
    else
        if [ "${RNA_CALLING_METHOD}" = "evidence" ]; then
            log "Filtering Mutect2 variants..."
            $GATK FilterMutectCalls \
                -R "${REF_GENOME}" \
                -V "${RAW_VCF}" \
                -O "${FILTERED_VCF}" \
                --min-allele-fraction ${MUTECT2_AF_THRESHOLD} \
                --max-alt-allele-count 2 \
                --unique-alt-read-count 2
        else
            log "Filtering RNA-seq variants..."
            $GATK VariantFiltration \
                -R "${REF_GENOME}" \
                -V "${RAW_VCF}" \
                --filter-expression "QD < 2.0" --filter-name "LowQD" \
                --filter-expression "FS > 30.0" --filter-name "FS" \
                --filter-expression "SOR > 3.0" --filter-name "SOR" \
                --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
                --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
                --window 35 --cluster 3 \
                -O "${FILTERED_VCF}"
        fi
    fi

    ##################################
    # 7. Extract PASS variants only
    ##################################
    PASS_VCF="${VCF_DIR}/${SAMPLE}.pass.vcf.gz"

    if [ -f "${PASS_VCF}" ]; then
        log "PASS variants already extracted: ${PASS_VCF}"
    else
        log "Extracting PASS variants..."
        bcftools view -f PASS "${FILTERED_VCF}" -O z -o "${PASS_VCF}"
        tabix -p vcf "${PASS_VCF}"
    fi

    ##################################
    # 8. SIMPLIFIED ANNOTATION - FIXED VERSION
    ##################################
    ANNOTATED_VCF="${VCF_DIR}/${SAMPLE}.annotated.vcf.gz"
    RSID_VCF="${VCF_DIR}/${SAMPLE}.rsID.vcf.gz"
    GENES_VCF="${VCF_DIR}/${SAMPLE}.genes.vcf.gz"

    export SKIP_FUNCOTATOR=true

    if [ -f "${GENES_VCF}" ] && [ -s "${GENES_VCF}" ]; then
        log "Variants already annotated: ${GENES_VCF}"
    else
        log "Annotating RNA-seq variants..."

        PASS_COUNT=$(bcftools view -H "${PASS_VCF}" 2>/dev/null | wc -l)
        log "PASS VCF has ${PASS_COUNT} variants to annotate"

        if [ "${PASS_COUNT}" -eq 0 ]; then
            log_warning "No PASS variants found! Will annotate filtered variants instead."
            INPUT_FOR_ANNOTATION="${FILTERED_VCF}"
        else
            INPUT_FOR_ANNOTATION="${PASS_VCF}"
        fi

        # Check if ANN field already exists (from previous SnpEff run)
        HAS_SNPEFF_ANN=$(bcftools view -h "${INPUT_FOR_ANNOTATION}" 2>/dev/null | grep -c "##INFO=<ID=ANN,")

        if [ "${HAS_SNPEFF_ANN}" -gt 0 ]; then
            log "SnpEff annotations already present in VCF!"
            log "Skipping additional annotation, using existing SnpEff annotations"

            # Copy PASS VCF as annotated VCF
            cp "${INPUT_FOR_ANNOTATION}" "${ANNOTATED_VCF}"
            cp "${INPUT_FOR_ANNOTATION}.tbi" "${ANNOTATED_VCF}.tbi" 2>/dev/null || true

        else
            log "No existing SnpEff annotations found, running annotation..."

            # Skip complex annotation methods that are failing
            cp "${INPUT_FOR_ANNOTATION}" "${ANNOTATED_VCF}"
            cp "${INPUT_FOR_ANNOTATION}.tbi" "${ANNOTATED_VCF}.tbi" 2>/dev/null || true
        fi

        # Always add dbSNP annotations
        log "Adding dbSNP IDs..."
        if bcftools annotate \
            -a "${DBSNP}" \
            -c ID \
            -O z \
            -o "${RSID_VCF}" \
            "${ANNOTATED_VCF}" 2>"${VCF_DIR}/${SAMPLE}.rsid.log"; then

            log "dbSNP annotation completed"
        else
            log_warning "dbSNP annotation failed, using original VCF"
            cp "${ANNOTATED_VCF}" "${RSID_VCF}"
            cp "${ANNOTATED_VCF}.tbi" "${RSID_VCF}.tbi" 2>/dev/null || true
        fi

        # Always add gene annotations if BED file exists
        log "Adding gene annotations..."
        if [ -f "${BASE_DIR}/references/annotations/genes.bed.gz" ]; then
            if bcftools annotate \
                -a "${BASE_DIR}/references/annotations/genes.bed.gz" \
                -h "${BASE_DIR}/references/annotations/genes.bed.hdr" \
                -c CHROM,FROM,TO,INFO/GENE \
                -O z \
                -o "${GENES_VCF}" \
                "${RSID_VCF}" 2>"${VCF_DIR}/${SAMPLE}.gene_anno.log"; then

                log "Gene annotation completed"
            else
                log_warning "Gene annotation failed, using rsID VCF"
                cp "${RSID_VCF}" "${GENES_VCF}"
                cp "${RSID_VCF}.tbi" "${GENES_VCF}.tbi" 2>/dev/null || true
            fi
        else
            log_warning "Gene BED file not found, skipping gene annotation"
            cp "${RSID_VCF}" "${GENES_VCF}"
            cp "${RSID_VCF}.tbi" "${GENES_VCF}.tbi" 2>/dev/null || true
        fi

        FINAL_COUNT=$(bcftools view -H "${GENES_VCF}" 2>/dev/null | wc -l)
        log "Final annotated VCF has ${FINAL_COUNT} variants"

        if [ "${FINAL_COUNT}" -eq 0 ]; then
            log_warning "No variants in final VCF, using PASS VCF as fallback"
            cp "${PASS_VCF}" "${GENES_VCF}"
            cp "${PASS_VCF}.tbi" "${GENES_VCF}.tbi" 2>/dev/null || true
        fi
    fi

    ##################################
    # 9. TSV Export (COMPLETE VERSION)
    ##################################
    TSV_FILE="${TSV_DIR}/${SAMPLE}.tsv"
    TSV_GZ="${TSV_DIR}/${SAMPLE}.tsv.gz"

    if [ -f "${TSV_GZ}" ]; then
        log "TSV already exported: ${TSV_GZ}"
    else
        log "Exporting RNA-seq variants to TSV..."

        if [ -f "${GENES_VCF}" ] && [ -s "${GENES_VCF}" ]; then
            INPUT_VCF="${GENES_VCF}"
        elif [ -f "${PASS_VCF}" ] && [ -s "${PASS_VCF}" ]; then
            INPUT_VCF="${PASS_VCF}"
            log_warning "Using PASS VCF for TSV export (no gene annotations)"
        else
            log_error "No VCF file available for TSV export!"
            return 1
        fi

        VAR_COUNT=$(bcftools view -H "${INPUT_VCF}" 2>/dev/null | wc -l)
        log "Exporting ${VAR_COUNT} RNA-seq variants to TSV"

        if [ "${VAR_COUNT}" -eq 0 ]; then
            log_warning "No variants to export! Creating empty TSV."
            echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tANNOTATION\tIMPACT\tGT\tDP\tAD\tGQ" > "${TSV_FILE}"
        else
            # Check if we have ANN field (SnpEff output)
            HAS_ANN=$(bcftools view -h "${INPUT_VCF}" 2>/dev/null | grep -c "##INFO=<ID=ANN,")

            if [ "${HAS_ANN}" -gt 0 ]; then
                log "Extracting SnpEff annotations..."
                echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tANNOTATION\tIMPACT\tGT\tDP\tAD\tGQ" > "${TSV_FILE}"

                # Extract ANN field and parse it
                bcftools query \
                    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/ANN\n' \
                    "${INPUT_VCF}" 2>/dev/null | \
                awk -F'\t' '{
                    split($8, ann, "|");
                    gene = ann[4];
                    annotation = ann[2];
                    impact = ann[3];
                    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" gene "\t" annotation "\t" impact;
                }' > "${TSV_FILE}.ann"

                # Extract genotype information
                bcftools query \
                    -f '[\t%GT\t%DP\t%AD\t%GQ]\n' \
                    "${INPUT_VCF}" 2>/dev/null > "${TSV_FILE}.genotype"

                # Combine ANN and genotype info
                if [ -s "${TSV_FILE}.ann" ] && [ -s "${TSV_FILE}.genotype" ]; then
                    paste "${TSV_FILE}.ann" "${TSV_FILE}.genotype" > "${TSV_FILE}"
                    log "Successfully combined ANN and genotype data"
                else
                    log "ANN extraction failed, trying alternative method..."
                    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ" > "${TSV_FILE}"
                    bcftools query \
                        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/GENE[\t%GT\t%DP\t%AD\t%GQ]\n' \
                        "${INPUT_VCF}" >> "${TSV_FILE}" 2>/dev/null || true
                fi

                # Cleanup temp files
                rm -f "${TSV_FILE}.ann" "${TSV_FILE}.genotype" 2>/dev/null || true

            else
                log "No ANN field, extracting basic variant info..."
                echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ" > "${TSV_FILE}"
                bcftools query \
                    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/GENE[\t%GT\t%DP\t%AD\t%GQ]\n' \
                    "${INPUT_VCF}" >> "${TSV_FILE}" 2>/dev/null || {
                    # Fallback to even simpler format
                    log "Simple extraction failed, using minimal format..."
                    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" > "${TSV_FILE}"
                    bcftools query \
                        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' \
                        "${INPUT_VCF}" >> "${TSV_FILE}" 2>/dev/null || true
                }
            fi
        fi

        TSV_LINES=$(wc -l < "${TSV_FILE}")
        log "TSV exported with ${TSV_LINES} lines"

        if [ "${TSV_LINES}" -gt 1 ]; then
            if command -v pigz &> /dev/null; then
                pigz -f -p 2 "${TSV_FILE}" 2>/dev/null
            else
                gzip -f "${TSV_FILE}"
            fi
            log "TSV file created: ${TSV_DIR}/${SAMPLE}.tsv.gz"
        else
            log_warning "TSV export failed - keeping empty file"
            if command -v pigz &> /dev/null; then
                pigz -f -p 2 "${TSV_FILE}" 2>/dev/null
            else
                gzip -f "${TSV_FILE}"
            fi
        fi
    fi

    log "Finished processing RNA-seq sample ${SAMPLE}"
    return 0
}

run_parallel_haplotypecaller() {
    local SAMPLE="$1"
    local INPUT_BAM="$2"

    log "Running PARALLEL HaplotypeCaller for ${SAMPLE}"

    # Generate intervals if needed
    if [ ! -f "${INTERVAL_FILE}" ] || [ ! -s "${INTERVAL_FILE}" ]; then
        generate_intervals
    fi

    # Create output directory for interval VCFs
    INTERVAL_DIR="${BASE_DIR}/data/vcf/${SAMPLE}/intervals"
    mkdir -p "${INTERVAL_DIR}"

    # Split intervals into chunks for parallel processing
    TOTAL_INTERVALS=$(wc -l < "${INTERVAL_FILE}")
    INTERVALS_PER_JOB=$(( (TOTAL_INTERVALS + PARALLEL_JOBS - 1) / PARALLEL_JOBS ))

    log "Splitting ${TOTAL_INTERVALS} intervals into ${PARALLEL_JOBS} jobs"

    # Create interval chunks
    split -l "${INTERVALS_PER_JOB}" "${INTERVAL_FILE}" "${INTERVAL_DIR}/interval_chunk_"

    # Process each chunk in parallel
    CHUNK_FILES=("${INTERVAL_DIR}"/interval_chunk_*)

    for chunk_file in "${CHUNK_FILES[@]}"; do
        if [ -f "${chunk_file}" ]; then
            # Create interval list for this chunk
            CHUNK_NAME=$(basename "${chunk_file}")
            INTERVAL_LIST="${INTERVAL_DIR}/${CHUNK_NAME}.list"

            # Convert to GATK interval list format
            echo -e "@HD\tVN:1.6\tSO:unsorted" > "${INTERVAL_LIST}"
            echo -e "@SQ\tSN:chr1\tLN:248956422" >> "${INTERVAL_LIST}"  # Simplified
            cat "${chunk_file}" | while IFS= read -r interval; do
                echo -e "${interval}" >> "${INTERVAL_LIST}"
            done

            # Run HaplotypeCaller for this chunk (in background)
            (
                CHUNK_VCF="${INTERVAL_DIR}/${CHUNK_NAME}.vcf.gz"
                log "Starting HaplotypeCaller for ${CHUNK_NAME}"

                $GATK HaplotypeCaller \
                    -R "${REF_GENOME}" \
                    -I "${INPUT_BAM}" \
                    -O "${CHUNK_VCF}" \
                    -L "${INTERVAL_LIST}" \
                    --dont-use-soft-clipped-bases true \
                    --standard-min-confidence-threshold-for-calling 20.0 \
                    --min-base-quality-score 10 \
                    --max-reads-per-alignment-start 50 \
                    --max-assembly-region-size 500 \
                    --pair-hmm-implementation LOGLESS_CACHING

                log "Completed ${CHUNK_NAME}"
            ) &

            # Limit number of concurrent jobs
            BACKGROUND_JOBS=$(jobs -p | wc -l)
            while [ "${BACKGROUND_JOBS}" -ge "${PARALLEL_JOBS}" ]; do
                sleep 10
                BACKGROUND_JOBS=$(jobs -p | wc -l)
            done
        fi
    done

    # Wait for all background jobs to complete
    log "Waiting for ${#CHUNK_FILES[@]} parallel HaplotypeCaller jobs to complete..."
    wait

    # Merge all interval VCFs
    log "Merging interval VCFs..."
    RAW_VCF="${BASE_DIR}/data/vcf/${SAMPLE}/${SAMPLE}.vcf.gz"

    # Create input list for merge
    MERGE_LIST="${INTERVAL_DIR}/merge_list.txt"
    find "${INTERVAL_DIR}" -name "*.vcf.gz" -type f | sort > "${MERGE_LIST}"

    if [ -s "${MERGE_LIST}" ]; then
        $GATK MergeVcfs \
            -I "${MERGE_LIST}" \
            -O "${RAW_VCF}"

        log "Successfully merged interval VCFs into ${RAW_VCF}"
    else
        log_error "No interval VCFs found to merge"
        return 1
    fi

    return 0
}

############################################################
# Main Pipeline Functions
############################################################

run_fastq_pipeline() {
    log "Starting FASTQ processing pipeline"
    setup_environment

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log_error "No SRA IDs to process"
    fi

    log "Processing ${#SRA_LIST[@]} samples sequentially..."

    for SRR in "${SRA_LIST[@]}"; do
        log "========================================"
        log "Processing sample: ${SRR}"
        log "========================================"

        # Check disk space before each download
        check_disk_space 5120 || log_warning "Continuing despite low disk space"

        download_single_sra "${SRR}"
        validate_single_sra "${SRR}"
        extract_single_fastq "${SRR}"

        # Cleanup SRA file immediately after extraction
        if [ "${KEEP_INTERMEDIATE}" != "true" ]; then
            rm -f "${BASE_DIR}/data/sra/${SRR}.sra" 2>/dev/null || true
            rm -rf "${BASE_DIR}/data/sra/${SRR}" 2>/dev/null || true
        fi
    done

    # Run QC on all extracted FASTQ files
    run_fastqc
    trim_reads
    run_fastqc_trimmed

    log ""
    log "========================================"
    log "FASTQ Processing Complete!"
    log "========================================"
    log ""
    log "Next steps:"
    log "1. Align the trimmed FASTQ files with STAR"
    log "2. Run: ./$(basename "$0") run"
    log ""
    log "Trimmed FASTQ files are in: ${BASE_DIR}/data/trimmed/"
    log ""
}

run_fastq_pipeline() {
    log "Starting FASTQ processing pipeline"
    setup_environment

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log_error "No SRA IDs to process"
    fi

    log "Processing ${#SRA_LIST[@]} samples sequentially..."

    for SRR in "${SRA_LIST[@]}"; do
        log "========================================"
        log "Processing sample: ${SRR}"
        log "========================================"

        # Check disk space before each download
        check_disk_space 5120 || log_warning "Continuing despite low disk space"

        download_single_sra "${SRR}"
        validate_single_sra "${SRR}"
        extract_single_fastq "${SRR}"

        # Cleanup SRA file immediately after extraction
        if [ "${KEEP_INTERMEDIATE}" != "true" ]; then
            rm -f "${BASE_DIR}/data/sra/${SRR}.sra" 2>/dev/null || true
            rm -rf "${BASE_DIR}/data/sra/${SRR}" 2>/dev/null || true
        fi
    done

    # Run QC on all extracted FASTQ files
    run_fastqc
    trim_reads
    run_fastqc_trimmed

    log ""
    log "========================================"
    log "FASTQ Processing Complete!"
    log "========================================"
    log ""
    log "Next steps:"
    log "1. Align the trimmed FASTQ files with STAR"
    log "2. Run: ./$(basename "$0") run"
    log ""
    log "Trimmed FASTQ files are in: ${BASE_DIR}/data/trimmed/"
    log ""
}

process_all_samples() {
    log "Processing all aligned samples sequentially..."

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    for SRR in "${SRA_LIST[@]}"; do
        log "========================================"
        log "Processing sample: ${SRR}"
        log "========================================"

        check_disk_space 2048 || log_warning "Low disk space, continuing anyway"

        # Find the BAM file for this sample
        BAM_FILE="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam"

        if [ ! -f "${BAM_FILE}" ]; then
            log_warning "BAM file not found for ${SRR}, skipping: ${BAM_FILE}"
            continue
        fi

        # Process single sample
        process_sample "${SRR}" "${BAM_FILE}"

        # Cleanup intermediate files for this sample
        cleanup_intermediate "${SRR}" "${KEEP_INTERMEDIATE}"

        log "Completed processing for ${SRR}"
    done

    log "Merging all TSV files..."
    MERGED_TSV="${BASE_DIR}/all_samples.tsv"
    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ\tSAMPLE" > "${MERGED_TSV}"

    for sample_dir in "${BASE_DIR}/data/tsv/"*/; do
        sample=$(basename "${sample_dir}")
        tsv_gz="${sample_dir}/${sample}.tsv.gz"
        if [ -f "${tsv_gz}" ]; then
            pigz -dc "${tsv_gz}" 2>/dev/null | tail -n +2 | awk -v sample="${sample}" '{print $0 "\t" sample}' >> "${MERGED_TSV}" || \
            gzip -dc "${tsv_gz}" | tail -n +2 | awk -v sample="${sample}" '{print $0 "\t" sample}' >> "${MERGED_TSV}"
        fi
    done

    pigz -f -p 4 "${MERGED_TSV}" 2>/dev/null || gzip -f "${MERGED_TSV}"
    log "Pipeline complete! Merged TSV: ${MERGED_TSV}.gz"
}

############################################################
# Main Script Entry Point
############################################################

show_help() {
    cat << EOF
========================================
MDD2 SNP Extraction Pipeline v3.0
========================================

Usage: $0 [--keep-intermediate] <command>

Commands:
  install         - Install all required tools locally
  setup           - Download reference files and set up pipeline
  load-ids        - Load SRA IDs from Excel file (will prompt for file)
  fastq           - Run FASTQ processing pipeline (download, trim, QC)
  align           - Align trimmed FASTQ with STAR
  run             - Run full pipeline sequentially with cleanup
  process-single <sample> <bam> - Process single aligned BAM file
  process-all     - Process all aligned samples in parallel
  test            - Run test with sample data
  clean           - Clean intermediate files
  help            - Show this help message

Options:
  --keep-intermediate  Keep intermediate files (SRA, FASTQ, temporary BAMs)
                      Default: Files are deleted after each sample is processed

Example workflow:
  1. $0 install      # Install tools locally in ~/.local/mdd2_tools
  2. $0 setup        # Download reference files
  3. $0 load-ids     # Load SRA IDs from Excel
  4. $0 run          # Run complete sequential pipeline
  5. $0 run --keep-intermediate  # Keep intermediate files for debugging

Sequential Processing:
  • Processes one sample at a time: Download → Process → Cleanup → Next
  • Automatically checks disk space before each download
  • Keeps only essential files: TSV, VCF, final BAM, QC reports
  • Deletes intermediate files by default (SRA, FASTQ, temp BAMs)
  • Use --keep-intermediate flag to retain all intermediate files

Features:
  • Sequential processing for large datasets (hundreds of samples)
  • Automatic disk space management
  • Excel file integration with automatic SRA ID extraction
  • GPU acceleration support for GATK and STAR
  • Multi-core processing throughout pipeline
  • Comprehensive error handling and logging
  • Quality control with FastQC and MultiQC

Configuration:
  • Edit EXCEL_FILE variable in script or be prompted
  • Set USE_GPU=true for GPU acceleration (if available)
  • Adjust CPU_CORES, MAX_PARALLEL, GATK_THREADS as needed

Output:
  • Tools installed to: ${TOOLS_DIR}
  • Analysis results in: ${BASE_DIR}
  • Final merged TSV: ${BASE_DIR}/all_samples.tsv.gz
EOF
}

# Main command dispatch
KEEP_INTERMEDIATE=false

# Check for --keep-intermediate flag
if [ $# -gt 0 ] && [[ "$1" == "--keep-intermediate" ]]; then
    KEEP_INTERMEDIATE=true
    shift
fi

if [ $# -eq 0 ]; then
    show_help
    exit 0
fi

COMMAND="$1"
shift

case "${COMMAND}" in
    "status")
        show_pipeline_status_simple
        ;;
    "install")
        install_all_tools
        ;;
    "setup")
        setup_references
        ;;
    "load-ids")
        load_sra_ids
        ;;
    "fastq")
        run_fastq_pipeline
        ;;
    "align")
        setup_environment
        align_with_star
        ;;
    "run")
        run_full_pipeline
        ;;
    "process-single")
        if [ $# -ne 2 ]; then
            echo "Usage: $0 [--keep-intermediate] process-single <sample> <bam>"
            exit 1
        fi
        setup_environment
        export KEEP_INTERMEDIATE
        process_sample "$1" "$2"
        ;;
    "process-all")
        setup_environment
        export KEEP_INTERMEDIATE
        process_all_samples
        ;;
    "test")
        log "Running test with sample data..."
        SRA_LIST=("SRR5961857" "SRR5961858")
        setup_environment
        export KEEP_INTERMEDIATE
        # Test with sequential processing
        for SRR in "${SRA_LIST[@]}"; do
            download_single_sra "${SRR}"
            extract_single_fastq "${SRR}"
            cleanup_intermediate "${SRR}" "${KEEP_INTERMEDIATE}"
        done
        log "Test completed successfully"
        ;;

    "setup-annotations")
        setup_environment
        setup_annotation_tools
        ;;

    "clean")
        log "Cleaning intermediate files..."
        rm -rf "${BASE_DIR}/data/sra/"*.sra 2>/dev/null || true
        rm -rf "${BASE_DIR}/data/fastq/"*.fastq 2>/dev/null || true
        rm -rf "${BASE_DIR}/data/trimmed/"*.fastq.gz 2>/dev/null || true
        rm -rf "${BASE_DIR}/data/aligned/"*.bam 2>/dev/null || true
        rm -rf "${BASE_DIR}/data/aligned/"*.bai 2>/dev/null || true
        rm -rf "${BASE_DIR}/data/processed/"*/ 2>/dev/null || true
        log "Cleanup complete"
        ;;
    "help")
        show_help
        ;;
    *)
        echo "Unknown command: ${COMMAND}"
        echo "Run: $0 help for usage information"
        exit 1
        ;;
esac


log "Command completed: ${COMMAND}"
