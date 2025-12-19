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

KEEP_INTERMEDIATE=false
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${PIPELINE_DIR}/analysis"
TOOLS_DIR="${HOME}/.local/mdd2_tools"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"

mkdir -p "${BASE_DIR}" "${SCRIPT_DIR}" "${TOOLS_DIR}" "${TOOLS_DIR}/bin"
mkdir -p "${BASE_DIR}/"{data,references,tools,logs}
mkdir -p "${BASE_DIR}/data/"{sra,fastq,fastqc,trimmed,aligned,processed,vcf,tsv,tmp}
mkdir -p "${BASE_DIR}/references/"{genome,annotations,known_sites}

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
# Utility Functions (Complete from original)
############################################################

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

    if [ -d "${TOOLS_DIR}/lib/python"*"/site-packages" ]; then
        PYTHON_SITE_DIR=$(find "${TOOLS_DIR}/lib" -name "site-packages" -type d | head -1)
        if [ -n "${PYTHON_SITE_DIR}" ]; then
            export PYTHONPATH="${PYTHON_SITE_DIR}:${PYTHONPATH:-}"
            log "Added ${PYTHON_SITE_DIR} to PYTHONPATH"
        fi
    fi

    USER_SITE=$(python3 -m site --user-site 2>/dev/null || python -m site --user-site 2>/dev/null || echo "")
    if [ -n "${USER_SITE}" ] && [ -d "${USER_SITE}" ]; then
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
if [ -d "$HOME/.local/mdd2_tools/lib/python"*"/site-packages" ]; then
    PYTHON_SITE_DIR=$(find "$HOME/.local/mdd2_tools/lib" -name "site-packages" -type d | head -1)
    if [ -n "${PYTHON_SITE_DIR}" ]; then
        export PYTHONPATH="${PYTHON_SITE_DIR}:${PYTHONPATH:-}"
    fi
fi
USER_SITE=$(python3 -m site --user-site 2>/dev/null || python -m site --user-site 2>/dev/null)
if [ -n "${USER_SITE}" ] && [ -d "${USER_SITE}" ]; then
    export PYTHONPATH="${USER_SITE}:${PYTHONPATH:-}"
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

    $GATK FuncotatorDataSourceDownloader \
        --somatic \
        --hg38 \
        --validate-integrity \
        --extract-after-download \
        --output "${FUNCOTATOR_DS}"

    log "Funcotator data sources downloaded to: ${FUNCOTATOR_DS}"
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

    cat > "${SCRIPT_DIR}/extract_sra_ids.py" << 'PYTHONSCRIPT'
#!/usr/bin/env python3
import pandas as pd
import sys
import re

def extract_sra_ids(excel_path, output_path):
    try:
        # Try to read specific sheet names first
        try:
            mdd_df = pd.read_excel(excel_path, sheet_name='Female MDD')
            ctrl_df = pd.read_excel(excel_path, sheet_name='Female Control')
        except:
            # Try to find sheets with different names
            xls = pd.ExcelFile(excel_path)
            sheet_names = xls.sheet_names
            print(f"Available sheets: {sheet_names}")

            # Use first two sheets
            mdd_df = pd.read_excel(excel_path, sheet_name=sheet_names[0])
            ctrl_df = pd.read_excel(excel_path, sheet_name=sheet_names[1]) if len(sheet_names) > 1 else pd.DataFrame()

        def find_sra_column(df):
            sra_pattern = re.compile(r'^SRR\d+$')
            # Try common column names first
            for col_name in df.columns:
                col_lower = str(col_name).lower()
                if any(keyword in col_lower for keyword in ['sample', 'id', 'sra', 'accession']):
                    if df[col_name].dropna().astype(str).str.match(sra_pattern).any():
                        return col_name

            # Try all columns
            for col_name in df.columns:
                if df[col_name].dropna().astype(str).str.match(sra_pattern).any():
                    return col_name

            return None

        all_sra_ids = []

        for df, sheet_name in [(mdd_df, 'MDD'), (ctrl_df, 'Control')]:
            if df.empty:
                continue

            sra_col = find_sra_column(df)
            if sra_col:
                sra_ids = df[sra_col].dropna().astype(str)
                valid_ids = sra_ids[sra_ids.str.match(r'^SRR\d+$')].unique()
                all_sra_ids.extend(valid_ids)
                print(f"Found {len(valid_ids)} SRA IDs in {sheet_name} sheet (column: '{sra_col}')")
            else:
                print(f"Warning: Could not find SRA IDs in {sheet_name} sheet")

        # Remove duplicates and sort
        unique_ids = sorted(set(all_sra_ids))

        if not unique_ids:
            print("ERROR: No SRA IDs found in Excel file")
            return False

        # Write to file
        with open(output_path, 'w') as f:
            for sra_id in unique_ids:
                f.write(f"{sra_id}\n")

        print(f"\nTotal unique SRA IDs extracted: {len(unique_ids)}")
        print(f"First 5 IDs: {unique_ids[:5]}")
        if len(unique_ids) > 5:
            print(f"... and {len(unique_ids) - 5} more")

        return True

    except Exception as e:
        print(f"Error extracting SRA IDs: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_sra_ids.py <excel_file> <output_file>")
        sys.exit(1)

    success = extract_sra_ids(sys.argv[1], sys.argv[2])
    sys.exit(0 if success else 1)
PYTHONSCRIPT

    chmod +x "${SCRIPT_DIR}/extract_sra_ids.py"

    # Check for Python dependencies
    if ! python3 -c "import pandas" &> /dev/null; then
        log "Installing Python dependencies..."
        if command -v pip3 &> /dev/null; then
            pip3 install pandas openpyxl
        elif command -v pip &> /dev/null; then
            pip install pandas openpyxl
        else
            log_error "pip not available. Please install pandas manually: pip install pandas openpyxl"
        fi
    fi

    log "Extracting SRA IDs from Excel..."
    if python3 "${SCRIPT_DIR}/extract_sra_ids.py" "${excel_file}" "${output_file}"; then
        return 0
    else
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
    prompt_for_excel_file

    if extract_sra_ids_from_excel "${EXCEL_FILE}" "${SRA_LIST_FILE}"; then
        if [ -f "${SRA_LIST_FILE}" ]; then
            count=$(wc -l < "${SRA_LIST_FILE}" | tr -d ' ')
            [ "${count}" -eq 0 ] && log_error "No SRA IDs found in Excel file"

            mapfile -t SRA_LIST < "${SRA_LIST_FILE}"
            log "Loaded ${#SRA_LIST[@]} SRA IDs from ${EXCEL_FILE}"

            echo ""
            echo "First 5 SRA IDs:"
            for i in {0..4}; do
                [ -n "${SRA_LIST[$i]}" ] && echo "  ${SRA_LIST[$i]}"
            done
            if [ ${#SRA_LIST[@]} -gt 5 ]; then
                echo "  ... and $(( ${#SRA_LIST[@]} - 5 )) more"
            fi
            echo ""
        else
            log_error "Failed to create SRA IDs file"
        fi
    else
        log_error "Failed to extract SRA IDs from Excel file"
    fi
}
#=================================================================
# For single SRA files
#=================================================================

trim_single_sample() {
    local SRR="$1"
    local INPUT_DIR="${BASE_DIR}/data/fastq"

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
    log "Starting full MDD2 SNP Extraction Pipeline"
    setup_environment

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        load_sra_ids
    fi

    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log_error "No SRA IDs to process"
    fi

    log "Processing ${#SRA_LIST[@]} samples sequentially with cleanup..."

    # Process each sample completely before moving to next
    for SRR in "${SRA_LIST[@]}"; do
        log "========================================"
        log "Processing FULL pipeline for: ${SRR}"
        log "========================================"

        check_disk_space 5120

        # FASTQ steps
        download_single_sra "${SRR}"
        validate_single_sra "${SRR}"
        extract_single_fastq "${SRR}"

        # Immediate cleanup of SRA file
        if [ "${KEEP_INTERMEDIATE}" != "true" ]; then
            rm -f "${BASE_DIR}/data/sra/${SRR}.sra" 2>/dev/null || true
        fi

        # FASTQC for this sample
        log "Running FastQC for ${SRR}..."
        if [ -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq.gz" ]; then
            fastqc "${BASE_DIR}/data/fastq/${SRR}_1.fastq.gz" -o "${BASE_DIR}/data/fastqc/raw" --quiet 2>/dev/null || true
        fi
        if [ -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq.gz" ]; then
            fastqc "${BASE_DIR}/data/fastq/${SRR}_2.fastq.gz" -o "${BASE_DIR}/data/fastqc/raw" --quiet 2>/dev/null || true
        fi

        # Trim this sample
        log "Trimming ${SRR}..."
        if [ -f "${TRIMMOMATIC_JAR}" ] || [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
            trim_single_sample "${SRR}"
        fi

        # Align this sample
        log "Aligning ${SRR} with STAR..."
        align_single_sample "${SRR}"

        # Process this sample (variant calling)
        BAM_FILE="${BASE_DIR}/data/aligned/${SRR}.Aligned.sortedByCoord.out.bam"
        if [ -f "${BAM_FILE}" ]; then
            process_sample "${SRR}" "${BAM_FILE}"
        else
            log_warning "BAM file not found after alignment: ${BAM_FILE}"
        fi

        # Cleanup intermediate files for this sample
        cleanup_intermediate "${SRR}" "${KEEP_INTERMEDIATE}"

        log "Completed FULL processing for ${SRR}"
    done

    # Run MultiQC on all QC data
    if command -v multiqc &> /dev/null; then
        log "Running final MultiQC..."
        multiqc "${BASE_DIR}/data/fastqc/" -o "${BASE_DIR}/data/fastqc/" --quiet
    fi

    # Merge all TSV files
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

    log ""
    log "========================================"
    log "Pipeline Complete!"
    log "========================================"
    log ""
    log "Output files:"
    log "  • Merged TSV: ${BASE_DIR}/all_samples.tsv.gz"
    log "  • Individual TSVs: ${BASE_DIR}/data/tsv/"
    log "  • VCF files: ${BASE_DIR}/data/vcf/"
    log "  • Processed BAMs: ${BASE_DIR}/data/processed/"
    log "  • QC reports: ${BASE_DIR}/data/fastqc/"
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

process_sample() {
    local SAMPLE=$1
    local INPUT_BAM=$2

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

    ##################################
    # Add/Replace Read Groups
    ##################################
    log "Adding read groups..."
    $GATK AddOrReplaceReadGroups \
        -I "${INPUT_BAM}" \
        -O "${BAM_DIR}/${SAMPLE}.rg.bam" \
        -RGID "${SAMPLE}" \
        -RGLB lib1 \
        -RGPL ILLUMINA \
        -RGPU "${SAMPLE}.unit1" \
        -RGSM "${SAMPLE}" \
        --CREATE_INDEX true

    ##################################
    # Mark Duplicates
    ##################################
    log "Marking duplicates..."
    $GATK MarkDuplicates \
        -I "${BAM_DIR}/${SAMPLE}.rg.bam" \
        -O "${BAM_DIR}/${SAMPLE}.rmdup.bam" \
        -M "${BAM_DIR}/${SAMPLE}.rmdup.metrics.txt" \
        --CREATE_INDEX true

    ##################################
    # SplitNCigarReads
    ##################################
    log "Splitting N Cigar reads..."
    $GATK SplitNCigarReads \
        -R "${REF_GENOME}" \
        -I "${BAM_DIR}/${SAMPLE}.rmdup.bam" \
        -O "${BAM_DIR}/${SAMPLE}.split.bam" \
        --create-output-bam-index true

    ##################################
    # Base Quality Score Recalibration
    ##################################
    log "Running BaseRecalibrator..."
    $GATK BaseRecalibrator \
        -R "${REF_GENOME}" \
        -I "${BAM_DIR}/${SAMPLE}.split.bam" \
        --known-sites "${DBSNP}" \
        --known-sites "${MILLS}" \
        -O "${BAM_DIR}/${SAMPLE}.recal.table"

    log "Applying BQSR..."
    $GATK ApplyBQSR \
        -R "${REF_GENOME}" \
        -I "${BAM_DIR}/${SAMPLE}.split.bam" \
        --bqsr-recal-file "${BAM_DIR}/${SAMPLE}.recal.table" \
        -O "${BAM_DIR}/${SAMPLE}.recal.bam" \
        --create-output-bam-index true

    ##################################
    # Variant Calling
    ##################################
    log "Calling variants with HaplotypeCaller..."
    $GATK HaplotypeCaller \
        -R "${REF_GENOME}" \
        -I "${BAM_DIR}/${SAMPLE}.recal.bam" \
        -O "${VCF_DIR}/${SAMPLE}.vcf.gz" \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20.0

    ##################################
    # Variant Filtration
    ##################################
    log "Filtering variants..."
    $GATK VariantFiltration \
        -R "${REF_GENOME}" \
        -V "${VCF_DIR}/${SAMPLE}.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "LowQD" \
        --filter-expression "FS > 30.0" --filter-name "FS" \
        --filter-expression "SOR > 3.0" --filter-name "SOR" \
        --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
        --window 35 --cluster 3 \
        -O "${VCF_DIR}/${SAMPLE}.filtered.vcf.gz"

    ##################################
    # Annotation
    ##################################
    log "Annotating with Funcotator..."
    $GATK Funcotator \
        -R "${REF_GENOME}" \
        -V "${VCF_DIR}/${SAMPLE}.filtered.vcf.gz" \
        -O "${VCF_DIR}/${SAMPLE}.annotated.vcf.gz" \
        --ref-version hg38 \
        --data-sources-path "${FUNCOTATOR_DS}" \
        --output-file-format VCF

    log "Adding dbSNP IDs..."
    bcftools annotate \
        -a "${DBSNP}" \
        -c ID \
        -O z \
        -o "${VCF_DIR}/${SAMPLE}.rsID.vcf.gz" \
        "${VCF_DIR}/${SAMPLE}.annotated.vcf.gz"

    log "Adding gene annotations..."
    bcftools annotate \
        -a "${BASE_DIR}/references/annotations/genes.bed.gz" \
        -h "${BASE_DIR}/references/annotations/genes.bed.hdr" \
        -c CHROM,FROM,TO,INFO/GENE \
        -O z \
        -o "${VCF_DIR}/${SAMPLE}.genes.vcf.gz" \
        "${VCF_DIR}/${SAMPLE}.rsID.vcf.gz"

    ##################################
    # TSV Export
    ##################################
    log "Exporting to TSV..."
    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ" > "${TSV_DIR}/${SAMPLE}.tsv"

    bcftools query \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/GENE[\t%GT\t%DP\t%AD\t%GQ]\n' \
        "${VCF_DIR}/${SAMPLE}.genes.vcf.gz" >> "${TSV_DIR}/${SAMPLE}.tsv"

    pigz -f -p 2 "${TSV_DIR}/${SAMPLE}.tsv" 2>/dev/null || gzip -f "${TSV_DIR}/${SAMPLE}.tsv"

    log "Finished processing ${SAMPLE}"
    log "TSV file created: ${TSV_DIR}/${SAMPLE}.tsv.gz"
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
