#!/bin/bash

##################################
# MDD2 SNP Extraction Pipeline
# Version: 1.0
# Description: Complete pipeline for downloading FASTQ files and
#              extracting SNP TSV files using GATK best practices
##################################

set -e  # Exit on error
set -u  # Treat unset variables as error

##################################
# Configuration Section
##################################

# Base directory setup
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${PIPELINE_DIR}/analysis"
TOOLS_DIR="${HOME}/.local/mdd2_tools"  # Install tools in user's home directory
SCRIPT_DIR="${PIPELINE_DIR}/scripts"

# Create directories if they don't exist
mkdir -p "${BASE_DIR}" "${SCRIPT_DIR}" "${TOOLS_DIR}"

# Directory structure
mkdir -p "${BASE_DIR}/"{data,references,tools,logs}
mkdir -p "${BASE_DIR}/data/"{sra,fastq,fastqc,trimmed,aligned,processed,vcf,tsv}
mkdir -p "${BASE_DIR}/references/"{genome,annotations,known_sites}

# Reference files
REF_GENOME="${BASE_DIR}/references/genome/GRCh38.primary_assembly.genome.fa"
DBSNP="${BASE_DIR}/references/known_sites/dbsnp_146.hg38.vcf.gz"
MILLS="${BASE_DIR}/references/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR_DS="${BASE_DIR}/references/annotations/funcotator_dataSources"

# Sample list (to be provided by user - EDIT THIS SECTION)
# Add your SRA IDs here
SRA_LIST=("SRR5961857" "SRR5961858")  # EXAMPLE - REPLACE WITH YOUR SAMPLES

##################################
# Function Definitions
##################################

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${BASE_DIR}/logs/pipeline.log"
}

check_tool() {
    local tool="$1"
    if command -v "${tool}" &> /dev/null; then
        log "✓ ${tool} is available"
        return 0
    fi

    # Check in tools directory
    if [ -x "${TOOLS_DIR}/bin/${tool}" ]; then
        log "✓ ${tool} is available in tools directory"
        return 0
    fi

    # Check for specific tool paths
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
            ;;
        "trimmomatic")
            if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
                log "✓ trimmomatic is available in tools directory"
                return 0
            fi
            ;;
        "prefetch"|"fasterq-dump"|"vdb-validate")
            if [ -x "${TOOLS_DIR}/sratoolkit/bin/${tool}" ]; then
                log "✓ ${tool} is available in tools directory"
                return 0
            fi
            ;;
    esac

    log "✗ ${tool} is not installed or not in PATH"
    return 1
}

check_file() {
    if [ ! -f "$1" ]; then
        log "ERROR: File not found: $1"
        return 1
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

setup_environment() {
    log "Setting up environment..."

    # Add tool directories to PATH
    add_to_path "${TOOLS_DIR}/bin"
    add_to_path "${TOOLS_DIR}/FastQC"
    add_to_path "${TOOLS_DIR}/sratoolkit/bin"

    # Check for GATK in multiple locations
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

    # Set Java path for Trimmomatic
    if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        export TRIMMOMATIC_JAR="${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"
    fi

    # Create environment file
    cat > "${BASE_DIR}/env.sh" << EOF
#!/bin/bash
# MDD2 Pipeline Environment Setup
export PATH="${TOOLS_DIR}/bin:${TOOLS_DIR}/FastQC:${TOOLS_DIR}/sratoolkit/bin:\${PATH}"
export PATH="${TOOLS_DIR}/gatk-4.6.2.0:\${PATH}"
export PATH="${TOOLS_DIR}/gatk-4.4.0.0:\${PATH}"
if [ -x "${TOOLS_DIR}/gatk" ]; then
    export GATK="${TOOLS_DIR}/gatk"
elif [ -x "${TOOLS_DIR}/gatk-4.6.2.0/gatk" ]; then
    export GATK="${TOOLS_DIR}/gatk-4.6.2.0/gatk"
elif [ -x "${TOOLS_DIR}/gatk-4.4.0.0/gatk" ]; then
    export GATK="${TOOLS_DIR}/gatk-4.4.0.0/gatk"
fi
if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
    export TRIMMOMATIC_JAR="${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar"
fi
EOF

    chmod +x "${BASE_DIR}/env.sh"
    log "Environment setup complete. Source with: source ${BASE_DIR}/env.sh"
}

download_reference() {
    log "Downloading reference genome..."
    mkdir -p "$(dirname "${REF_GENOME}")"

    wget -q -O "${REF_GENOME}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"

    if [ -f "${REF_GENOME}.gz" ]; then
        gunzip "${REF_GENOME}.gz"

        log "Indexing reference genome..."
        if check_tool "samtools"; then
            samtools faidx "${REF_GENOME}"
        else
            log "WARNING: samtools not available, skipping faidx"
        fi

        # Check if GATK is available for CreateSequenceDictionary
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
    # dbSNP
    mkdir -p "$(dirname "${DBSNP}")"
    wget -q -O "${DBSNP}" \
        "https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/dbsnp_146.hg38.vcf.gz"
    wget -q -O "${DBSNP}.tbi" \
        "https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/dbsnp_146.hg38.vcf.gz.tbi"

    # Mills INDELs
    wget -q -O "${MILLS}" \
        "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    wget -q -O "${MILLS}.tbi" \
        "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

    return 0
}

generate_gene_bed() {
    log "Generating gene BED file from Gencode annotation..."
    GTF="${BASE_DIR}/references/annotations/gencode.v49.annotation.gtf.gz"
    GENE_BED="${BASE_DIR}/references/annotations/genes.bed"

    mkdir -p "$(dirname "${GTF}")"

    wget -q -O "${GTF}" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"

    gunzip -c "${GTF}" | \
    awk 'BEGIN {OFS="\t"} $3 == "gene" {
        gene_name = "."
        if (match($0, /gene_name "([^"]+)"/, a)) gene_name = a[1]
        print $1, $4 - 1, $5, gene_name
    }' > "${GENE_BED}"

    if check_tool "bgzip"; then
        bgzip -c "${GENE_BED}" > "${GENE_BED}.gz"
    else
        gzip -c "${GENE_BED}" > "${GENE_BED}.gz"
    fi

    if check_tool "tabix"; then
        tabix -p bed "${GENE_BED}.gz"
    fi

    # Create header file for bcftools
    echo -e '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">' > "${GENE_BED}.hdr"

    log "Gene BED file created: ${GENE_BED}.gz"
}

download_funcotator() {
    log "Downloading Funcotator data sources..."

    # Check if GATK is available
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

##################################
# Tool Installation Functions
##################################

install_sratoolkit() {
    log "Installing SRA Toolkit..."

    if [ -d "${TOOLS_DIR}/sratoolkit" ] && [ -x "${TOOLS_DIR}/sratoolkit/bin/prefetch" ]; then
        log "SRA Toolkit already installed"
        return 0
    fi

    cd "${TOOLS_DIR}"

    # Try different versions for different architectures
    ARCH=$(uname -m)

    if [ "${ARCH}" == "x86_64" ]; then
        # Try Ubuntu version first (common on HPC)
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
            # Make all binaries executable
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

    # Check if already installed
    if [ -x "${TOOLS_DIR}/FastQC/fastqc" ]; then
        log "FastQC already installed"
        return 0
    fi

    mkdir -p "${TOOLS_DIR}/FastQC"
    cd "${TOOLS_DIR}"

    # Remove any existing zip file
    rm -f fastqc.zip

    # Download FastQC
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

    # Extract
    log "Extracting FastQC..."
    unzip -q fastqc.zip 2>/dev/null || {
        log "ERROR: Failed to extract FastQC zip file"
        return 1
    }

    # Look for the fastqc file in common extraction patterns
    FASTQC_FOUND=false

    # Pattern 1: Standard FastQC directory
    if [ -f "FastQC/fastqc" ]; then
        mv FastQC/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
        FASTQC_FOUND=true
    # Pattern 2: Versioned directory (v0.12.1)
    elif [ -f "fastqc_v0.12.1/fastqc" ]; then
        mv fastqc_v0.12.1/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
        FASTQC_FOUND=true
    # Pattern 3: Versioned directory (v0.11.9)
    elif [ -f "fastqc_v0.11.9/fastqc" ]; then
        mv fastqc_v0.11.9/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
        FASTQC_FOUND=true
    else
        # Search for any fastqc file
        FOUND_FILE=$(find . -type f -name "fastqc" 2>/dev/null | head -1)
        if [ -n "${FOUND_FILE}" ]; then
            # Get its directory
            FOUND_DIR=$(dirname "${FOUND_FILE}")
            mv "${FOUND_DIR}"/* "${TOOLS_DIR}/FastQC/" 2>/dev/null || true
            FASTQC_FOUND=true
        fi
    fi

    if ${FASTQC_FOUND}; then
        # Make sure the fastqc file is executable
        if [ -f "${TOOLS_DIR}/FastQC/fastqc" ]; then
            chmod +x "${TOOLS_DIR}/FastQC/fastqc"
            log "FastQC successfully installed to: ${TOOLS_DIR}/FastQC"
        else
            # The file might be inside a subdirectory
            SUB_FILE=$(find "${TOOLS_DIR}/FastQC" -type f -name "fastqc" 2>/dev/null | head -1)
            if [ -n "${SUB_FILE}" ]; then
                # Move it to the main directory
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

        # List extracted files for debugging
        log "Extracted files:"
        find . -type f -name "*" | head -20

        return 1
    fi

    # Clean up
    rm -f fastqc.zip
    rm -rf FastQC fastqc_v0.* 2>/dev/null || true

    return 0
}

install_multiqc() {
    log "Installing MultiQC..."

    # Check if already installed
    if command -v multiqc &> /dev/null; then
        log "MultiQC already available in PATH"
        return 0
    fi

    # Try different installation methods
    if command -v pip3 &> /dev/null; then
        pip3 install --prefix="${TOOLS_DIR}" multiqc 2>/dev/null && return 0
    fi

    if command -v pip &> /dev/null; then
        pip install --prefix="${TOOLS_DIR}" multiqc 2>/dev/null && return 0
    fi

    if command -v python3 &> /dev/null; then
        python3 -m pip install --prefix="${TOOLS_DIR}" multiqc 2>/dev/null && return 0
    fi

    # Try conda if available
    if command -v conda &> /dev/null; then
        conda install -c bioconda multiqc -y 2>/dev/null && return 0
    fi

    log "WARNING: Could not install MultiQC automatically"
    log "You may need to install it manually: pip install multiqc"
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

    # Try multiple download URLs
    wget -q --tries=3 --timeout=30 -O Trimmomatic-0.39.zip \
        "https://github.com/usadellab/Trimmomatic/files/5854849/Trimmomatic-0.39.zip" || \
    wget -q --tries=3 --timeout=30 -O Trimmomatic-0.39.zip \
        "https://github.com/usadellab/Trimmomatic/releases/download/v0.39/Trimmomatic-0.39.zip" || \
    wget -q --tries=3 --timeout=30 -O Trimmomatic-0.39.zip \
        "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"

    if [ -f "Trimmomatic-0.39.zip" ]; then
        unzip -q Trimmomatic-0.39.zip -d "${TOOLS_DIR}"

        # Check what was extracted
        if [ -d "Trimmomatic-0.39" ]; then
            # Already in correct location
            true
        elif [ -d "${TOOLS_DIR}/Trimmomatic-0.39" ]; then
            # Already moved
            true
        else
            # Try to find the extracted directory
            TRIMM_DIR=$(find "${TOOLS_DIR}" -name "*Trimmomatic*" -type d | head -1)
            if [ -n "${TRIMM_DIR}" ]; then
                mv "${TRIMM_DIR}" "${TOOLS_DIR}/Trimmomatic-0.39" 2>/dev/null || true
            fi
        fi

        rm -f Trimmomatic-0.39.zip

        # Verify installation
        if [ -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
            log "Trimmomatic installed to: ${TOOLS_DIR}/Trimmomatic-0.39"
        else
            # Try to find the jar file
            JAR_FILE=$(find "${TOOLS_DIR}/Trimmomatic-0.39" -name "*.jar" | head -1)
            if [ -n "${JAR_FILE}" ]; then
                # Rename to standard name
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

    # Check if already installed
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

    # Try multiple GATK versions and URLs
    GATK_DOWNLOADED=false

    # Try GATK 4.6.2.0 first
    if ! ${GATK_DOWNLOADED}; then
        wget -q --tries=3 --timeout=60 -O gatk-4.6.2.0.zip \
            "https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip" && \
        GATK_DOWNLOADED=true
    fi

    # Try GATK 4.4.0.0
    if ! ${GATK_DOWNLOADED}; then
        wget -q --tries=3 --timeout=60 -O gatk-4.4.0.0.zip \
            "https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip" && \
        GATK_DOWNLOADED=true
    fi

    # Try generic gatk download
    if ! ${GATK_DOWNLOADED}; then
        wget -q --tries=3 --timeout=60 -O gatk.zip \
            "https://github.com/broadinstitute/gatk/releases/latest/download/gatk.zip" && \
        GATK_DOWNLOADED=true
    fi

    if ! ${GATK_DOWNLOADED}; then
        log "ERROR: Failed to download GATK"
        return 1
    fi

    # Extract the downloaded file
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

    # Find and organize the GATK files
    log "Organizing GATK files..."

    # Look for the gatk executable
    GATK_BIN=$(find "${TOOLS_DIR}" -name "gatk" -type f | head -1)

    if [ -n "${GATK_BIN}" ]; then
        # Make it executable
        chmod +x "${GATK_BIN}"

        # Check if it's in a versioned directory
        GATK_DIR=$(dirname "${GATK_BIN}")
        GATK_BASE_DIR=$(dirname "${GATK_DIR}")

        if [[ "${GATK_DIR}" == *"gatk-4."* ]] && [ "${GATK_BASE_DIR}" == "${TOOLS_DIR}" ]; then
            # Already in a versioned directory under TOOLS_DIR
            log "GATK found in: ${GATK_DIR}"
        else
            # Move to a standard location
            if [[ "${GATK_DIR}" == *"gatk-4."* ]]; then
                # It's already in a versioned directory, move the whole directory
                mv "${GATK_DIR}" "${TOOLS_DIR}/" 2>/dev/null || true
            else
                # Create a version directory and move the binary
                mkdir -p "${TOOLS_DIR}/gatk-4.6.2.0"
                cp "${GATK_BIN}" "${TOOLS_DIR}/gatk-4.6.2.0/gatk"
                chmod +x "${TOOLS_DIR}/gatk-4.6.2.0/gatk"
            fi
        fi
    else
        # Try to find any Java jar file that might be GATK
        GATK_JAR=$(find "${TOOLS_DIR}" -name "*.jar" -type f | grep -i gatk | head -1)

        if [ -n "${GATK_JAR}" ]; then
            log "Found GATK JAR file: ${GATK_JAR}"
            # Create a wrapper script
            mkdir -p "${TOOLS_DIR}/gatk-4.6.2.0"
            cat > "${TOOLS_DIR}/gatk-4.6.2.0/gatk" << 'EOF'
#!/bin/bash
# GATK wrapper script
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
java -jar "${DIR}/../$(basename "$(find "${DIR}/.." -name "*.jar" -type f | grep -i gatk | head -1)")" "$@"
EOF
            chmod +x "${TOOLS_DIR}/gatk-4.6.2.0/gatk"

            # Also copy the jar file near the wrapper
            cp "${GATK_JAR}" "${TOOLS_DIR}/" 2>/dev/null || true
        else
            log "ERROR: Could not find GATK binary or JAR file after extraction"
            log "Contents of ${TOOLS_DIR}:"
            ls -la "${TOOLS_DIR}" 2>/dev/null || true
            return 1
        fi
    fi

    # Verify installation
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

    # Check if already installed
    if command -v samtools &> /dev/null && command -v bcftools &> /dev/null; then
        log "samtools and bcftools already available in PATH"
        return 0
    fi

    # Try to install via conda first (if available)
    if command -v conda &> /dev/null; then
        log "Installing via conda..."
        conda install -c bioconda samtools bcftools htslib -y 2>/dev/null && return 0
    fi

    # Create build directory
    BUILD_DIR="${TOOLS_DIR}/build"
    mkdir -p "${BUILD_DIR}"

    # Check for required build tools
    if ! command -v gcc &> /dev/null || ! command -v make &> /dev/null; then
        log "WARNING: gcc or make not found. Cannot compile tools from source."
        log "Please install samtools and bcftools manually."
        return 1
    fi

    # Install htslib (dependency for samtools and bcftools)
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

    # Install samtools
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

    # Install bcftools
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

    # Clean up
    rm -rf "${BUILD_DIR}" 2>/dev/null || true

    # Verify installation
    if [ -x "${TOOLS_DIR}/bin/samtools" ] && [ -x "${TOOLS_DIR}/bin/bcftools" ]; then
        log "Bioinformatics tools installed successfully"
    else
        log "WARNING: Some tools may not have installed correctly"
        log "You may need to install samtools and bcftools manually"
    fi
}

install_java() {
    log "Checking for Java..."

    if command -v java &> /dev/null; then
        JAVA_VERSION=$(java -version 2>&1 | head -1 | cut -d'"' -f2)
        log "Java is already available: ${JAVA_VERSION}"
        return 0
    fi

    log "Java not found in PATH. Checking if we can install locally..."

    # Try to find Java in common locations
    if [ -d "/usr/lib/jvm" ]; then
        JAVA_HOME=$(find /usr/lib/jvm -name "java*" -type d | grep -v debug | grep -v src | head -1)
        if [ -n "${JAVA_HOME}" ] && [ -x "${JAVA_HOME}/bin/java" ]; then
            export JAVA_HOME
            export PATH="${JAVA_HOME}/bin:${PATH}"
            log "Found Java at: ${JAVA_HOME}"
            return 0
        fi
    fi

    # Check if Java is in /usr/bin
    if [ -x "/usr/bin/java" ]; then
        export JAVA_HOME="/usr"
        log "Found Java at: /usr/bin/java"
        return 0
    fi

    log "WARNING: Java not found. Some tools (Trimmomatic, GATK) require Java."
    log "Please install Java manually or ask your HPC administrator."
    log "On Rocky Linux, you can try: module load java"
    return 1
}

##################################
# Pipeline Steps
##################################

download_sra_files() {
    log "Step 1: Downloading SRA files..."
    mkdir -p "${BASE_DIR}/data/sra/"

    for SRR in "${SRA_LIST[@]}"; do
        log "Downloading ${SRR}..."
        if command -v prefetch &> /dev/null; then
            prefetch --progress "${SRR}" -O "${BASE_DIR}/data/sra/"
        elif [ -x "${TOOLS_DIR}/sratoolkit/bin/prefetch" ]; then
            "${TOOLS_DIR}/sratoolkit/bin/prefetch" --progress "${SRR}" -O "${BASE_DIR}/data/sra/"
        else
            log "ERROR: prefetch command not found. Install SRA Toolkit first."
            return 1
        fi
    done
}

validate_sra_files() {
    log "Step 2: Validating SRA files..."
    for SRR in "${SRA_LIST[@]}"; do
        log "Validating ${SRR}..."
        SRA_FILE="${BASE_DIR}/data/sra/${SRR}.sra"
        if [ -f "${SRA_FILE}" ]; then
            if command -v vdb-validate &> /dev/null; then
                vdb-validate "${SRA_FILE}"
            elif [ -x "${TOOLS_DIR}/sratoolkit/bin/vdb-validate" ]; then
                "${TOOLS_DIR}/sratoolkit/bin/vdb-validate" "${SRA_FILE}"
            else
                log "WARNING: vdb-validate not found, skipping validation"
            fi
        else
            log "WARNING: ${SRA_FILE} not found, skipping validation"
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

        SRA_FILE="${BASE_DIR}/data/sra/${SRR}.sra"
        if [ ! -f "${SRA_FILE}" ]; then
            log "WARNING: ${SRA_FILE} not found, skipping"
            COUNT=$((COUNT + 1))
            continue
        fi

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

        # Compress FASTQ files
        if [ -f "${BASE_DIR}/data/fastq/${SRR}_1.fastq" ]; then
            gzip "${BASE_DIR}/data/fastq/${SRR}_1.fastq"
        fi
        if [ -f "${BASE_DIR}/data/fastq/${SRR}_2.fastq" ]; then
            gzip "${BASE_DIR}/data/fastq/${SRR}_2.fastq"
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

    # Run MultiQC if available
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

    # Check for Trimmomatic
    if [ ! -f "${TRIMMOMATIC_JAR}" ] && [ ! -f "${TOOLS_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        log "WARNING: Trimmomatic not found"
        log "Skipping trimming step"
        return 0
    fi

    # Use the jar file
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

    # Get list of R1 files
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

        # Check if R2 exists
        if [ ! -f "${r2}" ]; then
            log "WARNING: ${r2} not found, skipping ${base}"
            continue
        fi

        # Output files
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

align_with_star() {
    log "Step 7: Aligning with STAR..."
    mkdir -p "${BASE_DIR}/data/aligned"

    # Note: STAR 2-pass alignment requires genome indexing first
    # This is a simplified version - users should customize based on their STAR setup
    log "STAR alignment requires manual setup."
    log "Please configure STAR based on your system."
    log "Reference: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"

    # Example command structure (commented out)
    # for SAMPLE in "${SRA_LIST[@]}"; do
    #     log "Aligning ${SAMPLE}..."
    #     STAR --runThreadN 12 \
    #          --genomeDir /path/to/STAR_index \
    #          --readFilesIn "${BASE_DIR}/data/trimmed/${SAMPLE}_1.paired.fastq.gz" "${BASE_DIR}/data/trimmed/${SAMPLE}_2.paired.fastq.gz" \
    #          --readFilesCommand zcat \
    #          --outFileNamePrefix "${BASE_DIR}/data/aligned/${SAMPLE}" \
    #          --outSAMtype BAM SortedByCoordinate
    # done

    log "After STAR alignment, BAM files will be in: ${BASE_DIR}/data/aligned/"
}

##################################
# GATK Processing Functions
##################################

process_sample() {
    local SAMPLE=$1
    local INPUT_BAM=$2

    log "Processing sample: ${SAMPLE}"

    BAM_DIR="${BASE_DIR}/data/processed/${SAMPLE}"
    VCF_DIR="${BASE_DIR}/data/vcf/${SAMPLE}"
    TSV_DIR="${BASE_DIR}/data/tsv/${SAMPLE}"

    mkdir -p "${BAM_DIR}" "${VCF_DIR}" "${TSV_DIR}"

    # Check if GATK is available
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
    # Create header
    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tGT\tDP\tAD\tGQ" > "${TSV_DIR}/${SAMPLE}.tsv"

    # Extract data
    bcftools query \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/GENE[\t%GT\t%DP\t%AD\t%GQ]\n' \
        "${VCF_DIR}/${SAMPLE}.genes.vcf.gz" >> "${TSV_DIR}/${SAMPLE}.tsv"

    log "Finished processing ${SAMPLE}"
    log "TSV file created: ${TSV_DIR}/${SAMPLE}.tsv"
}

##################################
# Main Script Functions
##################################

install_tools() {
    log "Starting tool installation..."
    log "Tools will be installed to: ${TOOLS_DIR}"

    # Setup environment
    setup_environment

    # Install Java (check first)
    install_java

    # Install tools
    install_sratoolkit
    install_fastqc
    install_multiqc
    install_trimmomatic
    install_gatk
    install_bioinformatics_tools

    log ""
    log "========================================"
    log "Installation Complete!"
    log "========================================"
    log ""
    log "To use the tools, run:"
    log "  source ${BASE_DIR}/env.sh"
    log ""
    log "Or add to your ~/.bashrc:"
    log "  source ${BASE_DIR}/env.sh"
    log ""
    log "Then run: ./$(basename "$0") setup"
}

setup_pipeline() {
    log "Setting up pipeline..."

    # Setup environment
    setup_environment

    # Check for required tools
    echo "Checking for required tools..."

    REQUIRED_TOOLS=("java" "samtools" "bcftools" "bgzip" "tabix" "wget")

    for tool in "${REQUIRED_TOOLS[@]}"; do
        check_tool "${tool}"
    done

    # Check for GATK
    check_tool "gatk"

    # Download references
    echo ""
    echo "Downloading reference files..."
    download_reference
    generate_gene_bed

    echo ""
    echo "Setup complete!"
    echo "Edit the SRA_LIST array in the script to add your SRA IDs"
    echo "Then run: ./$(basename "$0") run"
}

run_pipeline() {
    log "Starting MDD2 SNP Extraction Pipeline"

    # Setup environment
    setup_environment

    # Check if SRA list is provided
    if [ ${#SRA_LIST[@]} -eq 0 ]; then
        log "ERROR: No SRA IDs provided. Please edit the SRA_LIST array in the script."
        log "Example: SRA_LIST=(\"SRR5961857\" \"SRR5961858\")"
        exit 1
    fi

    # Check for required tools
    if ! check_tool "prefetch" && [ ! -x "${TOOLS_DIR}/sratoolkit/bin/prefetch" ]; then
        log "ERROR: prefetch not found. Install SRA Toolkit first."
        log "Run: ./$(basename "$0") install"
        exit 1
    fi

    # Execute pipeline steps
    download_sra_files
    validate_sra_files
    extract_fastq
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
    log "2. Run: ./$(basename "$0") process <sample_name> <aligned_bam>"
    log ""
    log "Trimmed FASTQ files are in: ${BASE_DIR}/data/trimmed/"
    log ""
}

process_aligned_sample() {
    if [ $# -ne 2 ]; then
        echo "Usage: $0 process <sample_name> <aligned_bam>"
        echo "Example: $0 process SRR5961857 /path/to/aligned.bam"
        exit 1
    fi

    local SAMPLE=$1
    local ALIGNED_BAM=$2

    # Setup environment
    setup_environment

    # Check files
    if [ ! -f "${ALIGNED_BAM}" ]; then
        log "ERROR: Aligned BAM file not found: ${ALIGNED_BAM}"
        exit 1
    fi

    if [ ! -f "${REF_GENOME}" ]; then
        log "ERROR: Reference genome not found. Run setup first:"
        log "./$(basename "$0") setup"
        exit 1
    fi

    # Process the sample
    process_sample "${SAMPLE}" "${ALIGNED_BAM}"
}

##################################
# Main Script Entry Point
##################################

# Display help if no arguments
if [ $# -eq 0 ]; then
    echo "========================================"
    echo "MDD2 SNP Extraction Pipeline"
    echo "========================================"
    echo ""
    echo "Usage: $0 [command]"
    echo ""
    echo "Commands:"
    echo "  install     - Install all required tools locally"
    echo "  setup       - Download reference files and set up pipeline"
    echo "  run         - Run FASTQ processing pipeline (download, trim, QC)"
    echo "  process <sample> <bam> - Process aligned BAM file through GATK"
    echo "  help        - Show this help message"
    echo ""
    echo "Example workflow:"
    echo "  1. $0 install      # Install tools locally in ~/.local/mdd2_tools"
    echo "  2. $0 setup        # Download reference files"
    echo "  3. Edit SRA_LIST in the script"
    echo "  4. $0 run          # Process FASTQ files"
    echo "  5. Align with STAR (manual step)"
    echo "  6. $0 process SRR5961857 /path/to/aligned.bam"
    echo ""
    exit 0
fi

# Parse command line arguments
COMMAND="$1"
shift

case "${COMMAND}" in
    "install")
        install_tools
        ;;
    "setup")
        setup_pipeline
        ;;
    "run")
        run_pipeline
        ;;
    "process")
        process_aligned_sample "$@"
        ;;
    "help")
        # Show help (already displayed above)
        exec "$0"
        ;;
    *)
        echo "Unknown command: ${COMMAND}"
        echo "Run: $0 help for usage information"
        exit 1
        ;;
esac

log "Command completed: ${COMMAND}"
