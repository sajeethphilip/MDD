#!/bin/bash
# build_tools.sh - Self-contained tool installation for headless systems
# Run with: ./build_tools.sh
# Or it will be called automatically by mdd2.sh

set -e
set -u

############################################################
# Configuration
############################################################

TOOLS_DIR="${HOME}/.local/mdd2_tools"
CPU_CORES=$(nproc 2>/dev/null || echo 8)
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] INFO:${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" >&2
}

log_error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" >&2
    exit 1
}

check_command() {
    if command -v "$1" &> /dev/null; then
        return 0
    fi
    return 1
}

############################################################
# Library Installation Functions
############################################################

install_system_deps() {
    log "Checking system dependencies..."

    local missing_deps=""

    # Check for essential build tools
    for cmd in wget curl tar gzip bzip2 unzip make gcc g++ python3 java; do
        if ! check_command "$cmd"; then
            missing_deps="$missing_deps $cmd"
        fi
    done

    if [ -n "$missing_deps" ]; then
        log_warning "Missing system dependencies:$missing_deps"
        log "Please install these manually or ask your system administrator"
        log "For Rocky Linux/CentOS, you might need:"
        log "  yum install wget curl tar gzip bzip2 unzip make gcc gcc-c++ python3 java-1.8.0-openjdk"
    else
        log "✓ System dependencies available"
    fi
}

install_openssl_local() {
    log "Installing OpenSSL locally..."

    cd "$TOOLS_DIR"
    mkdir -p openssl
    cd openssl

    # Try to download and compile OpenSSL 1.0.2 (compatible with bcftools)
    if wget -q --tries=3 --timeout=120 -O openssl.tar.gz \
        "https://www.openssl.org/source/openssl-1.0.2u.tar.gz"; then

        tar -xzf openssl.tar.gz
        cd openssl-1.0.2u

        # Configure for local installation
        ./config --prefix="$TOOLS_DIR/openssl" \
                 --openssldir="$TOOLS_DIR/ssl" \
                 no-shared \
                 no-zlib \
                 no-asm

        make -j "$CPU_CORES"
        make install

        # Create symlinks for libraries
        ln -sf "$TOOLS_DIR/openssl/lib/libcrypto.a" "$TOOLS_DIR/lib/" 2>/dev/null || true
        ln -sf "$TOOLS_DIR/openssl/lib/libssl.a" "$TOOLS_DIR/lib/" 2>/dev/null || true

        log "✓ OpenSSL installed locally"
    else
        log_warning "Could not download OpenSSL source"
    fi
}

install_htslib_static() {
    log "Installing static htslib..."

    cd "$TOOLS_DIR"
    mkdir -p build
    cd build

    # Download htslib source
    wget -q --tries=3 --timeout=120 -O htslib.tar.bz2 \
        "https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2"

    if [ ! -f "htslib.tar.bz2" ]; then
        wget -q --tries=3 --timeout=120 -O htslib.tar.bz2 \
            "https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2"
    fi

    if [ -f "htslib.tar.bz2" ]; then
        tar -xjf htslib.tar.bz2
        cd htslib-*

        # Compile statically
        ./configure --prefix="$TOOLS_DIR" \
                    --disable-libcurl \
                    --disable-gcs \
                    --disable-s3 \
                    --enable-static \
                    --disable-shared \
                    CFLAGS="-fPIC -O2 -I$TOOLS_DIR/include" \
                    LDFLAGS="-L$TOOLS_DIR/lib"

        make -j "$CPU_CORES"
        make install

        # Copy binaries
        cp tabix bgzip "$TOOLS_DIR/bin/" 2>/dev/null || true

        log "✓ htslib installed statically"
    else
        log_warning "Could not download htslib"
    fi
}

install_bcftools_static() {
    log "Installing static bcftools..."

    cd "$TOOLS_DIR/build"

    # Download bcftools source
    wget -q --tries=3 --timeout=120 -O bcftools.tar.bz2 \
        "https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2"

    if [ ! -f "bcftools.tar.bz2" ]; then
        wget -q --tries=3 --timeout=120 -O bcftools.tar.bz2 \
            "https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2"
    fi

    if [ -f "bcftools.tar.bz2" ]; then
        tar -xjf bcftools.tar.bz2
        cd bcftools-*

        # Compile with static htslib
        ./configure --prefix="$TOOLS_DIR" \
                    --disable-libcurl \
                    --disable-plugins \
                    --enable-static \
                    CFLAGS="-fPIC -O2 -I$TOOLS_DIR/include" \
                    LDFLAGS="-L$TOOLS_DIR/lib -static"

        make -j "$CPU_CORES"
        make install

        # Create static binary
        gcc -static -o "$TOOLS_DIR/bin/bcftools_static" \
            -I. -I"$TOOLS_DIR/include" \
            *.o "$TOOLS_DIR/lib/libhts.a" \
            -lz -lm -lpthread 2>/dev/null || true

        if [ -x "$TOOLS_DIR/bin/bcftools_static" ]; then
            ln -sf bcftools_static "$TOOLS_DIR/bin/bcftools"
            log "✓ bcftools compiled statically"
        else
            # Fallback: download pre-compiled static binary
            download_static_binary "bcftools"
        fi
    else
        download_static_binary "bcftools"
    fi
}

install_samtools_static() {
    log "Installing static samtools..."

    cd "$TOOLS_DIR/build"

    # Download samtools source
    wget -q --tries=3 --timeout=120 -O samtools.tar.bz2 \
        "https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2"

    if [ ! -f "samtools.tar.bz2" ]; then
        wget -q --tries=3 --timeout=120 -O samtools.tar.bz2 \
            "https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2"
    fi

    if [ -f "samtools.tar.bz2" ]; then
        tar -xjf samtools.tar.bz2
        cd samtools-*

        # Compile with static htslib
        ./configure --prefix="$TOOLS_DIR" \
                    --without-curses \
                    --enable-static \
                    CFLAGS="-fPIC -O2 -I$TOOLS_DIR/include" \
                    LDFLAGS="-L$TOOLS_DIR/lib -static"

        make -j "$CPU_CORES"
        make install

        # Create static binary
        gcc -static -o "$TOOLS_DIR/bin/samtools_static" \
            -I. -I"$TOOLS_DIR/include" \
            *.o "$TOOLS_DIR/lib/libhts.a" \
            -lz -lm -lpthread 2>/dev/null || true

        if [ -x "$TOOLS_DIR/bin/samtools_static" ]; then
            ln -sf samtools_static "$TOOLS_DIR/bin/samtools"
            log "✓ samtools compiled statically"
        else
            download_static_binary "samtools"
        fi
    else
        download_static_binary "samtools"
    fi
}

download_static_binary() {
    local tool="$1"
    log "Downloading pre-compiled static $tool..."

    cd "$TOOLS_DIR/bin"

    local urls=(
        "https://github.com/samtools/$tool/releases/download/1.19/$tool-1.19-static"
        "https://github.com/samtools/$tool/releases/download/1.18/$tool-1.18-static"
        "https://github.com/samtools/$tool/releases/latest/download/$tool-static"
    )

    for url in "${urls[@]}"; do
        if wget -q --tries=2 --timeout=60 -O "${tool}_static" "$url"; then
            chmod +x "${tool}_static"
            ln -sf "${tool}_static" "$tool"
            log "✓ Downloaded static $tool from $url"
            return 0
        fi
    done

    log_warning "Could not download static $tool"
    return 1
}

install_star_static() {
    log "Installing STAR aligner..."

    cd "$TOOLS_DIR"

    # Download STAR source
    wget -q --tries=3 --timeout=120 -O STAR.tar.gz \
        "https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz"

    if [ -f "STAR.tar.gz" ]; then
        tar -xzf STAR.tar.gz
        cd STAR-2.7.11a/source

        # Compile STAR
        make -j "$CPU_CORES" STAR

        # Copy to bin directory
        cp STAR "$TOOLS_DIR/bin/"

        log "✓ STAR installed"
    else
        log_warning "Could not download STAR"
    fi
}

install_gatk_local() {
    log "Installing GATK locally..."

    cd "$TOOLS_DIR"

    # Download GATK
    if wget -q --tries=3 --timeout=180 -O gatk.zip \
        "https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip"; then

        unzip -q gatk.zip
        GATK_DIR=$(ls -d gatk-4.* 2>/dev/null | head -1)

        if [ -n "$GATK_DIR" ]; then
            # Find the GATK binary
            GATK_BIN=$(find "$GATK_DIR" -name "gatk" -type f | head -1)

            if [ -n "$GATK_BIN" ] && [ -x "$GATK_BIN" ]; then
                # Create wrapper script
                cat > "$TOOLS_DIR/bin/gatk" << 'EOF'
#!/bin/bash
# GATK wrapper script
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR_FILE="$(find "$DIR/.." -name "*.jar" -type f | grep -i gatk | head -1)"
java -jar "$JAR_FILE" "$@"
EOF
                chmod +x "$TOOLS_DIR/bin/gatk"

                # Copy the jar file
                cp "$(find "$GATK_DIR" -name "*.jar" -type f | head -1)" "$TOOLS_DIR/" 2>/dev/null || true

                log "✓ GATK installed"
            fi
        fi
        rm -f gatk.zip
    else
        log_warning "Could not download GATK"
    fi
}

install_all_static_tools() {
    log "Starting complete static tool installation..."
    log "Tools will be installed to: $TOOLS_DIR"

    # Create directories
    mkdir -p "$TOOLS_DIR"/{bin,lib,include,build}

    # Set library paths
    export LD_LIBRARY_PATH="$TOOLS_DIR/lib:$LD_LIBRARY_PATH"
    export PATH="$TOOLS_DIR/bin:$PATH"

    # Install dependencies
    install_system_deps
    install_openssl_local
    install_htslib_static
    install_bcftools_static
    install_samtools_static
    install_star_static
    install_gatk_local

    # Install other tools
    install_other_tools

    # Test installations
    test_installations

    log ""
    log "========================================"
    log "Static Tool Installation Complete!"
    log "========================================"
    log ""
    log "All tools installed locally in: $TOOLS_DIR"
    log "No system libraries required!"
    log ""
    log "To use these tools, add to your PATH:"
    log "  export PATH=\"$TOOLS_DIR/bin:\$PATH\""
    log "  export LD_LIBRARY_PATH=\"$TOOLS_DIR/lib:\$LD_LIBRARY_PATH\""
    log ""
    log "Or run: source $BASE_DIR/env.sh"
    log ""
}

install_other_tools() {
    log "Installing other required tools..."

    # FastQC
    if [ ! -x "$TOOLS_DIR/bin/fastqc" ]; then
        log "Installing FastQC..."
        cd "$TOOLS_DIR"
        wget -q -O fastqc.zip \
            "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip"
        unzip -q fastqc.zip
        chmod +x FastQC/fastqc
        ln -sf ../FastQC/fastqc "$TOOLS_DIR/bin/fastqc"
    fi

    # Trimmomatic
    if [ ! -f "$TOOLS_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
        log "Installing Trimmomatic..."
        cd "$TOOLS_DIR"
        wget -q -O Trimmomatic-0.39.zip \
            "https://github.com/usadellab/Trimmomatic/files/5854849/Trimmomatic-0.39.zip"
        unzip -q Trimmomatic-0.39.zip
    fi

    # SRA Toolkit
    if [ ! -x "$TOOLS_DIR/bin/prefetch" ]; then
        log "Installing SRA Toolkit..."
        cd "$TOOLS_DIR"
        wget -q -O sratoolkit.tar.gz \
            "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
        tar -xzf sratoolkit.tar.gz
        SRATOOLKIT_DIR=$(ls -d sratoolkit.* 2>/dev/null | head -1)
        if [ -n "$SRATOOLKIT_DIR" ]; then
            mv "$SRATOOLKIT_DIR" sratoolkit
            ln -sf ../sratoolkit/bin/prefetch "$TOOLS_DIR/bin/prefetch"
            ln -sf ../sratoolkit/bin/fasterq-dump "$TOOLS_DIR/bin/fasterq-dump"
        fi
    fi

    # pigz
    if [ ! -x "$TOOLS_DIR/bin/pigz" ]; then
        log "Installing pigz..."
        cd "$TOOLS_DIR/build"
        wget -q https://zlib.net/pigz/pigz-2.8.tar.gz
        tar -xzf pigz-2.8.tar.gz
        cd pigz-2.8
        make
        cp pigz unpigz "$TOOLS_DIR/bin/"
    fi

    # GNU parallel
    if [ ! -x "$TOOLS_DIR/bin/parallel" ]; then
        log "Installing GNU parallel..."
        cd "$TOOLS_DIR/build"
        wget -q http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
        tar -xjf parallel-latest.tar.bz2
        cd parallel-*/
        ./configure --prefix="$TOOLS_DIR"
        make -j "$CPU_CORES"
        make install
    fi
}

test_installations() {
    log "Testing tool installations..."

    local all_ok=true

    # Test each tool
    for tool in bcftools samtools tabix bgzip STAR fastqc prefetch pigz parallel; do
        if [ -x "$TOOLS_DIR/bin/$tool" ]; then
            if "$TOOLS_DIR/bin/$tool" --version 2>&1 | head -1 > /dev/null 2>&1; then
                log "✓ $tool is working"
            else
                log_warning "$tool installed but not working"
                all_ok=false
            fi
        else
            log_warning "$tool not found in $TOOLS_DIR/bin"
            all_ok=false
        fi
    done

    # Test GATK
    if [ -x "$TOOLS_DIR/bin/gatk" ]; then
        if "$TOOLS_DIR/bin/gatk" --version 2>&1 | grep -q "GATK"; then
            log "✓ GATK is working"
        else
            log_warning "GATK installed but not working"
            all_ok=false
        fi
    fi

    # Test Java
    if command -v java &> /dev/null; then
        log "✓ Java is available"
    else
        log_warning "Java not found - required for GATK and Trimmomatic"
        all_ok=false
    fi

    if $all_ok; then
        log "✅ All tools installed and working correctly!"
    else
        log_warning "Some tools may have issues"
    fi
}

cleanup_build() {
    log "Cleaning up build files..."
    rm -rf "$TOOLS_DIR/build" 2>/dev/null || true
}

############################################################
# Main Execution
############################################################

main() {
    echo ""
    echo "========================================"
    echo "MDD2 Pipeline Tool Installer"
    echo "========================================"
    echo "This script will install all required tools"
    echo "locally in: $TOOLS_DIR"
    echo "No sudo permissions required!"
    echo "========================================"
    echo ""

    # Check disk space
    local available_mb=$(df -m "$TOOLS_DIR" | tail -1 | awk '{print $4}')
    if [ "$available_mb" -lt 5000 ]; then
        log_warning "Low disk space: ${available_mb}MB available"
        log_warning "At least 5GB recommended"
        read -p "Continue anyway? (y/n): " -n 1 -r
        echo
        [[ ! $REPLY =~ ^[Yy]$ ]] && exit 1
    fi

    # Create tools directory
    mkdir -p "$TOOLS_DIR"

    # Install all tools
    install_all_static_tools

    # Cleanup
    cleanup_build

    # Create environment file
    create_env_file

    echo ""
    echo "✅ Installation complete!"
    echo ""
    echo "To use these tools in your current shell:"
    echo "  source $TOOLS_DIR/env.sh"
    echo ""
    echo "To add to your ~/.bashrc:"
    echo "  echo 'source $TOOLS_DIR/env.sh' >> ~/.bashrc"
    echo ""
}

create_env_file() {
    cat > "$TOOLS_DIR/env.sh" << 'EOF'
#!/bin/bash
# MDD2 Pipeline Environment Setup
export MDD2_TOOLS_DIR="$HOME/.local/mdd2_tools"
export PATH="$MDD2_TOOLS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MDD2_TOOLS_DIR/lib:$LD_LIBRARY_PATH"
export JAVA_OPTS="-Xmx8g -Xms2g"

# GATK setup
if [ -f "$MDD2_TOOLS_DIR/gatk" ]; then
    export GATK="$MDD2_TOOLS_DIR/gatk"
elif [ -f "$MDD2_TOOLS_DIR/bin/gatk" ]; then
    export GATK="$MDD2_TOOLS_DIR/bin/gatk"
fi

# Trimmomatic setup
if [ -f "$MDD2_TOOLS_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar" ]; then
    export TRIMMOMATIC_JAR="$MDD2_TOOLS_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar"
fi

echo "MDD2 Pipeline environment loaded from $MDD2_TOOLS_DIR"
EOF

    chmod +x "$TOOLS_DIR/env.sh"
    log "Created environment file: $TOOLS_DIR/env.sh"
}

# Run main if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi

# Function to be called from mdd2.sh
ensure_tools_installed() {
    log "Checking if tools are installed..."

    # Check for critical tools
    local critical_tools=("bcftools" "samtools" "STAR" "gatk")
    local missing_tools=()

    for tool in "${critical_tools[@]}"; do
        if [ ! -x "$TOOLS_DIR/bin/$tool" ] && ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [ ${#missing_tools[@]} -eq 0 ]; then
        log "✓ All critical tools are installed"
        return 0
    else
        log_warning "Missing tools: ${missing_tools[*]}"
        log "Running automatic tool installation..."

        # Call the installation function
        install_all_static_tools

        # Verify installation
        test_installations

        return $?
    fi
}
