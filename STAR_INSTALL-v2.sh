#!/bin/bash
# STAR RNA-seq Pipeline Setup Script

echo "Setting up RNA-seq pipeline..."

# Install required Python packages without --user in virtualenv
if [ -n "$VIRTUAL_ENV" ]; then
    echo "Installing packages in virtual environment..."
    pip install PyYAML gdown
else
    echo "Installing packages for user..."
    pip install --user PyYAML gdown
fi

# Add ~/.local/bin to PATH if not already there
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
    export PATH="$HOME/.local/bin:$PATH"
    echo "Added ~/.local/bin to PATH"
fi

# Create project directory structure
if [ ! -d "rna_seq_project" ]; then
    mkdir rna_seq_project
fi
cd rna_seq_project
mkdir -p samples results downloads temp local_tools/bin local_tools/lib local_tools/include BAM_Files

# Get download links file from user
echo "Enter the path to the download links text file (or press Enter to use default Google Drive):"
read download_links_file

if [ -z "$download_links_file" ]; then
    download_links_file="default"
    echo "Using default Google Drive URL"
else
    if [ ! -f "$download_links_file" ]; then
        echo "Warning: File '$download_links_file' not found. Using default Google Drive."
        download_links_file="default"
    else
        echo "Using download links from: $download_links_file"
        # Copy the links file to project directory for reference
        cp "$download_links_file" "download_links.txt"
    fi
fi

# Create config.yaml file with local paths
cat > config.yaml << EOF
# STAR RNA-seq Pipeline Configuration
project:
  name: "brain_rna_study"
  organism: "human"  # human, mouse, rat, etc.
  version: "auto"    # "auto" for latest, or specific version like "49"

paths:
  install_dir: "./local_tools/STAR"
  genome_dir: "./genome_index"
  input_dir: "./samples"
  output_dir: "./results"
  temp_dir: "./temp"
  download_dir: "./downloads"
  local_tools: "./local_tools"
  bam_files_dir: "./BAM_Files"

resources:
  threads: 8
  memory_gb: 16
  bam_sort_ram: 4000000000

data:
  google_drive_url: "https://drive.google.com/drive/folders/1Yec-8GAchIhwFeVa4LtY0BfxlY3Ttt5b?usp=sharing"
  download_files: true
  auto_organize: true
  download_links_file: "$download_links_file"

star:
  sjdbOverhang: 99
  genomeSAindexNbases: 14
  genomeChrBinNbits: 18
  quantMode: "GeneCounts"
  outSAMtype: "BAM SortedByCoordinate"
  twoPass: true

alignment:
  read_files_pattern: "*_{1,2}.paired.fastq.gz"
  strandedness: "reverse"  # reverse, forward, or unstranded

upload:
  enabled: false
  destination_dir: ""
  drive_upload: false

email:
  notifications: false
  address: "ninansajeethphilip@gmail.com"

# Advanced settings - usually do not need modification
advanced:
  genome_base_urls:
    human: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
    mouse: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse"
  cleanup_temp: true
  keep_intermediate: false
EOF

# Create the Python pipeline script with robust STAR installation
cat > run_pipeline.py << 'EOF'
#!/usr/bin/env python3
"""
Professional STAR RNA-seq Pipeline with Automatic Data Download
Handles multiple samples with automatic latest genome version detection
"""

import os
import sys
import yaml
import glob
import argparse
import logging
import subprocess
import shutil
import re
from pathlib import Path
from datetime import datetime
import urllib.request
import json
from typing import Dict, List, Tuple, Optional
import gdown
import zipfile
import tarfile

class STARPipeline:
    def __init__(self, config_path: str = "config.yaml"):
        self.config = self.load_config(config_path)
        self.setup_directories()
        self.setup_logging()
        self.samples = {}
        self.star_binary = None
        self.samtools_available = False
        self.local_bin_path = self.config['paths']['local_tools'] / 'bin'

    def load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file"""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Expand user paths and convert to Path objects
        for path_key in ['install_dir', 'genome_dir', 'input_dir', 'output_dir', 'temp_dir', 'download_dir', 'local_tools', 'bam_files_dir']:
            config['paths'][path_key] = Path(config['paths'][path_key]).expanduser().absolute()

        return config

    def setup_directories(self):
        """Create necessary directories"""
        dirs = [
            self.config['paths']['install_dir'],
            self.config['paths']['genome_dir'],
            self.config['paths']['input_dir'],
            self.config['paths']['output_dir'],
            self.config['paths']['temp_dir'],
            self.config['paths']['download_dir'],
            self.config['paths']['local_tools'] / 'bin',
            self.config['paths']['local_tools'] / 'lib',
            self.config['paths']['local_tools'] / 'include',
            self.config['paths']['bam_files_dir']
        ]

        for directory in dirs:
            directory.mkdir(parents=True, exist_ok=True)

        # Add local tools to PATH for this session
        local_bin_path = str(self.config['paths']['local_tools'] / 'bin')
        if local_bin_path not in os.environ['PATH']:
            os.environ['PATH'] = local_bin_path + ':' + os.environ['PATH']

    def setup_logging(self):
        """Setup comprehensive logging"""
        log_dir = self.config['paths']['output_dir'] / "logs"
        log_dir.mkdir(exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"star_pipeline_{timestamp}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger("STAR_Pipeline")

        self.logger.info(f"STAR Pipeline initialized for project: {self.config['project']['name']}")

    def check_system_dependencies(self):
        """Check and install system dependencies"""
        self.logger.info("Checking system dependencies...")

        # Check for essential build tools
        build_tools = ['gcc', 'g++', 'make', 'git', 'wget']
        missing_tools = []

        for tool in build_tools:
            try:
                subprocess.run(f"which {tool}", shell=True, check=True, capture_output=True)
            except subprocess.CalledProcessError:
                missing_tools.append(tool)

        if missing_tools:
            self.logger.warning(f"Missing build tools: {missing_tools}")
            self.logger.info("Please install missing tools with: sudo yum install -y {' '.join(missing_tools)}")
            return False

        self.logger.info("All build tools are available")
        return True

    def install_star(self):
        """Install STAR with multiple fallback methods"""
        star_dir = self.config['paths']['install_dir']
        self.star_binary = star_dir / "source" / "STAR"

        # Try method 1: Compile from source with system libraries
        if self.install_star_compile_system():
            return

        # Try method 2: Download older compatible binary
        if self.install_star_older_binary():
            return

        # Try method 3: Use conda if available
        if self.install_star_conda():
            return

        # Final fallback: Try simple compilation
        if self.install_star_simple():
            return

        raise RuntimeError("All STAR installation methods failed")

    def install_star_compile_system(self):
        """Compile STAR from source using system libraries"""
        try:
            self.logger.info("Attempting to compile STAR from source with system libraries...")

            star_dir = self.config['paths']['install_dir']

            # Clean environment to avoid Anaconda conflicts
            env = os.environ.copy()
            # Remove Anaconda from environment
            for key in list(env.keys()):
                if 'anaconda' in env.get(key, '') or 'conda' in env.get(key, ''):
                    if key != 'PATH':
                        env.pop(key, None)

            # Clean PATH
            env['PATH'] = ':'.join([
                p for p in env['PATH'].split(':')
                if 'anaconda' not in p and 'conda' not in p
            ])

            # Use system libraries explicitly
            env['CXX'] = 'g++'
            env['CC'] = 'gcc'

            commands = [
                f"cd {star_dir.parent} && rm -rf STAR && git clone https://github.com/alexdobin/STAR.git",
                f"cd {star_dir}/source && make -j {min(4, self.config['resources']['threads'])} STAR"
            ]

            for cmd in commands:
                result = subprocess.run(
                    cmd, shell=True, check=True, capture_output=True, text=True, env=env
                )

            # Verify the binary works
            if self.verify_star_binary():
                self.logger.info("STAR compiled successfully with system libraries")
                return True
            else:
                self.logger.warning("STAR compilation completed but binary verification failed")
                return False

        except Exception as e:
            self.logger.warning(f"STAR compilation failed: {e}")
            return False

    def install_star_older_binary(self):
        """Download an older compatible STAR binary"""
        try:
            self.logger.info("Trying older compatible STAR binary...")

            star_dir = self.config['paths']['install_dir']
            self.star_binary = star_dir / "STAR_2.7.0"

            # Download STAR 2.7.0 which has better compatibility
            download_url = "https://github.com/alexdobin/STAR/releases/download/2.7.0f/STAR_Linux_x86_64"
            urllib.request.urlretrieve(download_url, self.star_binary)
            self.star_binary.chmod(0o755)

            if self.verify_star_binary():
                self.logger.info("Older STAR binary installed successfully")
                return True
            else:
                return False

        except Exception as e:
            self.logger.warning(f"Older binary installation failed: {e}")
            return False

    def install_star_conda(self):
        """Install STAR using conda if available"""
        try:
            self.logger.info("Trying conda installation...")

            # Check if conda is available
            result = subprocess.run("which conda", shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.info("Conda not available, skipping conda installation")
                return False

            star_dir = self.config['paths']['install_dir']
            self.star_binary = star_dir / "STAR"

            # Create conda environment
            env_name = f"star_env_{os.getpid()}"
            commands = [
                f"conda create -y -n {env_name} star=2.7.0",
                f"conda run -n {env_name} which STAR"
            ]

            for cmd in commands:
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.returncode != 0:
                    self.logger.warning(f"Conda command failed: {cmd}")
                    return False

            # Get the STAR binary path from conda
            result = subprocess.run(
                f"conda run -n {env_name} which STAR",
                shell=True, capture_output=True, text=True
            )

            if result.returncode == 0:
                conda_star_path = result.stdout.strip()
                # Copy to our local directory
                shutil.copy2(conda_star_path, self.star_binary)
                self.star_binary.chmod(0o755)

                if self.verify_star_binary():
                    self.logger.info("STAR installed successfully via conda")
                    return True

            return False

        except Exception as e:
            self.logger.warning(f"Conda installation failed: {e}")
            return False

    def install_star_simple(self):
        """Simple STAR installation as last resort"""
        try:
            self.logger.info("Trying simple STAR installation...")

            # Try to use system STAR if available
            result = subprocess.run("which STAR", shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                system_star = result.stdout.strip()
                star_dir = self.config['paths']['install_dir']
                self.star_binary = star_dir / "STAR"
                shutil.copy2(system_star, self.star_binary)

                if self.verify_star_binary():
                    self.logger.info("Using system STAR binary")
                    return True

            return False

        except Exception as e:
            self.logger.warning(f"Simple installation failed: {e}")
            return False

    def verify_star_binary(self):
        """Verify that STAR binary works"""
        try:
            result = subprocess.run(
                f"{self.star_binary} --version",
                shell=True, check=True, capture_output=True, text=True
            )
            self.logger.info(f"STAR verification successful: {result.stdout.strip()}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"STAR binary verification failed: {e}")
            return False

    def download_and_organize_data(self):
        """Download data and organize it"""
        if not self.config['data']['download_files']:
            self.logger.info("Data download disabled in config")
            return

        download_dir = self.config['paths']['download_dir']
        input_dir = self.config['paths']['input_dir']

        # Download files based on configuration
        if self.config['data']['download_links_file'] != "default":
            self.download_from_links_file(download_dir)
        else:
            self.download_from_google_drive(download_dir)

        # Organize downloaded files
        self.organize_downloaded_files(download_dir, input_dir)

    def download_from_links_file(self, download_dir: Path):
        """Download files from a text file containing URLs"""
        links_file = self.config['data']['download_links_file']
        if links_file == "default":
            self.logger.info("Using default Google Drive download method")
            return self.download_from_google_drive(download_dir)

        if not Path(links_file).exists():
            self.logger.warning(f"Download links file not found: {links_file}")
            self.logger.info("Falling back to default Google Drive download")
            return self.download_from_google_drive(download_dir)

        self.logger.info(f"Downloading files from links in: {links_file}")

        try:
            with open(links_file, 'r') as f:
                urls = [line.strip() for line in f if line.strip() and not line.startswith('#')]

            downloaded_count = 0
            for url in urls:
                try:
                    filename = url.split('/')[-1]
                    if not filename:
                        filename = f"downloaded_file_{downloaded_count}"

                    output_path = download_dir / filename

                    self.logger.info(f"Downloading: {filename}")

                    if 'drive.google.com' in url:
                        # Google Drive link
                        file_id = self.extract_google_drive_id(url)
                        if file_id:
                            gdown.download(id=file_id, output=str(output_path), quiet=False)
                            downloaded_count += 1
                        else:
                            self.logger.warning(f"Could not extract Google Drive ID from: {url}")
                    else:
                        # Regular URL
                        urllib.request.urlretrieve(url, output_path)
                        downloaded_count += 1

                except Exception as e:
                    self.logger.error(f"Failed to download {url}: {e}")

            self.logger.info(f"Downloaded {downloaded_count} files from links file")

            # Extract any downloaded archives
            self.extract_zip_files(download_dir)

        except Exception as e:
            self.logger.error(f"Error processing download links file: {e}")
            self.logger.info("Falling back to default Google Drive download")
            self.download_from_google_drive(download_dir)

    def extract_google_drive_id(self, url: str) -> str:
        """Extract file ID from Google Drive URL"""
        patterns = [
            r'/file/d/([^/]+)',
            r'id=([^&]+)',
            r'/folders/([^/?]+)'
        ]

        for pattern in patterns:
            match = re.search(pattern, url)
            if match:
                return match.group(1)
        return None

    def download_from_google_drive(self, download_dir: Path):
        """Download data from Google Drive"""
        self.logger.info("Downloading data from Google Drive...")

        # Google Drive folder ID extracted from the URL
        folder_id = "1Yec-8GAchIhwFeVa4LtY0BfxlY3Ttt5b"

        try:
            # Download the entire folder
            self.logger.info("Starting download from Google Drive...")
            gdown.download_folder(
                id=folder_id,
                output=str(download_dir),
                quiet=False,
                use_cookies=False
            )

            self.logger.info("Download completed successfully")

            # Extract ZIP files if any
            self.extract_zip_files(download_dir)

        except Exception as e:
            self.logger.error(f"Failed to download data from Google Drive: {e}")
            self.logger.info("Attempting alternative download method...")
            self.alternative_download_method(download_dir)

    def extract_zip_files(self, download_dir: Path):
        """Extract all ZIP files found in download directory"""
        self.logger.info("Looking for ZIP files to extract...")

        zip_files = list(download_dir.rglob("*.zip"))
        tar_files = list(download_dir.rglob("*.tar.gz")) + list(download_dir.rglob("*.tgz"))

        all_archives = zip_files + tar_files

        if not all_archives:
            self.logger.info("No archive files found to extract")
            return

        for archive_path in all_archives:
            self.logger.info(f"Extracting archive: {archive_path.name}")
            try:
                # Create extraction directory with archive file name
                extract_dir = download_dir / archive_path.stem
                extract_dir.mkdir(exist_ok=True)

                if archive_path.suffix == '.zip':
                    with zipfile.ZipFile(archive_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_dir)
                        self.logger.info(f"Extracted {len(zip_ref.namelist())} files to {extract_dir}")

                elif archive_path.suffix in ['.gz', '.tgz'] or '.tar.gz' in archive_path.name:
                    with tarfile.open(archive_path, 'r:gz') as tar_ref:
                        tar_ref.extractall(extract_dir)
                        self.logger.info(f"Extracted files to {extract_dir}")

            except (zipfile.BadZipFile, tarfile.TarError) as e:
                self.logger.error(f"Bad archive file: {archive_path} - {e}")
            except Exception as e:
                self.logger.error(f"Failed to extract {archive_path}: {e}")

    def alternative_download_method(self, download_dir: Path):
        """Alternative method for downloading data"""
        try:
            # Try direct download with gdown for specific files
            self.logger.info("Trying direct download with gdown...")

            # List of potential file IDs (you might need to get these from your Google Drive)
            file_ids = {
                "data.zip": "1Yec-8GAchIhwFeVa4LtY0BfxlY3Ttt5b"
            }

            for filename, file_id in file_ids.items():
                output_path = download_dir / filename
                if not output_path.exists():
                    self.logger.info(f"Downloading {filename}...")
                    gdown.download(id=file_id, output=str(output_path), quiet=False)

                    # Extract if it's an archive file
                    self.extract_zip_files(download_dir)

        except Exception as e:
            self.logger.error(f"Alternative download also failed: {e}")
            self.logger.info("Please manually download files and place in samples/ directory")

    def organize_downloaded_files(self, download_dir: Path, input_dir: Path):
        """Organize downloaded files into the samples directory"""
        self.logger.info("Organizing downloaded files...")

        # Find all downloaded FASTQ files (including in extracted directories)
        fastq_patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
        downloaded_files = []

        for pattern in fastq_patterns:
            # Search in download_dir and all subdirectories
            downloaded_files.extend(download_dir.rglob(pattern))

        self.logger.info(f"Found {len(downloaded_files)} FASTQ files")

        if not downloaded_files:
            self.logger.warning("No FASTQ files found after extraction")
            self.logger.info("Available files in download directory:")
            for item in download_dir.rglob("*"):
                if item.is_file():
                    self.logger.info(f"  - {item.relative_to(download_dir)}")
            return

        # Organize files based on naming patterns
        paired_end_files = self.identify_paired_end_files(downloaded_files)

        # Move and rename files to standard format
        moved_count = 0
        for sample_name, (read1, read2) in paired_end_files.items():
            # Create standardized names
            new_read1 = input_dir / f"{sample_name}_1.paired.fastq.gz"
            new_read2 = input_dir / f"{sample_name}_2.paired.fastq.gz"

            # Move files (copy instead of move to keep originals)
            try:
                shutil.copy2(read1, new_read1)
                shutil.copy2(read2, new_read2)
                moved_count += 2
                self.logger.info(f"Organized: {sample_name} -> {new_read1.name}, {new_read2.name}")
            except Exception as e:
                self.logger.error(f"Failed to copy files for {sample_name}: {e}")

        if moved_count > 0:
            self.logger.info(f"Successfully organized {moved_count} files into {len(paired_end_files)} sample pairs")
        else:
            self.logger.warning("No files were organized. Manual organization may be needed.")
            self.logger.info("Please ensure FASTQ files follow naming patterns like:")
            self.logger.info("  - sample_R1.fastq.gz, sample_R2.fastq.gz")
            self.logger.info("  - sample_1.fastq.gz, sample_2.fastq.gz")
            self.logger.info("  - sample.read1.fastq.gz, sample.read2.fastq.gz")

    def identify_paired_end_files(self, file_list: List[Path]) -> Dict[str, Tuple[Path, Path]]:
        """Identify paired-end files from a list of files"""
        paired_files = {}

        # Common paired-end naming patterns
        patterns = [
            (r'(.+)_[Rr]1[._].*', r'(.+)_[Rr]2[._].*'),  # sample_R1.fastq.gz, sample_R2.fastq.gz
            (r'(.+)_1[._].*', r'(.+)_2[._].*'),          # sample_1.fastq.gz, sample_2.fastq.gz
            (r'(.+)[._]1[._].*', r'(.+)[._]2[._].*'),    # sample.1.fastq.gz, sample.2.fastq.gz
            (r'(.+)_read1[._].*', r'(.+)_read2[._].*'),  # sample_read1.fastq.gz, sample_read2.fastq.gz
            (r'(.+)[._]read1[._].*', r'(.+)[._]read2[._].*'),  # sample.read1.fastq.gz, sample.read2.fastq.gz
            (r'(.+)[._]R1[._].*', r'(.+)[._]R2[._].*'),  # sample.R1.fastq.gz, sample.R2.fastq.gz
        ]

        for pattern1, pattern2 in patterns:
            read1_files = {}
            read2_files = {}

            for file_path in file_list:
                filename = file_path.name

                # Check for read1 pattern
                match1 = re.match(pattern1, filename)
                if match1:
                    sample_name = match1.group(1)
                    # Clean up sample name (remove common suffixes)
                    sample_name = re.sub(r'[._](fastq|fq|gz)$', '', sample_name)
                    read1_files[sample_name] = file_path
                    continue

                # Check for read2 pattern
                match2 = re.match(pattern2, filename)
                if match2:
                    sample_name = match2.group(1)
                    # Clean up sample name (remove common suffixes)
                    sample_name = re.sub(r'[._](fastq|fq|gz)$', '', sample_name)
                    read2_files[sample_name] = file_path

            # Find matching pairs
            common_samples = set(read1_files.keys()) & set(read2_files.keys())
            for sample_name in common_samples:
                if sample_name not in paired_files:
                    paired_files[sample_name] = (read1_files[sample_name], read2_files[sample_name])

        self.logger.info(f"Identified {len(paired_files)} paired-end samples")
        return paired_files

    def get_latest_gencode_version(self) -> str:
        """Automatically detect the latest GENCODE version"""
        organism = self.config['project']['organism']
        base_url = self.config['advanced']['genome_base_urls'][organism]

        self.logger.info(f"Detecting latest GENCODE version for {organism}...")

        try:
            # Try to get the listing page
            with urllib.request.urlopen(base_url) as response:
                html = response.read().decode('utf-8')

            # Extract version numbers from directory listings
            versions = re.findall(r'release_(\d+)', html)
            if versions:
                latest_version = max(int(v) for v in versions)
                self.logger.info(f"Latest version detected: {latest_version}")
                return str(latest_version)
            else:
                self.logger.warning("Could not auto-detect version, using fallback")
                return "49"  # Fallback

        except Exception as e:
            self.logger.warning(f"Version detection failed: {e}, using fallback")
            return "49"  # Fallback

    def discover_samples(self) -> Dict[str, Tuple[Path, Path]]:
        """Discover sample pairs from input directory"""
        input_dir = self.config['paths']['input_dir']
        pattern = self.config['alignment']['read_files_pattern']

        # Find all read1 files
        read1_files = glob.glob(str(input_dir / pattern.replace("{1,2}", "1")))
        samples = {}

        for read1_path in read1_files:
            read1 = Path(read1_path)
            sample_name = read1.name.replace("_1.paired.fastq.gz", "").replace("_1.fastq.gz", "")

            # Find corresponding read2
            read2_pattern = read1_path.replace("_1.paired.fastq.gz", "_2.paired.fastq.gz").replace("_1.fastq.gz", "_2.fastq.gz")
            read2 = Path(read2_pattern)

            if read2.exists():
                samples[sample_name] = (read1, read2)
                self.logger.info(f"Found sample: {sample_name} -> {read1.name}, {read2.name}")
            else:
                self.logger.warning(f"Missing read2 for sample {sample_name}")

        if not samples:
            self.logger.error("No sample pairs found! Check input directory and file patterns.")
            self.logger.info("Available files in input directory:")
            for file in input_dir.iterdir():
                self.logger.info(f"  - {file.name}")
            sys.exit(1)

        self.logger.info(f"Discovered {len(samples)} sample pairs")
        return samples

    def download_reference_files(self) -> Tuple[Path, Path]:
        """Download reference genome and annotations"""
        organism = self.config['project']['organism']
        version = self.config['project']['version']

        if version == "auto":
            version = self.get_latest_gencode_version()

        self.logger.info(f"Using GENCODE version: {version} for {organism}")

        base_url = self.config['advanced']['genome_base_urls'][organism]
        release_url = f"{base_url}/release_{version}"

        # Define file names based on organism
        if organism == "human":
            genome_file = f"GRCh38.primary_assembly.genome.fa"
            annotation_file = f"gencode.v{version}.annotation.gtf"
        elif organism == "mouse":
            genome_file = f"GRCm39.primary_assembly.genome.fa"
            annotation_file = f"gencode.v{version}.annotation.gtf"
        else:
            self.logger.error(f"Unsupported organism: {organism}")
            sys.exit(1)

        genome_path = self.config['paths']['genome_dir'] / genome_file
        annotation_path = self.config['paths']['genome_dir'] / annotation_file

        # Download files if they don't exist
        if not genome_path.exists():
            self.logger.info(f"Downloading genome: {genome_file}")
            genome_url = f"{release_url}/{genome_file}.gz"
            self.download_file(genome_url, genome_path)

        if not annotation_path.exists():
            self.logger.info(f"Downloading annotations: {annotation_file}")
            annotation_url = f"{release_url}/{annotation_file}.gz"
            self.download_file(annotation_url, annotation_path)

        return genome_path, annotation_path

    def download_file(self, url: str, output_path: Path):
        """Download and extract gzipped file"""
        try:
            # Download
            temp_path = output_path.with_suffix('.gz')
            self.logger.info(f"Downloading from: {url}")
            urllib.request.urlretrieve(url, temp_path)

            # Extract
            self.logger.info(f"Extracting to: {output_path}")
            subprocess.run(f"gunzip -c {temp_path} > {output_path}", shell=True, check=True)
            temp_path.unlink()

            self.logger.info(f"Downloaded and extracted: {output_path.name}")

        except Exception as e:
            self.logger.error(f"Failed to download {url}: {e}")
            sys.exit(1)

    def run_star_command(self, star_args: str, description: str) -> subprocess.CompletedProcess:
        """Run STAR command with clean environment"""
        self.logger.info(f"{description}")

        # Clean up the arguments
        star_args = ' '.join(star_args.split())

        # Use clean environment to avoid library conflicts
        env = os.environ.copy()
        # Remove problematic environment variables
        env.pop('LD_LIBRARY_PATH', None)
        env.pop('LD_PRELOAD', None)
        # Clean PATH from conda/anaconda
        env['PATH'] = ':'.join([
            p for p in env['PATH'].split(':')
            if 'anaconda' not in p and 'conda' not in p
        ])

        cmd = f"{self.star_binary} {star_args}"

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                executable='/bin/bash',
                env=env
            )
            return result
        except subprocess.CalledProcessError as e:
            self.logger.error(f"STAR command failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise

    def generate_genome_index(self, genome_path: Path, annotation_path: Path):
        """Generate STAR genome index"""
        genome_dir = self.config['paths']['genome_dir']

        # Check if index already exists
        if (genome_dir / "SAindex").exists():
            self.logger.info("Genome index already exists")
            return

        self.logger.info("Generating genome index...")

        # Build the STAR command
        star_args = (
            f"--runThreadN {min(4, self.config['resources']['threads'])} "
            f"--runMode genomeGenerate "
            f"--genomeDir {genome_dir} "
            f"--genomeFastaFiles {genome_path} "
            f"--sjdbGTFfile {annotation_path} "
            f"--sjdbOverhang {self.config['star']['sjdbOverhang']} "
            f"--genomeSAindexNbases {self.config['star']['genomeSAindexNbases']} "
            f"--genomeChrBinNbits {self.config['star']['genomeChrBinNbits']} "
            f"--limitGenomeGenerateRAM {self.config['resources']['memory_gb'] * 1024**3}"
        )

        self.run_star_command(star_args, "Generating genome index")
        self.logger.info("Genome index generation completed")

    def run_sample_alignment(self, sample_name: str, read1: Path, read2: Path):
        """Run two-pass alignment for a single sample"""
        self.logger.info(f"Processing sample: {sample_name}")

        sample_output_dir = self.config['paths']['output_dir'] / sample_name
        sample_output_dir.mkdir(exist_ok=True)

        genome_dir = self.config['paths']['genome_dir']

        # Use STAR's built-in BAM sorting
        self.logger.info("Using STAR's built-in BAM sorting")
        self.run_alignment_star_bam(sample_name, read1, read2, genome_dir, sample_output_dir)

        # Generate summary statistics
        self.generate_sample_summary(sample_name, sample_output_dir)

        # Copy BAM file to BAM_Files directory
        self.copy_bam_to_collection(sample_name, sample_output_dir)

        self.logger.info(f"Completed processing sample: {sample_name}")

    def copy_bam_to_collection(self, sample_name: str, sample_output_dir: Path):
        """Copy final BAM file to BAM_Files directory"""
        bam_files_dir = self.config['paths']['bam_files_dir']

        # Find the final BAM file
        bam_patterns = [
            f"{sample_name}_final_Aligned.sortedByCoord.out.bam",
            f"*_final_Aligned.sortedByCoord.out.bam"
        ]

        bam_file = None
        for pattern in bam_patterns:
            matches = list(sample_output_dir.glob(pattern))
            if matches:
                bam_file = matches[0]
                break

        if bam_file and bam_file.exists():
            dest_bam = bam_files_dir / f"{sample_name}.sorted.bam"
            dest_bai = bam_files_dir / f"{sample_name}.sorted.bam.bai"

            # Copy BAM file
            shutil.copy2(bam_file, dest_bam)
            self.logger.info(f"Copied BAM file to: {dest_bam}")

            # Copy BAM index if exists
            bai_file = bam_file.with_suffix('.bam.bai')
            if not bai_file.exists():
                bai_file = Path(str(bam_file) + '.bai')

            if bai_file.exists():
                shutil.copy2(bai_file, dest_bai)
                self.logger.info(f"Copied BAM index to: {dest_bai}")
        else:
            self.logger.warning(f"Could not find BAM file for sample: {sample_name}")

    def run_alignment_star_bam(self, sample_name: str, read1: Path, read2: Path, genome_dir: Path, sample_output_dir: Path):
        """Run alignment using STAR's built-in BAM sorting"""
        # First pass
        self.logger.info(f"First pass for {sample_name}")
        star_args1 = (
            f"--runThreadN {min(4, self.config['resources']['threads'])} "
            f"--genomeDir {genome_dir} "
            f"--readFilesIn {read1} {read2} "
            f"--readFilesCommand zcat "
            f"--outSAMtype BAM Unsorted "
            f"--outTmpDir {self.config['paths']['temp_dir'] / f'{sample_name}_pass1'} "
            f"--outFileNamePrefix {sample_output_dir / f'{sample_name}_pass1_'}"
        )
        self.run_star_command(star_args1, f"First pass - {sample_name}")

        # Second pass with sorted BAM output
        self.logger.info(f"Second pass for {sample_name}")
        star_args2 = (
            f"--runThreadN {min(4, self.config['resources']['threads'])} "
            f"--genomeDir {genome_dir} "
            f"--readFilesIn {read1} {read2} "
            f"--readFilesCommand zcat "
            f"--sjdbFileChrStartEnd {sample_output_dir / f'{sample_name}_pass1_SJ.out.tab'} "
            f"--quantMode {self.config['star']['quantMode']} "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outTmpDir {self.config['paths']['temp_dir'] / f'{sample_name}_pass2'} "
            f"--outSAMattrRGline ID:{sample_name} SM:{sample_name} "
            f"--outFileNamePrefix {sample_output_dir / f'{sample_name}_final_'}"
        )
        self.run_star_command(star_args2, f"Second pass - {sample_name}")

    def generate_sample_summary(self, sample_name: str, sample_dir: Path):
        """Generate summary statistics for a sample"""
        log_file = sample_dir / f"{sample_name}_final_Log.final.out"
        summary_file = sample_dir / f"{sample_name}_alignment_summary.txt"

        if log_file.exists():
            shutil.copy2(log_file, summary_file)

    def run_command(self, command: str, description: str) -> subprocess.CompletedProcess:
        """Run a shell command with logging and error handling"""
        self.logger.info(f"{description}")

        try:
            result = subprocess.run(
                command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                executable='/bin/bash'
            )
            return result
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise

    def generate_project_report(self):
        """Generate a comprehensive project report"""
        report_file = self.config['paths']['output_dir'] / "project_report.md"

        with open(report_file, 'w') as f:
            f.write("# STAR RNA-seq Alignment Report\n\n")
            f.write(f"**Project:** {self.config['project']['name']}\n")
            f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Organism:** {self.config['project']['organism']}\n")
            f.write(f"**Samples Processed:** {len(self.samples)}\n\n")

            f.write("## Sample Summary\n\n")
            f.write("| Sample | Status | Output Directory |\n")
            f.write("|--------|--------|------------------|\n")

            for sample_name in self.samples:
                sample_dir = self.config['paths']['output_dir'] / sample_name
                bam_file = sample_dir / f"{sample_name}_final_Aligned.sortedByCoord.out.bam"
                status = "✓ Complete" if bam_file.exists() else "✗ Failed"
                f.write(f"| {sample_name} | {status} | {sample_dir} |\n")

            f.write("\n## BAM Files Collection\n")
            bam_files_dir = self.config['paths']['bam_files_dir']
            bam_files = list(bam_files_dir.glob("*.bam"))
            if bam_files:
                f.write("The following BAM files are available in the BAM_Files directory:\n")
                for bam_file in bam_files:
                    f.write(f"- {bam_file.name}\n")
            else:
                f.write("No BAM files found in BAM_Files directory.\n")

        self.logger.info(f"Project report generated: {report_file}")

    def cleanup(self):
        """Clean up temporary files"""
        if self.config['advanced']['cleanup_temp']:
            self.logger.info("Cleaning up temporary files...")
            if self.config['paths']['temp_dir'].exists():
                shutil.rmtree(self.config['paths']['temp_dir'])

    def run_pipeline(self):
        """Main pipeline execution"""
        try:
            # Check system dependencies first
            if not self.check_system_dependencies():
                self.logger.warning("Some system dependencies may be missing, but continuing...")

            # Pipeline steps
            self.download_and_organize_data()
            self.install_star()

            genome_path, annotation_path = self.download_reference_files()

            self.generate_genome_index(genome_path, annotation_path)

            self.samples = self.discover_samples()

            # Process each sample
            for sample_name, (read1, read2) in self.samples.items():
                self.run_sample_alignment(sample_name, read1, read2)

            # Generate reports
            self.generate_project_report()

            # Cleanup
            self.cleanup()

            self.logger.info("STAR pipeline completed successfully!")

        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(description="Professional STAR RNA-seq Pipeline")
    parser.add_argument("--config", default="config.yaml", help="Path to config file")
    parser.add_argument("--samples-only", action="store_true", help="Only process samples (skip installation and genome setup)")
    parser.add_argument("--download-only", action="store_true", help="Only download and organize data")

    args = parser.parse_args()

    pipeline = STARPipeline(args.config)

    if args.download_only:
        pipeline.download_and_organize_data()
    elif args.samples_only:
        pipeline.samples = pipeline.discover_samples()
        for sample_name, (read1, read2) in pipeline.samples.items():
            pipeline.run_sample_alignment(sample_name, read1, read2)
        pipeline.generate_project_report()
    else:
        pipeline.run_pipeline()

if __name__ == "__main__":
    main()
EOF

# Make the Python script executable
chmod +x run_pipeline.py

# Create a manual installation script as backup
cat > manual_install_star.sh << 'EOF'
#!/bin/bash
echo "Manual STAR installation script"
echo "This script will help you install STAR manually if the automated methods fail"

# Check if we're in the right directory
if [ ! -f "config.yaml" ]; then
    echo "Please run this script from the rna_seq_project directory"
    exit 1
fi

echo "Installing STAR manually..."

# Method 1: Try system package manager
echo "Method 1: Trying system package manager..."
if command -v dnf &> /dev/null; then
    sudo dnf install -y STAR
elif command -v yum &> /dev/null; then
    sudo yum install -y STAR
elif command -v apt &> /dev/null; then
    sudo apt update && sudo apt install -y star
else
    echo "No supported package manager found"
fi

# Check if STAR is now available
if command -v STAR &> /dev/null; then
    echo "STAR installed successfully via package manager"
    STAR_PATH=$(which STAR)
    mkdir -p local_tools/STAR/source
    cp "$STAR_PATH" local_tools/STAR/source/STAR
    echo "STAR copied to local_tools/STAR/source/STAR"
    exit 0
fi

# Method 2: Download and compile
echo "Method 2: Downloading and compiling STAR..."
mkdir -p star_build
cd star_build
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
make -j 4 STAR
if [ -f "STAR" ]; then
    cp STAR ../../../local_tools/STAR/source/STAR
    echo "STAR compiled and copied successfully"
    cd ../..
    rm -rf star_build
    exit 0
fi

echo "Manual installation failed. Please install STAR manually and place it in local_tools/STAR/source/STAR"
echo "You can download pre-compiled binaries from: https://github.com/alexdobin/STAR/releases"
EOF

chmod +x manual_install_star.sh

echo "Setup completed successfully!"
echo ""
echo "This version includes multiple STAR installation methods:"
echo "1. Compilation from source with system libraries"
echo "2. Older compatible binary (2.7.0)"
echo "3. Conda installation (if available)"
echo "4. System STAR binary"
echo ""
echo "If automatic installation fails, run: ./manual_install_star.sh"
echo ""
echo "To run the pipeline:"
echo "  cd rna_seq_project && python3 run_pipeline.py"
