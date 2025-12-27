#!/usr/bin/env python3

"""
MDD2 Differential SNP Analyzer - Multi-Directory Version
Handles separate Control and MDD directories with nested TSV structure
AND performs actual differential SNP analysis
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path
import time
from datetime import datetime
import argparse
import subprocess

def discover_project_structure(base_path):
    """Discover the project structure with separate Control/MDD directories"""
    print("\n" + "="*70)
    print("PROJECT STRUCTURE DISCOVERY")
    print("="*70)

    base_path = Path(base_path)

    # Look for common directory patterns
    found_dirs = []

    # Common patterns for MDD/Case directories
    mdd_patterns = ['analysis2', 'analysis', 'mdd', 'case', 'MDD', 'patients']
    ctrl_patterns = ['analysisCtrl', 'control', 'ctrl', 'controls', 'healthy']

    # Search for directories
    print("Searching for project directories...")

    # First, check current directory
    for item in base_path.iterdir():
        if item.is_dir():
            found_dirs.append(item)

    # Also check one level up
    parent_dir = base_path.parent
    for item in parent_dir.iterdir():
        if item.is_dir() and item not in found_dirs:
            found_dirs.append(item)

    # Categorize directories
    mdd_dirs = []
    ctrl_dirs = []
    other_dirs = []

    for dir_path in found_dirs:
        dir_name = dir_path.name.lower()

        # Check if it's an analysis directory by looking for data/tsv structure
        has_tsv_structure = False
        for tsv_path in dir_path.glob("**/tsv"):
            if tsv_path.is_dir():
                # Check if it has any .tsv.gz files in subdirectories
                tsv_files = list(tsv_path.glob("**/*.tsv.gz"))
                if tsv_files:
                    has_tsv_structure = True
                    break

        if has_tsv_structure:
            # Try to categorize
            if any(pattern in dir_name for pattern in mdd_patterns):
                mdd_dirs.append(dir_path)
            elif any(pattern in dir_name for pattern in ctrl_patterns):
                ctrl_dirs.append(dir_path)
            else:
                other_dirs.append(dir_path)

    print(f"\nFound {len(mdd_dirs)} potential MDD directories:")
    for d in mdd_dirs:
        print(f"  • {d}")

    print(f"\nFound {len(ctrl_dirs)} potential Control directories:")
    for d in ctrl_dirs:
        print(f"  • {d}")

    if other_dirs:
        print(f"\nFound {len(other_dirs)} other analysis directories:")
        for d in other_dirs[:5]:  # Show first 5
            print(f"  • {d}")
        if len(other_dirs) > 5:
            print(f"  ... and {len(other_dirs)-5} more")

    return mdd_dirs, ctrl_dirs, other_dirs

def collect_tsv_files_from_directory(directory, group_type):
    """Collect all TSV files from a directory with subdirectory structure"""
    tsv_files = []

    # Look for data/tsv pattern
    tsv_dirs = []

    # Try common patterns
    patterns = [
        directory / "data" / "tsv",
        directory / "tsv",
        directory / "analysis" / "data" / "tsv",
    ]

    for pattern in patterns:
        if pattern.exists():
            tsv_dirs.append(pattern)

    # Also search recursively
    for tsv_dir in directory.rglob("tsv"):
        if tsv_dir.is_dir() and tsv_dir not in tsv_dirs:
            tsv_dirs.append(tsv_dir)

    # Collect files from all tsv directories
    for tsv_dir in tsv_dirs:
        # Files directly in tsv directory
        tsv_files.extend(tsv_dir.glob("*.tsv.gz"))
        tsv_files.extend(tsv_dir.glob("*.tsv"))

        # Files in subdirectories (SRA ID directories)
        for subdir in tsv_dir.iterdir():
            if subdir.is_dir():
                tsv_files.extend(subdir.glob("*.tsv.gz"))
                tsv_files.extend(subdir.glob("*.tsv"))

    return tsv_files

def extract_sample_id(file_path):
    """Extract sample ID from file path"""
    # Try to get from filename
    filename = file_path.name
    if filename.endswith('.tsv.gz'):
        return filename[:-7]
    elif filename.endswith('.tsv'):
        return filename[:-4]
    else:
        # Try to get from parent directory name
        parent = file_path.parent.name
        if parent and parent not in ['tsv', 'data']:
            return parent
        else:
            # Generate a unique ID
            import hashlib
            return f"SAMPLE_{hashlib.md5(str(file_path).encode()).hexdigest()[:8]}"

def extract_snps_from_tsv(tsv_file):
    """Extract SNPs from a TSV file - ACTUAL IMPLEMENTATION"""
    snps = []

    try:
        if tsv_file.suffix == '.gz':
            with gzip.open(tsv_file, 'rt') as f:
                # Skip header
                next(f)
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                        gene = parts[7] if len(parts) > 7 else 'UNKNOWN'

                        # Create unique SNP key
                        snp_key = f"{chrom}:{pos}:{ref}:{alt}"
                        snps.append((snp_key, chrom, pos, rsid, ref, alt, gene))
        else:
            with open(tsv_file, 'r') as f:
                next(f)
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                        gene = parts[7] if len(parts) > 7 else 'UNKNOWN'
                        snp_key = f"{chrom}:{pos}:{ref}:{alt}"
                        snps.append((snp_key, chrom, pos, rsid, ref, alt, gene))

    except Exception as e:
        print(f"Error reading {tsv_file}: {e}")

    return snps

def perform_differential_analysis(master_dir, threshold=95, min_samples=3, top_n=15, logger=None):
    """Perform actual differential SNP analysis"""

    print("\n" + "="*70)
    print("MDD2 DIFFERENTIAL SNP ANALYSIS")
    print("="*70)

    if logger:
        logger.log("Starting differential SNP analysis...")

    master_path = Path(master_dir)

    # Paths
    tsv_dir = master_path / "data" / "tsv"
    project_map = master_path / "project_mapping.txt"
    output_dir = master_path / "diff_analysis"
    output_dir.mkdir(exist_ok=True)

    # Load project mapping
    sample_groups = {}
    with open(project_map, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sample_groups[parts[0]] = parts[1]

    print(f"Loaded {len(sample_groups)} sample mappings")
    if logger:
        logger.log(f"Loaded {len(sample_groups)} sample mappings")

    # Collect all TSV files
    tsv_files = list(tsv_dir.glob("*.tsv.gz")) + list(tsv_dir.glob("*.tsv"))
    print(f"Found {len(tsv_files)} TSV files")
    if logger:
        logger.log(f"Found {len(tsv_files)} TSV files")

    # Process each sample
    print("\nProcessing samples...")
    if logger:
        logger.log("Processing samples...")

    control_snps = defaultdict(list)  # SNP key -> list of samples that have it
    mdd_snps = defaultdict(list)
    snp_info = {}  # SNP key -> (chrom, pos, rsid, ref, alt, gene)

    for tsv_file in tsv_files:
        sample_id = extract_sample_id(tsv_file)
        group = sample_groups.get(sample_id, 'UNKNOWN')

        print(f"  {sample_id} ({group})...")
        if logger:
            logger.log(f"  Processing {sample_id} ({group})...")

        snps = extract_snps_from_tsv(tsv_file)

        for snp_key, chrom, pos, rsid, ref, alt, gene in snps:
            if group == 'Control':
                control_snps[snp_key].append(sample_id)
            elif group == 'MDD':
                mdd_snps[snp_key].append(sample_id)

            # Store SNP info
            if snp_key not in snp_info:
                snp_info[snp_key] = (chrom, pos, rsid, ref, alt, gene)

    # Calculate statistics
    print("\nCalculating frequencies...")
    if logger:
        logger.log("Calculating frequencies...")

    control_samples = len([g for g in sample_groups.values() if g == 'Control'])
    mdd_samples = len([g for g in sample_groups.values() if g == 'MDD'])

    min_control = max(min_samples, int(np.ceil(control_samples * threshold / 100)))
    min_mdd = max(min_samples, int(np.ceil(mdd_samples * threshold / 100)))

    print(f"Control samples: {control_samples}, threshold: ≥{min_control}")
    print(f"MDD samples: {mdd_samples}, threshold: ≥{min_mdd}")

    if logger:
        logger.log(f"Control samples: {control_samples}, threshold: ≥{min_control}")
        logger.log(f"MDD samples: {mdd_samples}, threshold: ≥{min_mdd}")

    # Find differential SNPs
    print("\nIdentifying differential SNPs...")
    if logger:
        logger.log("Identifying differential SNPs...")

    differential_snps = []

    all_snp_keys = set(list(control_snps.keys()) + list(mdd_snps.keys()))

    for snp_key in all_snp_keys:
        control_count = len(control_snps.get(snp_key, []))
        mdd_count = len(mdd_snps.get(snp_key, []))

        control_freq = control_count / control_samples if control_samples > 0 else 0
        mdd_freq = mdd_count / mdd_samples if mdd_samples > 0 else 0

        specific_to = ""

        # Control-specific: present in ≥threshold% of Control, ≤1 MDD
        if control_count >= min_control and mdd_count <= 1:
            specific_to = "CONTROL"
        # MDD-specific: present in ≥threshold% of MDD, ≤1 Control
        elif mdd_count >= min_mdd and control_count <= 1:
            specific_to = "MDD"
        else:
            continue

        chrom, pos, rsid, ref, alt, gene = snp_info[snp_key]

        differential_snps.append({
            'SNP_KEY': snp_key,
            'CHROM': chrom,
            'POS': pos,
            'RSID': rsid if rsid and rsid != '.' else 'NA',
            'REF': ref,
            'ALT': alt,
            'GENE': gene,
            'CONTROL_COUNT': control_count,
            'MDD_COUNT': mdd_count,
            'CONTROL_FREQ': round(control_freq, 4),
            'MDD_FREQ': round(mdd_freq, 4),
            'DIFFERENCE': round(abs(control_freq - mdd_freq), 4),
            'SPECIFIC_TO': specific_to
        })

    print(f"Found {len(differential_snps)} differential SNPs")
    if logger:
        logger.log(f"Found {len(differential_snps)} differential SNPs")

    # Convert to DataFrame
    df = pd.DataFrame(differential_snps)

    # Save results
    print("\nSaving results...")
    if logger:
        logger.log("Saving results...")

    # 1. All differential SNPs
    all_snps_file = output_dir / "all_differential_snps.tsv"
    df.to_csv(all_snps_file, sep='\t', index=False)
    print(f"  ✓ All SNPs: {all_snps_file}")
    if logger:
        logger.log(f"  ✓ All SNPs: {all_snps_file}")

    # 2. Top Control-specific SNPs
    control_df = pd.DataFrame()
    mdd_df = pd.DataFrame()

    if 'CONTROL' in df['SPECIFIC_TO'].values:
        control_df = df[df['SPECIFIC_TO'] == 'CONTROL'].sort_values('DIFFERENCE', ascending=False)
        control_top_file = output_dir / f"top_{top_n}_control_specific_snps.tsv"
        control_df.head(top_n).to_csv(control_top_file, sep='\t', index=False)
        print(f"  ✓ Top Control-specific: {control_top_file}")
        if logger:
            logger.log(f"  ✓ Top Control-specific: {control_top_file}")

    # 3. Top MDD-specific SNPs
    if 'MDD' in df['SPECIFIC_TO'].values:
        mdd_df = df[df['SPECIFIC_TO'] == 'MDD'].sort_values('DIFFERENCE', ascending=False)
        mdd_top_file = output_dir / f"top_{top_n}_mdd_specific_snps.tsv"
        mdd_df.head(top_n).to_csv(mdd_top_file, sep='\t', index=False)
        print(f"  ✓ Top MDD-specific: {mdd_top_file}")
        if logger:
            logger.log(f"  ✓ Top MDD-specific: {mdd_top_file}")

    # 4. Create Excel report
    try:
        excel_file = output_dir / "differential_snps_report.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            # Control-specific
            if not control_df.empty:
                control_df.head(top_n).to_excel(writer, sheet_name='Control_Specific', index=False)

            # MDD-specific
            if not mdd_df.empty:
                mdd_df.head(top_n).to_excel(writer, sheet_name='MDD_Specific', index=False)

            # All SNPs
            df.to_excel(writer, sheet_name='All_Differential_SNPs', index=False)

            # Summary
            summary_data = {
                'Metric': ['Analysis Date', 'Threshold (%)', 'Min Samples', 'Top N',
                          'Control Samples', 'MDD Samples', 'Total SNPs',
                          'Control-specific', 'MDD-specific'],
                'Value': [datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                         threshold, min_samples, top_n,
                         control_samples, mdd_samples, len(df),
                         len(control_df) if not control_df.empty else 0,
                         len(mdd_df) if not mdd_df.empty else 0]
            }
            pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

        print(f"  ✓ Excel report: {excel_file}")
        if logger:
            logger.log(f"  ✓ Excel report: {excel_file}")

    except Exception as e:
        print(f"  ✗ Excel report failed: {e}")
        if logger:
            logger.warning(f"  Excel report failed: {e}")

    # 5. Create comprehensive summary
    summary_file = output_dir / "analysis_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("MDD2 Differential SNP Analysis Results\n")
        f.write("="*70 + "\n\n")

        f.write("ANALYSIS PARAMETERS\n")
        f.write("-"*70 + "\n")
        f.write(f"Threshold: {threshold}%\n")
        f.write(f"Minimum samples per group: {min_samples}\n")
        f.write(f"Top N results: {top_n}\n\n")

        f.write("SAMPLE INFORMATION\n")
        f.write("-"*70 + "\n")
        f.write(f"Control samples: {control_samples}\n")
        f.write(f"MDD samples: {mdd_samples}\n")
        f.write(f"Total samples: {control_samples + mdd_samples}\n\n")

        f.write("RESULTS\n")
        f.write("-"*70 + "\n")
        f.write(f"Total differential SNPs: {len(df)}\n")
        f.write(f"Control-specific SNPs: {len(control_df) if not control_df.empty else 0}\n")
        f.write(f"MDD-specific SNPs: {len(mdd_df) if not mdd_df.empty else 0}\n\n")

        if not control_df.empty:
            f.write("TOP CONTROL-SPECIFIC SNPS\n")
            f.write("-"*70 + "\n")
            for _, row in control_df.head(5).iterrows():
                f.write(f"{row['CHROM']}:{row['POS']} {row['REF']}>{row['ALT']} ({row['GENE']})\n")
                f.write(f"  Control: {row['CONTROL_COUNT']}/{control_samples} ({row['CONTROL_FREQ']:.1%})\n")
                f.write(f"  MDD: {row['MDD_COUNT']}/{mdd_samples} ({row['MDD_FREQ']:.1%})\n\n")

        if not mdd_df.empty:
            f.write("TOP MDD-SPECIFIC SNPS\n")
            f.write("-"*70 + "\n")
            for _, row in mdd_df.head(5).iterrows():
                f.write(f"{row['CHROM']}:{row['POS']} {row['REF']}>{row['ALT']} ({row['GENE']})\n")
                f.write(f"  Control: {row['CONTROL_COUNT']}/{control_samples} ({row['CONTROL_FREQ']:.1%})\n")
                f.write(f"  MDD: {row['MDD_COUNT']}/{mdd_samples} ({row['MDD_FREQ']:.1%})\n\n")

    print(f"  ✓ Text summary: {summary_file}")
    if logger:
        logger.log(f"  ✓ Text summary: {summary_file}")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nResults saved in: {output_dir}")

    if logger:
        logger.log("Analysis complete!")
        logger.log(f"Results saved in: {output_dir}")

    return True

def setup_multi_directory_analysis():
    """Setup analysis with separate Control and MDD directories"""
    print("\n" + "="*70)
    print("MULTI-DIRECTORY ANALYSIS SETUP")
    print("="*70)
    print("\nYour project appears to have separate directories for Control and MDD samples.")
    print("Let's configure them properly.")

    base_path = Path.cwd()

    # Discover directories
    mdd_dirs, ctrl_dirs, other_dirs = discover_project_structure(base_path)

    # Configure MDD directory
    print("\n" + "-"*70)
    print("STEP 1: SELECT MDD/CASE DIRECTORY")
    print("-"*70)

    mdd_dir = select_directory_interactive(mdd_dirs, other_dirs, "MDD/Case")

    # Configure Control directory
    print("\n" + "-"*70)
    print("STEP 2: SELECT CONTROL DIRECTORY")
    print("-"*70)

    ctrl_dir = select_directory_interactive(ctrl_dirs, other_dirs, "Control")

    # Collect files
    print("\n" + "-"*70)
    print("COLLECTING FILES")
    print("-"*70)

    mdd_files = collect_tsv_files_from_directory(mdd_dir, "MDD")
    ctrl_files = collect_tsv_files_from_directory(ctrl_dir, "Control")

    print(f"\nFound {len(mdd_files)} MDD files in {mdd_dir}")
    for f in mdd_files[:3]:
        print(f"  • {f.relative_to(mdd_dir)}")
    if len(mdd_files) > 3:
        print(f"  ... and {len(mdd_files)-3} more")

    print(f"\nFound {len(ctrl_files)} Control files in {ctrl_dir}")
    for f in ctrl_files[:3]:
        print(f"  • {f.relative_to(ctrl_dir)}")
    if len(ctrl_files) > 3:
        print(f"  ... and {len(ctrl_files)-3} more")

    # Create unified project structure
    print("\n" + "-"*70)
    print("CREATING UNIFIED PROJECT STRUCTURE")
    print("-"*70)

    # Create a master analysis directory
    master_dir = base_path / "MDD2_Analysis_Master"
    master_dir.mkdir(exist_ok=True)

    # Create symlinks or copies
    unified_tsv_dir = master_dir / "data" / "tsv"
    unified_tsv_dir.mkdir(parents=True, exist_ok=True)

    # Process MDD files
    print(f"\nProcessing MDD files...")
    for src_file in mdd_files:
        # Extract sample ID from filename or directory
        sample_id = extract_sample_id(src_file)
        dest_file = unified_tsv_dir / f"{sample_id}.tsv.gz"

        if not dest_file.exists():
            try:
                dest_file.symlink_to(src_file)
                print(f"  Linked MDD: {sample_id}")
            except:
                import shutil
                shutil.copy2(src_file, dest_file)
                print(f"  Copied MDD: {sample_id}")

    # Process Control files
    print(f"\nProcessing Control files...")
    for src_file in ctrl_files:
        sample_id = extract_sample_id(src_file)
        dest_file = unified_tsv_dir / f"{sample_id}.tsv.gz"

        if not dest_file.exists():
            try:
                dest_file.symlink_to(src_file)
                print(f"  Linked Control: {sample_id}")
            except:
                import shutil
                shutil.copy2(src_file, dest_file)
                print(f"  Copied Control: {sample_id}")

    # Create project mapping
    print(f"\nCreating project mapping...")
    project_map = master_dir / "project_mapping.txt"

    with open(project_map, 'w') as f:
        # Write MDD samples
        for src_file in mdd_files:
            sample_id = extract_sample_id(src_file)
            f.write(f"{sample_id}\tMDD\n")

        # Write Control samples
        for src_file in ctrl_files:
            sample_id = extract_sample_id(src_file)
            f.write(f"{sample_id}\tControl\n")

    print(f"Created project map with {len(mdd_files) + len(ctrl_files)} entries")

    # Show preview
    print(f"\nProject map preview:")
    with open(project_map, 'r') as f:
        for i, line in enumerate(f):
            if i < 5:
                print(f"  {line.strip()}")
            else:
                print(f"  ... and {len(mdd_files) + len(ctrl_files) - 5} more")
                break

    # Create output directory
    output_dir = master_dir / "diff_analysis"
    output_dir.mkdir(exist_ok=True)

    # Return paths
    paths = {
        'base': master_dir,
        'tsv': unified_tsv_dir,
        'output': output_dir,
        'project_map': project_map,
        'mdd_original': mdd_dir,
        'ctrl_original': ctrl_dir,
        'mdd_files': mdd_files,
        'ctrl_files': ctrl_files,
    }

    print(f"\n" + "="*70)
    print("SETUP COMPLETE")
    print("="*70)
    print(f"\nMaster directory: {master_dir}")
    print(f"MDD files: {len(mdd_files)} from {mdd_dir}")
    print(f"Control files: {len(ctrl_files)} from {ctrl_dir}")
    print(f"Total samples: {len(mdd_files) + len(ctrl_files)}")
    print(f"Project map: {project_map}")
    print(f"Output directory: {output_dir}")

    return paths

def select_directory_interactive(preferred_dirs, other_dirs, dir_type):
    """Interactively select a directory"""

    all_options = []

    # Add preferred directories first
    for i, dir_path in enumerate(preferred_dirs, 1):
        all_options.append((i, dir_path, "preferred"))

    # Add other directories
    start_idx = len(preferred_dirs) + 1
    for i, dir_path in enumerate(other_dirs[:10], start_idx):  # Limit to 10
        all_options.append((i, dir_path, "other"))

    # Add custom path option
    custom_idx = len(all_options) + 1
    all_options.append((custom_idx, "Enter custom path", "custom"))

    # Display options
    print(f"\nSelect {dir_type} directory:")
    print("-" * 50)

    for idx, dir_path, dir_type in all_options:
        if dir_path == "Enter custom path":
            print(f"{idx}. {dir_path}")
        else:
            # Count TSV files in this directory
            tsv_count = len(collect_tsv_files_from_directory(dir_path, dir_type))
            print(f"{idx}. {dir_path} ({tsv_count} TSV files)")

    print("-" * 50)

    # Get selection
    while True:
        try:
            choice = int(input(f"\nSelect option [1-{len(all_options)}]: "))
            if 1 <= choice <= len(all_options):
                selected = all_options[choice-1]

                if selected[1] == "Enter custom path":
                    # Get custom path
                    while True:
                        custom_path = input(f"Enter path to {dir_type} directory: ").strip()
                        if not custom_path:
                            print("Please enter a path.")
                            continue

                        custom_path = Path(custom_path)
                        if custom_path.exists():
                            # Check if it has TSV files
                            tsv_files = collect_tsv_files_from_directory(custom_path, dir_type)
                            if tsv_files:
                                print(f"✓ Found {len(tsv_files)} TSV files")
                                return custom_path
                            else:
                                print(f"✗ Directory exists but has no TSV files")
                                use_anyway = input("Use anyway? [y/N]: ").lower()
                                if use_anyway == 'y':
                                    return custom_path
                        else:
                            print(f"✗ Directory does not exist: {custom_path}")

                else:
                    return selected[1]

            else:
                print(f"Please enter a number between 1 and {len(all_options)}")

        except ValueError:
            print("Please enter a valid number.")

class Logger:
    """Simple logger class"""
    def __init__(self):
        self.log_file = None

    def log(self, msg):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        message = f"[{timestamp}] INFO: {msg}"
        print(message, flush=True)
        if self.log_file:
            with open(self.log_file, 'a') as f:
                f.write(message + '\n')

    def warning(self, msg):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        message = f"[{timestamp}] WARNING: {msg}"
        print(message, flush=True)
        if self.log_file:
            with open(self.log_file, 'a') as f:
                f.write(message + '\n')

    def error(self, msg, exit_code=1):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        message = f"[{timestamp}] ERROR: {msg}"
        print(message, flush=True)
        if self.log_file:
            with open(self.log_file, 'a') as f:
                f.write(message + '\n')
        sys.exit(exit_code)

def analyze_tsv_files(files, group_type, logger):
    """Analyze TSV files for a specific group"""
    logger.log(f"Analyzing {group_type} files...")

    sample_snps = {}

    for file_path in files:
        sample_id = extract_sample_id(file_path)
        logger.log(f"  Processing {sample_id}...")

        try:
            # Read and count SNPs
            if file_path.suffix == '.gz':
                with gzip.open(file_path, 'rt') as f:
                    # Skip header
                    next(f)
                    snp_count = sum(1 for line in f)
            else:
                with open(file_path, 'r') as f:
                    next(f)
                    snp_count = sum(1 for line in f)

            sample_snps[sample_id] = snp_count

        except Exception as e:
            logger.warning(f"    Error reading {file_path.name}: {e}")
            sample_snps[sample_id] = 0

    return sample_snps

def create_comprehensive_report(paths, mdd_snps, ctrl_snps, logger):
    """Create a comprehensive analysis report"""

    report_dir = paths['output']

    # 1. Create sample summary
    samples_df = pd.DataFrame({
        'Sample_ID': list(mdd_snps.keys()) + list(ctrl_snps.keys()),
        'Group': ['MDD'] * len(mdd_snps) + ['Control'] * len(ctrl_snps),
        'SNP_Count': list(mdd_snps.values()) + list(ctrl_snps.values()),
        'Source_File': [str(f) for f in paths['mdd_files']] + [str(f) for f in paths['ctrl_files']]
    })

    samples_file = report_dir / "sample_summary.csv"
    samples_df.to_csv(samples_file, index=False)

    # 2. Create analysis configuration
    config = f"""MDD2 Differential SNP Analysis Configuration
{'='*70}

ANALYSIS PARAMETERS
{'='*70}
Threshold: 95%
Minimum samples per group: 3
Top N results: 15

PROJECT STRUCTURE
{'='*70}
Master Directory: {paths['base']}
MDD Source: {paths['mdd_original']}
Control Source: {paths['ctrl_original']}
Output Directory: {paths['output']}

SAMPLE STATISTICS
{'='*70}
Total Samples: {len(mdd_snps) + len(ctrl_snps)}
MDD Samples: {len(mdd_snps)}
Control Samples: {len(ctrl_snps)}

{'='*70}
READY FOR ANALYSIS: {'✅ YES' if len(mdd_snps) >= 3 and len(ctrl_snps) >= 3 else '⚠  NEED MORE SAMPLES'}
{'='*70}
"""

    config_file = report_dir / "analysis_configuration.txt"
    with open(config_file, 'w') as f:
        f.write(config)

    logger.log(f"\nSetup reports created:")
    logger.log(f"  • {config_file}")
    logger.log(f"  • {samples_file}")

def main():
    """Main entry point"""

    print("\n" + "="*70)
    print("MDD2 Differential SNP Analyzer - Multi-Directory Version")
    print("="*70)
    print("\nThis tool handles projects with separate Control and MDD directories.")
    print("It will:")
    print("  1. Discover your project structure")
    print("  2. Help you select Control and MDD directories")
    print("  3. Create a unified analysis setup")
    print("  4. Perform ACTUAL differential SNP analysis")
    print("  5. Generate comprehensive results")
    print("\n" + "="*70)

    # Setup logger
    logger = Logger()

    try:
        # Setup multi-directory analysis
        paths = setup_multi_directory_analysis()

        # Setup logger file
        log_file = paths['output'] / "analysis.log"
        logger.log_file = log_file

        # Analyze TSV files (just for summary)
        mdd_snps = analyze_tsv_files(paths['mdd_files'], "MDD", logger)
        ctrl_snps = analyze_tsv_files(paths['ctrl_files'], "Control", logger)

        # Create initial report
        create_comprehensive_report(paths, mdd_snps, ctrl_snps, logger)

        # Ask user if they want to run analysis
        print("\n" + "="*70)
        print("READY FOR DIFFERENTIAL SNP ANALYSIS")
        print("="*70)

        logger.log(f"\nSample counts:")
        logger.log(f"  MDD samples: {len(mdd_snps)}")
        logger.log(f"  Control samples: {len(ctrl_snps)}")
        logger.log(f"  Total samples: {len(mdd_snps) + len(ctrl_snps)}")

        if len(mdd_snps) < 3 or len(ctrl_snps) < 3:
            logger.warning(f"\nWARNING: Insufficient samples for reliable analysis!")
            logger.warning(f"Minimum recommended: 3 per group")
            logger.warning(f"Current: MDD={len(mdd_snps)}, Control={len(ctrl_snps)}")
            proceed = input("\nProceed anyway? [y/N]: ").lower()
            if proceed != 'y':
                logger.log("Analysis cancelled by user.")
                return
        else:
            proceed = input("\nProceed with differential SNP analysis? [Y/n]: ").lower()
            if proceed == 'n':
                logger.log("Analysis cancelled by user.")
                return

        # Run the actual differential analysis
        logger.log("Starting differential SNP analysis...")

        success = perform_differential_analysis(
            master_dir=paths['base'],
            threshold=95,
            min_samples=3,
            top_n=15,
            logger=logger
        )

        if success:
            print("\n" + "="*70)
            print("ANALYSIS COMPLETE")
            print("="*70)
            print(f"\nDifferential SNP analysis has been completed successfully!")
            print(f"\nResults saved in: {paths['output']}")
            print(f"\nFiles created:")

            for file in sorted(paths['output'].glob("*")):
                if file.is_file():
                    size_kb = file.stat().st_size / 1024
                    print(f"  • {file.name} ({size_kb:.1f} KB)")

            print(f"\nKey results:")
            print(f"  1. all_differential_snps.tsv - All differential SNPs")
            print(f"  2. top_*_specific_snps.tsv - Top specific SNPs")
            print(f"  3. differential_snps_report.xlsx - Excel report")
            print(f"  4. analysis_summary.txt - Text summary")
            print(f"\nLog file: {log_file}")
            print("\n" + "="*70)

    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user.")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
