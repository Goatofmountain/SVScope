#!/usr/bin/env python3

import os
import re
import time
import sys

# 1. Path to the Minda executable
minda = '/NAS/wg_zql/SoftWare/minda-main/minda.py'

# 2. Working directory
work_dir = '/NAS/wg_lyc/project/TDscope/fanxiu/Stat/Tools'

print(f"Working Directory: {work_dir}")
print(f"Minda Path: {minda}")
print("\n--- Searching for all (CellLine_Tech).tsv files... ---")

try:
    all_files_in_dir = os.listdir(work_dir)
    # Filtering:
    # 1. Must end with .tsv
    # 2. Must not be _input.tsv (that is the merged output from the previous script)
    tsvs_to_process = [
        f for f in all_files_in_dir 
        if f.endswith('.tsv') and not f.endswith('_input.tsv')
    ]
    
    if not tsvs_to_process:
        print(f"Warning: No (CellLine_Tech).tsv files found in {work_dir}.")
        sys.exit(1)
        
except FileNotFoundError:
    print(f"Error: Working directory not found: {work_dir}")
    sys.exit(1)

print(f"Found {len(tsvs_to_process)} TSV files ready for individual processing.")


# --- Part Two: Loop and Run Minda ---

print(f"\n--- Starting Minda runs individually... ---")

# Loop through each TSV file we found
for tsv_filename in tsvs_to_process:
    
    # 1. --- Derive Names ---
    
    # SampleID is the filename without .tsv (e.g., "COLO829.Nano_ont" or "H1437_hifi")
    SampleID = tsv_filename.replace('.tsv', '')
    
    # CellName is the base cell line name (e.g., "COLO829" or "H1437")
    # We split using regex at the first . or _ and take the first part
    CellName = re.split(r'[._]', SampleID)[0]
    
    print(f"\n--- Processing Minda for: {SampleID} (Base Name: {CellName}) ---")

    # 2. --- Define Paths ---
    
    # Benchmark path
    if CellName == 'COLO829':
        benchmark = '/NAS/wg_tkl/SVScope_Data/Tools/COLO829_Nano/truthset_somaticSVs_COLO829_hg38lifted.chrom.vcf'
    elif CellName == 'HG008':
        benchmark = '/NAS/wg_lyc/project/TDscope/fanxiu/HG008_Golden/ensemble/GRCh38_HG008-T-V0.4_somatic-stvar_PASS.draftbenchmark.vcf'
    else:
        benchmark = f'/NAS/wg_tkl/SVScope_Data/GoldenSet/Filter_2Tech5Tool/{CellName}_minda_ensemble.Filter.vcf'
    
    # Minda input TSV (the file we are currently looping through)
    TSV_file = os.path.join(work_dir, tsv_filename)
    
    # Minda output directory (as requested: CellLine_Tech_minda)
    OUT = os.path.join(work_dir, f'{SampleID}_minda')

    # 3. --- Check Files ---
    if not os.path.exists(benchmark):
        print(f"  [!!ERROR!!] Benchmark file not found, skipping: {benchmark}")
        continue
    # TSV_file definitely exists, as we started the loop with it
    
    print(f"  Benchmark: {benchmark}")
    print(f"  TSV Input: {TSV_file}")
    print(f"  Output Dir: {OUT}")
    
    # 4. --- Execute Minda ---
    minda_cmd = f'{minda} truthset --base {benchmark} --tsv {TSV_file} --out_dir {OUT} --min_size 50 --multimatch'
    print(f"  Executing Minda command...")
    
    start_time = time.time()
    os.system(minda_cmd)
    os.system(f'chmod -R 775 {OUT}')
    end_time = time.time()
    
    print(f"  [SUCCESS] Minda run complete. Time elapsed: {end_time - start_time:.2f} seconds.")

print("\n--- Script execution complete ---")
