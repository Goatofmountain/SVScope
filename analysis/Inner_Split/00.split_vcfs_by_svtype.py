#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re 
import time
import subprocess
import shutil
import gzip
import sys
import glob

# --- Configuration ---

# 1. Your project base path
BASE_DIR = "/NAS/wg_lyc/project/TDscope/fanxiu"

# 2. New output directory
OUTPUT_DIR = "/NAS/wg_lyc/project/TDscope/fanxiu/Stat/Split_Inner"

# 3. Cell line and technology directories to iterate (from your list)
CELL_LINE_DIRS = [
    'COLO829_hifi', 'H1437_ont', 'H2009_ont', 'HCC1937_hifi', 'HCC1954.Nano_ont', 'HG008_hifi',
    'COLO829.Nano_ont', 'H2009_hifi', 'HCC1395_hifi', 'HCC1937_ont', 'HCC1954_ont', 'HG008_ont',
    'H1437_hifi', 'H2009.Nano_ont', 'HCC1395_ont', 'HCC1954_hifi'
]

# 4. Caller method names and corresponding file path patterns
CALLER_DATA = [
    # (Method Name, Path Pattern)
    ('nanomonsv', 'nanomonsv/*/*.nanomonsv.result.vcf'),
    ('svision-pro', 'svision-pro/{tech_placeholder}.svision_pro_v1.8.s3.somatic_s1.vcf'),
    ('savana', 'savana/{tech_placeholder}.classified.somatic.vcf'),
    ('severus', 'severus_clair3_out/somatic_SVs/severus_somatic.vcf'),
    ('sniffles2_multiple_supp3', 'sniffles2_multiple_supp3/{tech_placeholder}_somatic_sv.vcf')
]

# --- Helper Functions ---

def get_sv_class(vcf_cols):
    """
    Classifies a VCF line as 'inner' or 'split' based on the INFO and ALT fields.
    
    Inner: DEL, INS, or complex types consisting only of DEL/INS
    Split: DUP, INV, TRA, BND, or any complex types including these
    """
    if len(vcf_cols) < 8:
        return 'split' # Malformed line, classify as split for safety
    
    info_field = vcf_cols[7]
    alt_field = vcf_cols[4]
    
    svtype = None
    info_parts = info_field.split(';')
    for part in info_parts:
        if part.startswith('SVTYPE='):
            svtype = part.replace('SVTYPE=', '').upper()
            break
    
    if svtype:
        # 1. Simple types
        if svtype == 'INS' or svtype == 'DEL':
            return 'inner'
        
        if svtype in ['INV', 'DUP', 'TRA', 'BND']:
            return 'split'
        
        # 2. Complex types (e.g., "DEL,INS" or "DEL,DUP")
        if ',' in svtype:
            components = re.split(r'[,+]', svtype)
            has_split_type = False
            for comp in components:
                if comp in ['INV', 'DUP', 'TRA', 'BND']:
                    has_split_type = True
                    break
            
            if has_split_type:
                return 'split'  # Contains any of INV/DUP/TRA/BND
            else:
                return 'inner' # Consists only of INS and/or DEL

    # 3. Fallback for BNDs (if VCF lacks SVTYPE=BND tag)
    if '[' in alt_field or ']' in alt_field:
        return 'split'

    # 4. All other cases (including unknown SVTYPEs, <CNV>, etc.) default to split
    return 'split'

def process_vcf_file(vcf_path, output_subdir, tool_name):
    """
    Reads a VCF file and splits it into two new files based on 'inner' and 'split' classification.
    """
    out_inner_path = os.path.join(output_subdir, "inner.vcf")
    out_split_path = os.path.join(output_subdir, "split.vcf")
    
    if not os.path.exists(vcf_path):
        print(f"  [WARNING] VCF file not found, skipping: {vcf_path}")
        return False

    # Cache check: Skip if both files already exist
    if os.path.exists(out_inner_path) and os.path.exists(out_split_path):
        print(f"  [Cache] inner/split files for {tool_name} already exist, skipping.")
        return True

    print(f"  > Splitting: {tool_name} ...")
    os.makedirs(output_subdir, exist_ok=True)
    
    try:
        # Use 'latin-1' encoding to avoid UnicodeDecodeError
        if vcf_path.endswith('.gz'):
            fin = gzip.open(vcf_path, 'rt', encoding='latin-1', errors='ignore')
        else:
            fin = open(vcf_path, 'r', encoding='latin-1', errors='ignore')
        
        with open(out_inner_path, 'w') as fout_inner, open(out_split_path, 'w') as fout_split:
            for line in fin:
                if line.startswith('##'):
                    fout_inner.write(line)
                    fout_split.write(line)
                    continue
                
                if line.startswith('#CHROM'):
                    fout_inner.write(line)
                    fout_split.write(line)
                    continue
                
                if not line.strip(): 
                    continue
                    
                cols = line.strip().split('\t')
                if len(cols) < 8: 
                    continue
                
                sv_class = get_sv_class(cols)
                
                if sv_class == 'inner':
                    fout_inner.write(line)
                else: # 'split'
                    fout_split.write(line)
                    
        fin.close()
        return True
        
    except FileNotFoundError:
        print(f"  [!!ERROR!!] File lost during processing: {vcf_path}")
        # Clean up incomplete files
        if os.path.exists(out_inner_path): os.remove(out_inner_path)
        if os.path.exists(out_split_path): os.remove(out_split_path)
        return False
    except Exception as e:
        print(f"  [!!ERROR!!] Failed while processing {vcf_path}: {e}")
        return False


# --- Script Main Flow ---
def main():
    print(f"Script started...")
    print(f"Base Directory: {BASE_DIR}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output Directory: {OUTPUT_DIR}")

    # Stores (vcf_prefix, path) tuples
    benchmarks_to_process = set() 

    # Iterate through each cell line directory
    for SampleID in CELL_LINE_DIRS:
        CaseDir = os.path.join(BASE_DIR, SampleID) 
        
        if not os.path.isdir(CaseDir):
            print(f"\n[WARNING] Sample directory not found, skipping: {CaseDir}")
            continue

        print(f"\n========================================================")
        print(f"Starting check for cell line: {SampleID}")
        print(f"========================================================")

        MethodDict = {} # Stores (MethodName -> FilePath)

        # 1. Parse cell line, technology, and base name
        try:
            if 'hifi' in SampleID:
                tech_short = 'PB'
                cell_line_base = SampleID.replace('_hifi', '')
            elif 'ont' in SampleID:
                tech_short = 'ONT'
                cell_line_base = SampleID.replace('_ont', '')
            else:
                print(f"  [WARNING] Cannot recognize technology (hifi/ont) in {SampleID}, skipping VCF lookup.")
                continue
        except Exception as e:
            print(f"  [ERROR] Error parsing {SampleID} name: {e}")
            continue

        vcf_prefix = cell_line_base
        if '.Nano' in vcf_prefix:
            vcf_prefix = vcf_prefix.split('.')[0] # 'HCC1954.Nano' -> 'HCC1954'
        
        # 2. Iterate through CALLER_DATA (Standard Methods)
        for method_name, path_pattern in CALLER_DATA:
            formatted_pattern = path_pattern.replace('{tech_placeholder}', SampleID)
            full_path_pattern = os.path.join(CaseDir, formatted_pattern)
            
            try:
                found_files = glob.glob(full_path_pattern)
                if found_files:
                    abs_file_path = os.path.abspath(found_files[0])
                    MethodDict[method_name] = abs_file_path
            except Exception as e:
                print(f"  [ERROR] Error searching pattern: {full_path_pattern} ({e})")
        
        # 3. Handle special paths for abPOA and spoa
        vcf_filename = "HG008T.mergedSomatic.vcf" if vcf_prefix == "HG008" else f"{vcf_prefix}.mergedSomatic.vcf"
        
        special_paths_to_check = []
        if tech_short == 'PB':
            special_paths_to_check = [
                (f"{cell_line_base}_abPOA", os.path.join(CaseDir, f"{vcf_prefix}_abPOA_PacBioSpecific", vcf_filename)),
                (f"{cell_line_base}_sPOA" if vcf_prefix == "HG008" else f"{cell_line_base}_spoa", 
                 os.path.join(CaseDir, f"{vcf_prefix}_spoa_PacBioSpecific", vcf_filename))
            ]
        elif tech_short == 'ONT':
            ont_spoa_dir_name = {
                'COLO829': 'COLO829_spoa', 'H1437': 'H1437_sPOA', 'H2009': 'H2009_sPOA',
                'HCC1395': 'HCC1395_sPOA', 'HCC1937': 'HCC1937_sPOA', 'HCC1954': 'HCC1954_sPOA',
                'HG008': 'HG008_spoa'
            }.get(vcf_prefix)

            if '.Nano' in cell_line_base:
                ont_spoa_dir_name = {
                    'COLO829.Nano': 'COLO829_spoa',
                    'H2009.Nano': 'H2009_spoa',
                    'HCC1954.Nano': 'HCC1954_spoa'
                }.get(cell_line_base)
            
            if ont_spoa_dir_name:
                special_paths_to_check = [
                    (f"{cell_line_base}_abPOA", os.path.join(CaseDir, f"{vcf_prefix}_abPOA", vcf_filename)),
                    (f"{cell_line_base}_spoa", os.path.join(CaseDir, ont_spoa_dir_name, vcf_filename))
                ]
        
        for method_name, full_path in special_paths_to_check:
            if os.path.exists(full_path):
                MethodDict[method_name] = os.path.abspath(full_path)
        
        # 4. Collect Benchmark paths
        # Check if benchmark for this vcf_prefix has already been collected
        if vcf_prefix not in [p[0] for p in benchmarks_to_process]:
            benchmark_path = ""
            if vcf_prefix == 'COLO829':
                benchmark_path = '/NAS/wg_tkl/SVScope_Data/Tools/COLO829_Nano/truthset_somaticSVs_COLO829_hg38lifted.chrom.vcf'
            elif vcf_prefix == 'HG008':
                benchmark_path = '/NAS/wg_lyc/project/TDscope/fanxiu/HG008_Golden/ensemble/GRCh38_HG008-T-V0.4_somatic-stvar_PASS.draftbenchmark.vcf'
            else:
                benchmark_path = f'/NAS/wg_tkl/SVScope_Data/GoldenSet/Filter_2Tech5Tool/{vcf_prefix}_minda_ensemble.Filter.vcf'
            
            if os.path.exists(benchmark_path):
                # Add (base_name, path) to the set
                benchmarks_to_process.add((vcf_prefix, benchmark_path))
            else:
                print(f"  [!!WARNING!!] Benchmark file not found: {benchmark_path}")

        # 5. --- Loop through all found *Tool* VCFs ---
        print(f"  > Found {len(MethodDict)} *Tool* VCF files for {SampleID}. Starting processing...")
        
        if not MethodDict:
            print(f"  [WARNING] No tool VCF files found for {SampleID}.")
            
        for tool_name, vcf_path in MethodDict.items():
            # Create output subdirectory for this file
            output_subdir = os.path.join(OUTPUT_DIR, SampleID, tool_name)
            process_vcf_file(vcf_path, output_subdir, tool_name)
            
    # --- Process Benchmark Files Separately ---
    print(f"\n========================================================")
    print(f"Starting to process {len(benchmarks_to_process)} unique benchmark files")
    print(f"========================================================")
    
    for vcf_prefix, benchmark_path in benchmarks_to_process:
        # Define new output path: .../Split_Inner/benchmark/COLO829/
        output_subdir = os.path.join(OUTPUT_DIR, "benchmark", vcf_prefix)
        
        # Define tool name (for logging and caching)
        tool_name = f"{vcf_prefix}_benchmark"
        
        # Call the same processing function
        process_vcf_file(benchmark_path, output_subdir, tool_name)

    print("\n--- Script execution complete ---")
    print(f"All VCF splitting tasks are complete.")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    if not os.path.isdir(BASE_DIR):
        print(f"Error: Base directory (BASE_DIR) does not exist: {BASE_DIR}")
        print("Please modify the BASE_DIR variable in the script.")
        sys.exit(1)
        
    main()
