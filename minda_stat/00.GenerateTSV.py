#!/usr/bin/env python3

import os
import glob
import sys

# 1. Your project base path
BASE_DIR = "/NAS/wg_lyc/project/TDscope/fanxiu"

# 2. Cell line and technology directories to iterate (from your list)
CELL_LINE_DIRS = [
    'COLO829_hifi', 'H1437_ont', 'H2009_ont', 'HCC1937_hifi', 'HCC1954.Nano_ont', 'HG008_hifi',
    'COLO829.Nano_ont', 'H2009_hifi', 'HCC1395_hifi', 'HCC1937_ont', 'HCC1954_ont', 'HG008_ont',
    'H1437_hifi', 'H2009.Nano_ont', 'HCC1395_ont', 'HCC1954_hifi'
]

# 3. Caller method names and corresponding file path patterns
CALLER_DATA = [
    # (Method Name, Path Pattern)
    ('nanomonsv', 'nanomonsv/*/*.nanomonsv.result.vcf'),
    ('svision-pro', 'svision-pro/{tech_placeholder}.svision_pro_v1.8.s3.somatic_s1.vcf'),
    ('savana', 'savana/{tech_placeholder}.classified.somatic.vcf'),
    ('severus_new', 'severus_clair3_out/somatic_SVs/severus_somatic.vcf'),
    ('sniffles2_multiple_supp3', 'sniffles2_multiple_supp3/{tech_placeholder}_somatic_sv.vcf')
]

# --- Script Main Body ---

def main():
    print(f"Script started...")
    print(f"Base Directory: {BASE_DIR}")
    # Define output directory
    output_dir = os.path.join(BASE_DIR, 'Stat', 'Tools')
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output Directory: {output_dir}")

    total_tsvs_generated = 0
    total_files_found = 0

    # Iterate through each cell line directory
    for cell_line_full in CELL_LINE_DIRS:
        
        if cell_line_full in EXCLUDE_DIRS:
            print(f"\nSkipping (Excluded): {cell_line_full}")
            continue

        print(f"\n--- Processing: {cell_line_full} ---")

        # 1. Parse cell line, technology, and base name
        try:
            if 'hifi' in cell_line_full:
                tech_short = 'PB'
                cell_line_base = cell_line_full.replace('_hifi', '')
            elif 'ont' in cell_line_full:
                tech_short = 'ONT'
                cell_line_base = cell_line_full.replace('_ont', '') # e.g., 'HCC1954.Nano' or 'H1437'
            else:
                print(f"  [WARNING] Cannot identify technology (hifi/ont) in {cell_line_full}, skipping.")
                continue
        except Exception as e:
            print(f"  [ERROR] Error parsing {cell_line_full} name: {e}")
            continue

        # 2. Prepare to write new file for this cell line
        rows_for_this_file = []
        output_filename = os.path.join(output_dir, f"{cell_line_full}.tsv")
        files_found_for_this_cell_line = 0

        # 3. Iterate through all CALLER_DATA (Standard Methods)
        for method_name, path_pattern in CALLER_DATA:
            formatted_pattern = path_pattern.replace('{tech_placeholder}', cell_line_full)
            full_path_pattern = os.path.join(BASE_DIR, cell_line_full, formatted_pattern)
            
            try:
                found_files = glob.glob(full_path_pattern)
            except Exception as e:
                print(f"  [ERROR] Error searching pattern: {full_path_pattern} ({e})")
                continue
                
            if found_files:
                sample_method = f"{cell_line_base}_{method_name}"
                for file_path in found_files:
                    abs_file_path = os.path.abspath(file_path)
                    row = f"{abs_file_path}\t{sample_method}\t{tech_short}"
                    rows_for_this_file.append(row)
                    files_found_for_this_cell_line += 1
                    total_files_found += 1
                    print(f"  [FOUND] {method_name} -> {abs_file_path}")

        # 4. Handle special paths for abPOA and spoa
        
        # vcf_prefix is the true base name used to construct paths and filenames (e.g., 'HCC1954')
        vcf_prefix = cell_line_base
        if '.Nano' in vcf_prefix:
            vcf_prefix = vcf_prefix.split('.')[0] # 'HCC1954.Nano' -> 'HCC1954'
        
        # vcf_filename is now based on vcf_prefix
        vcf_filename = "HG008T.mergedSomatic.vcf" if vcf_prefix == "HG008" else f"{vcf_prefix}.mergedSomatic.vcf"
        
        special_paths_to_check = []
        
        if tech_short == 'PB':
            # PacBio / hifi paths
            # Use vcf_prefix to construct subdirectory
            special_paths_to_check = [
                # Method name uses cell_line_base (e.g., HCC1954_abPOA), path uses vcf_prefix
                (f"{cell_line_base}_abPOA", f"{cell_line_full}/{vcf_prefix}_abPOA_PacBioSpecific/{vcf_filename}"),
                (f"{cell_line_base}_sPOA" if vcf_prefix == "HG008" else f"{cell_line_base}_spoa", 
                 f"{cell_line_full}/{vcf_prefix}_spoa_PacBioSpecific/{vcf_filename}")
            ]
        elif tech_short == 'ONT':
            # Use vcf_prefix as the map key
            ont_spoa_dir_name = {
                'COLO829': 'COLO829_spoa',
                'H1437': 'H1437_sPOA',
                'H2009': 'H2009_sPOA', # Applies to H2009_ont
                'HCC1395': 'HCC1395_sPOA',
                'HCC1937': 'HCC1937_sPOA',
                'HCC1954': 'HCC1954_sPOA', # Applies to HCC1954_ont
                'HG008': 'HG008_spoa'
            }.get(vcf_prefix)

            # Add separate mapping for .Nano samples (e.g., HCC1954.Nano_ont)
            if '.Nano' in cell_line_base:
                ont_spoa_dir_name = {
                    'COLO829.Nano': 'COLO829_spoa',
                    'H2009.Nano': 'H2009_spoa',
                    'HCC1954.Nano': 'HCC1954_spoa' 
                }.get(cell_line_base)

            if ont_spoa_dir_name:
                special_paths_to_check = [
                    # Method name uses cell_line_base (e.g., HCC1954.Nano_abPOA), path uses vcf_prefix
                    (f"{cell_line_base}_abPOA", f"{cell_line_full}/{vcf_prefix}_abPOA/{vcf_filename}"),
                    # Method name uses cell_line_base (e.g., HCC1954.Nano_spoa)
                    (f"{cell_line_base}_spoa", f"{cell_line_full}/{ont_spoa_dir_name}/{vcf_filename}") 
                ]
            else:
                 print(f"  [WARNING] Could not find ONT spoa mapping for {cell_line_base} (vcf_prefix: {vcf_prefix}), skipping special paths.")
        
        # General check logic (for special paths)
        for method_name, path_suffix in special_paths_to_check:
            full_path = os.path.join(BASE_DIR, path_suffix)
            
            if os.path.exists(full_path):
                abs_file_path = os.path.abspath(full_path)
                row = f"{abs_file_path}\t{method_name}\t{tech_short}"
                rows_for_this_file.append(row)
                files_found_for_this_cell_line += 1
                total_files_found += 1
                print(f"  [FOUND] {method_name} (Special) -> {abs_file_path}")
            else:
                # Print paths not found to aid debugging
                # print(f"  -> Not found (Special): {full_path}")
                pass

        # 5. Write TSV file
        if len(rows_for_this_file) > 1: # Check if more than one row of data was found
            try:
                with open(output_filename, 'w', encoding='utf-8') as f:
                    f.write('\n'.join(rows_for_this_file))
                print(f"  [SUCCESS] Generated {output_filename} (containing {files_found_for_this_cell_line} entries)")
                total_tsvs_generated += 1
            except Exception as e:
                print(f"  [!!ERROR!!] Could not write file {output_filename}: {e}")
        else:
            print(f"  [NOTE] No VCF files found for {cell_line_full}, no TSV generated.")

    print("\n--- Script execution complete ---")
    print(f"Total {total_tsvs_generated} TSV files generated.")
    print(f"Total {total_files_found} VCF file paths found.")

if __name__ == "__main__":
    if not os.path.isdir(BASE_DIR):
        print(f"Error: Base directory (BASE_DIR) does not exist: {BASE_DIR}")
        print("Please modify the BASE_DIR variable in the script.")
        sys.exit(1)
        
    main()
