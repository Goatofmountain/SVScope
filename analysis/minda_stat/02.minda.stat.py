#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
import sys

# 1. Root directory for Minda statistics
MINDA_STAT_DIR = '/NAS/wg_lyc/project/TDscope/fanxiu/Stat/Tools'
# 2. Final summary output filename
OUTPUT_FILE = os.path.join(MINDA_STAT_DIR, 'minda_all_samples_summary.tsv')
# 3. Minda result filename to search for
SUMMARY_FILENAME = 'None_minda_results.txt'

# --- Script Main Flow ---
if __name__ == '__main__':
    # This will match */SUMMARY_FILENAME
    search_pattern = os.path.join(MINDA_STAT_DIR, '*', SUMMARY_FILENAME)
    # Find all matching summary files
    summary_files = glob.glob(search_pattern)
    
    if not summary_files:
        print(f"[ERROR] No '{SUMMARY_FILENAME}' files found in '{search_pattern}'.")
        print(f"Please check if MINDA_STAT_DIR path is correct.")
        sys.exit(1)
        
    print(f"Found {len(summary_files)} {SUMMARY_FILENAME} files, starting parsing...")
    
    all_results = []
    
    # Define header
    HEADER = [
        "Caller", 
        "True Positives", 
        "False Negatives", 
        "False Positives", 
        "Precision", 
        "Recall", 
        "F1 Score", 
        "CellLine"
    ]
    all_results.append(HEADER)
    
    # Loop through each found file
    for summary_file in summary_files:
        try:
            # Extract cell line name from the file path
            # Path structure example: .../{CellLine_Tech}_minda/None_minda_results.txt
            cell_line_dir = os.path.dirname(summary_file)
            cell_line = os.path.basename(cell_line_dir)
            
            print(f"  > Processing: {cell_line}")
            
            in_overall_section = False
            with open(summary_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                        
                    # Check if entering the OVERALL section
                    if line.startswith('OVERALL'):
                        in_overall_section = True
                        continue
                        
                    # Stop parsing if we reach the next section
                    if line.startswith('SV TYPE RESULTS'):
                        in_overall_section = False
                        break # Finished with this file
                        
                    # Process data if we are in the OVERALL section
                    if in_overall_section:
                        parts = line.split()
                        if not parts:
                            continue
                            
                        caller = parts[0]
                        
                        # Skip header rows and the 'base' benchmark row
                        if caller == 'True' or caller == 'base':
                            continue
                            
                        # Ensure the row has complete data (at least 7 columns: Caller + 6 metrics)
                        if len(parts) >= 7:
                            tp = parts[1]
                            fn = parts[2]
                            fp = parts[3]
                            precision = parts[4]
                            recall = parts[5]
                            f1 = parts[6]
                            
                            # Add to results list
                            all_results.append([
                                caller, tp, fn, fp, precision, recall, f1, cell_line
                            ])
                            
        except IOError as e:
            print(f"  [WARNING] Could not read file {summary_file}: {e}")
        except Exception as e:
            print(f"  [ERROR] Error parsing file {summary_file}: {e}")
            
    # --- Write Summary File ---
    if len(all_results) > 1: # Ensure there is data besides the header
        try:
            with open(OUTPUT_FILE, 'w') as f_out:
                for row in all_results:
                    # Use Tab Separated Values (\t)
                    f_out.write("\t".join(row) + "\n")
                    
            print(f"\n--- Processing Complete ---")
            print(f"Successfully summarized {len(all_results) - 1} records.")
            print(f"Summary file saved to: {OUTPUT_FILE}")
            
        except IOError as e:
            print(f"\n[CRITICAL ERROR] Could not write the final summary file: {e}")
            
    else:
        print(f"\n--- Processing Complete ---")
        print("[WARNING] Found no data to summarize.")
