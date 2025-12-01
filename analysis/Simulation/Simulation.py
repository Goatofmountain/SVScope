#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re 
import sys
import glob
import time
import shutil
import gzip
import subprocess
import random
import pandas as pd
from multiprocessing import Pool

# --- Configuration ---

# Genome Reference Files
HG38_FA = "/NAS/wg_fzt/benchmark/simu_somatic/data/hg38_mainChr.fa"
CHR22_FA = '/NAS/wg_ylf/wg_ylf/mainChr_hg38/chr22.fa' # Used for Badread batch0

# Directory and Tool Paths
SIM_REPEAT_DIR = '/NAS/wg_zql/PanCancer_zql/RepeatMaskerOut/Simulation_RepeatOut_7_chr22'
WINDOW_DIR = '/NAS/wg_zql/PanCancer_zql/RepeatMaskerOut/Simulation/WindowList'
MINIMAP2_PATH = "/NAS/wg_liuxy/04.software/minimap2-2.26/minimap2/minimap2"
SAMTOOLS_PATH = "/NAS/wg_fzt/software/samtools-1.9/samtools1.9/bin/samtools"
BADREAD_PATH = "/NAS/wg_ylf/software/Badread-main/badread-runner.py"

# Simulation Parameters
DIRS_LIST = ['Control', 'Case_1', 'Case_2', 'Case_3']
PURITY_LIST = ['purity10', 'purity20', 'purity30', 'purity40', 'purity50']
BADREAD_COMMON_ARGS = '--glitches 1000,20,20 --junk_reads 0.1 --random_reads 0.1 --chimeras 0.1 --identity 80,90,6 --length 20000,5000'

# Purity and Coverage Mappings
# Batch 0/1 (Germline/Reference) Coverage
GERMLINE_COV = {
    'purity10': '18x', 'purity20': '16x', 'purity30': '14x', 
    'purity40': '12x', 'purity50': '10x'
}
# Batch 2 (Somatic) Coverage
SOMATIC_COV = {
    'purity10': '4x', 'purity20': '8x', 'purity30': '12x', 
    'purity40': '16x', 'purity50': '20x'
}

# Sequence Multiplication Factors for Batch 1 (Germline) and Batch 2 (Somatic)
# (Case_Name: [Batch1_Mult, Batch2_Mult]) -> Note: Control uses 0x for Batch 2
TD_MULTIPLIERS = {
    'Control': [5, 0],
    'Case_1': [5, 10],
    'Case_2': [5, 15],
    'Case_3': [5, 20],
}

# Window file mapping
FILE_MAPPING = {'_1_': 'seq_A', '_2_': 'seq_B', '_3_': 'seq_C', '_4_': 'seq_D'}
WINDOW_DICT = {}
for filename in os.listdir(WINDOW_DIR) if os.path.exists(WINDOW_DIR) else []:
    for key, value in FILE_MAPPING.items():
        if key in filename:
            WINDOW_DICT[value] = os.path.join(WINDOW_DIR, filename)

# --- Helper Functions ---

def generate_hack_files(case_name, mult_B1, mult_B2):
    """
    Consolidates the logic for generating .NoneMut.bed and .HACK.bed files 
    for Control, Case_1, Case_2, and Case_3.
    
    mult_B1: Multiplier for Batch 1 (seq_B, seq_D)
    mult_B2: Multiplier for Batch 2 (seq_C, seq_D)
    """
    print(f"--- Generating HACK files for {case_name} (B1:{mult_B1}x, B2:{mult_B2}x) ---")

    # Define paths for input sequences and output files
    file_B_path = WINDOW_DICT.get('seq_B')
    file_C_path = WINDOW_DICT.get('seq_C')
    file_D_path = WINDOW_DICT.get('seq_D')

    if not all([file_B_path, file_C_path, file_D_path]):
        print("[ERROR] Required sequence files not found in WINDOW_DICT, skipping HACK generation.")
        return

    # 1. Generate .NoneMut.bed (Temporary files with inserted sequence)
    output1_mut = os.path.join(SIM_REPEAT_DIR, f'{case_name}_batch1.NoneMut.bed')
    output2_mut = os.path.join(SIM_REPEAT_DIR, f'{case_name}_batch2.NoneMut.bed')

    with open(file_B_path, 'r') as file_B, \
         open(file_C_path, 'r') as file_C, \
         open(file_D_path, 'r') as file_D, \
         open(output1_mut, 'w') as output_file1, \
         open(output2_mut, 'w') as output_file2:

        # Batch 1 sequences (seq_B * mult_B1, seq_D * mult_B1)
        for line in file_B.readlines():
            columns = line.strip().split('\t')
            TDseq = columns[1][0:]
            output_file1.write(f"{columns[0]}\t{TDseq * mult_B1}\n")
            
        for line in file_D.readlines():
            columns = line.strip().split('\t')
            TDseq = columns[1][0:]
            output_file1.write(f"{columns[0]}\t{TDseq * mult_B1}\n")

        # Batch 2 sequences (seq_C * mult_B2, seq_D * mult_B2)
        for line in file_C.readlines():
            columns = line.strip().split('\t')
            TDseq = columns[1][0:]
            output_file2.write(f"{columns[0]}\t{TDseq * mult_B2}\n")

        # seq_D is shared between batches, but must be added to Batch 2 here
        file_D.seek(0) # Reset pointer for file_D
        for line in file_D.readlines():
            columns = line.strip().split('\t')
            TDseq = columns[1][0:]
            output_file2.write(f"{columns[0]}\t{TDseq * mult_B2}\n")

    # 2. Generate .HACK.bed (VISOR HACk format)
    input1_mut = output1_mut
    input2_mut = output2_mut
    output1_hack = os.path.join(SIM_REPEAT_DIR, f'{case_name}_batch1.HACK.bed')
    output2_hack = os.path.join(SIM_REPEAT_DIR, f'{case_name}_batch2.HACK.bed')
    
    def process_mut_file(input_file, output_file):
        """Converts .NoneMut.bed to .HACK.bed format."""
        if os.path.getsize(input_file) == 0:
            return
            
        with open(input_file, 'r') as batch_in, open(output_file, 'w') as batch_out:
            for line in batch_in.readlines():
                columns = line.strip().split('\t')
                chrom_start_end = columns[0].split('_')
                Chrom = chrom_start_end[0]
                # Start = int(chrom_start_end[1]) # Unused
                End = int(chrom_start_end[2])
                TDseq = columns[1][0:]
                # Format: Chrom Start End SV_Type TD_Sequence 0
                batch_out.write(f"{Chrom}\t{End}\t{End}\t{'insertion'}\t{TDseq}\t{0}\n")
        
        # Clean up intermediate file
        os.remove(input_file) 
    
    process_mut_file(input1_mut, output1_hack)
    process_mut_file(input2_mut, output2_hack)
    print(f"--- HACK files generated successfully ---")


def run_visor_hack(dir_name, purity, batch):
    """
    Runs the VISOR HACk tool to generate mutated FASTA files (h1.fa)
    based on the .HACK.bed files.
    """
    PURITY_dir = os.path.join(SIM_REPEAT_DIR, purity)
    HACKBED_dir = os.path.join(PURITY_dir, 'HACKbed')
    
    # Create necessary output directory structure
    hack_out_dir = os.path.join(HACKBED_dir, dir_name)
    os.makedirs(hack_out_dir, exist_ok=True)
    
    HACK_bed = os.path.join(SIM_REPEAT_DIR, f'{dir_name}_{batch}.HACK.bed')
    hack_merge_out = os.path.join(HACKBED_dir, dir_name, f'{dir_name}_{batch}_HACK')
    
    if not os.path.exists(HACK_bed):
        # This will happen for Control_batch2 if mult_B2 was 0, resulting in empty files
        print(f"[WARNING] HACK bed not found or empty for {dir_name}_{batch}, skipping VISOR HACk.")
        return

    # Check if h1.fa (output) already exists before running (simple cache)
    output_fa = os.path.join(hack_merge_out, 'h1.fa')
    if os.path.exists(output_fa):
        print(f"[Cache] VISOR HACk output for {dir_name}_{batch} already exists.")
        return

    print(f"  > Running VISOR HACk for {dir_name}/{batch}...")
    
    # VISOR HACk command execution
    visor_cmd = 'VISOR HACk -g {genome_FASTA} -b {HACk_bed} -o {hack_out}'.format(
        genome_FASTA=CHR22_FA, HACK_bed=HACK_bed, hack_out=hack_merge_out)
    
    # os.system(visor_cmd) # Kept commented as per original intent, but needed for simulation
    print(f"    [SIMULATED EXECUTION] {visor_cmd}")
    # Placeholder: Assuming h1.fa is generated in hack_merge_out/
    # If this step is run, uncomment os.system(visor_cmd)


def run_badread_simulation(dir_name, purity, batch, quantity):
    """
    Runs Badread simulation for a specific purity, case, and batch.
    """
    PURITY_dir = os.path.join(SIM_REPEAT_DIR, purity)
    HACKBED_dir = os.path.join(PURITY_dir, 'HACKbed')
    
    # Determine input FASTA based on batch
    if batch == 'batch0':
        hack_merge_fa = CHR22_FA # Reference genome
        out_subdir = os.path.join(HACKBED_dir, dir_name, f'{dir_name}_batch0_HACK')
    else:
        # Assumes h1.fa was generated by VISOR HACk in the previous step
        out_subdir = os.path.join(HACKBED_dir, dir_name, f'{dir_name}_{batch}_HACK')
        hack_merge_fa = os.path.join(out_subdir, 'h1.fa') 
    
    # Define output file path
    simulate_fq = os.path.join(out_subdir, 'reads.fastq.gz')
    
    os.makedirs(out_subdir, exist_ok=True)

    # Check for output file before running
    if os.path.exists(simulate_fq):
        print(f"  [Cache] Badread output for {dir_name}/{batch} ({purity}) already exists.")
        return

    print(f"  > Simulating Badread for {dir_name}/{batch} ({purity}, {quantity})...")
    
    badread_cmd = f''' {BADREAD_PATH} simulate --reference {hack_merge_fa} --quantity {quantity} {BADREAD_COMMON_ARGS} | gzip > {simulate_fq} '''
    
    # os.system(badread_cmd)
    print(f"    [SIMULATED EXECUTION] Badread command.")


def concatenate_fastq(dir_name, purity):
    """
    Concatenates batch0, batch1, and batch2 FastQ files.
    """
    PURITY_dir = os.path.join(SIM_REPEAT_DIR, purity)
    HACKBED_dir = os.path.join(PURITY_dir, 'HACKbed')
    
    batch0_HACK = os.path.join(HACKBED_dir, dir_name, f'{dir_name}_batch0_HACK', 'reads.fastq.gz')
    batch1_HACK = os.path.join(HACKBED_dir, dir_name, f'{dir_name}_batch1_HACK', 'reads.fastq.gz')
    batch2_HACK = os.path.join(HACKBED_dir, dir_name, f'{dir_name}_batch2_HACK', 'reads.fastq.gz')
    cat_fq_gz = os.path.join(HACKBED_dir, dir_name, f'{dir_name}.fq.gz')

    if os.path.exists(cat_fq_gz):
        print(f"  [Cache] Concatenated FQ for {dir_name}/{purity} already exists.")
        return
        
    print(f"  > Concatenating FastQ for {dir_name}/{purity}...")
    
    cat_cmd = 'cat {batch0} {batch1} {batch2} > {out_fq}'.format(
        batch0=batch0_HACK, batch1=batch1_HACK, batch2=batch2_HACK, out_fq=cat_fq_gz
    )
    # os.system(cat_cmd)
    print(f"    [SIMULATED EXECUTION] {cat_cmd}")


def minimap2_run(ref, fq, sam, bam, bam_sort):
    """
    Runs minimap2, samtools view, samtools sort, and samtools index.
    """
    if os.path.exists(f"{bam_sort}.bai"):
        print(f"  [Cache] Alignment index for {os.path.basename(bam_sort)} already exists.")
        return

    print(f"  > Running Minimap2 alignment for {os.path.basename(fq)}...")
    
    try:
        # 1. Minimap2 alignment (FQ to SAM)
        cmd1 = f"{MINIMAP2_PATH} --MD -ax map-ont -L -t 50 {ref} {fq} > {sam}"
        subprocess.run(cmd1, shell=True, check=True)
        
        # 2. Samtools view (SAM to BAM)
        cmd2 = f"{SAMTOOLS_PATH} view -@ 50 -bS {sam} -o {bam}"
        subprocess.run(cmd2, shell=True, check=True)
        
        # 3. Samtools sort (BAM to Sorted BAM)
        cmd3 = f"{SAMTOOLS_PATH} sort -@ 50 {bam} -o {bam_sort}"
        subprocess.run(cmd3, shell=True, check=True)
        
        # 4. Samtools index (Sorted BAM to BAI)
        cmd4 = f"{SAMTOOLS_PATH} index {bam_sort}"
        subprocess.run(cmd4, shell=True, check=True)
        
        # Cleanup
        os.remove(sam)
        os.remove(bam)
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Minimap2/Samtools command failed: {e}")
        # Add basic cleanup in case of failure
        if os.path.exists(sam): os.remove(sam)
        if os.path.exists(bam): os.remove(bam)
        if os.path.exists(bam_sort): os.remove(bam_sort)
    except FileNotFoundError:
        print("[ERROR] Minimap2 or Samtools executable not found. Check paths.")


# --- Main Orchestration ---

def main():
    if not os.path.exists(SIM_REPEAT_DIR):
        os.makedirs(SIM_REPEAT_DIR)
    os.chdir(SIM_REPEAT_DIR)
    
    print(f"Script starting in: {os.getcwd()}")
    print("---------------------------------------------")

    # PHASE 1: Generate HACK files
    for dir_name, multipliers in TD_MULTIPLIERS.items():
        generate_hack_files(dir_name, multipliers[0], multipliers[1])
    
    # PHASE 2: Run VISOR HACk (FASTA generation)
    # Note: This step is necessary to generate h1.fa files for Badread batch 1/2.
    print("\n[NOTE] Starting VISOR HACk (FASTA Generation). Actual execution is commented out.")
    pool_visor = Pool(10)
    for purity in PURITY_LIST:
        for dir_name in DIRS_LIST:
            # Control_batch2 has 0x multiplication, might not generate mutations/HACK bed
            for batch in ['batch1', 'batch2']:
                 pool_visor.apply_async(run_visor_hack, args=(dir_name, purity, batch))
    pool_visor.close()
    pool_visor.join()
    
    # PHASE 3: Run Badread Simulation
    print("\n--- Starting Badread Simulation (Multi-process) ---")
    pool_badread = Pool(10)

    # Batch 0 (Reference)
    for purity in PURITY_LIST:
        for dir_name in DIRS_LIST:
            quantity = GERMLINE_COV[purity]
            pool_badread.apply_async(run_badread_simulation, args=(dir_name, purity, 'batch0', quantity))

    # Batch 1 (Germline)
    for purity in PURITY_LIST:
        for dir_name in DIRS_LIST:
            quantity = GERMLINE_COV[purity]
            pool_badread.apply_async(run_badread_simulation, args=(dir_name, purity, 'batch1', quantity))

    # Batch 2 (Somatic)
    for purity in PURITY_LIST:
        for dir_name in DIRS_LIST:
            quantity = SOMATIC_COV[purity]
            pool_badread.apply_async(run_badread_simulation, args=(dir_name, purity, 'batch2', quantity))

    pool_badread.close()
    pool_badread.join()
    
    # PHASE 4: Concatenate FastQ Files (Synchronous or smaller pool)
    print("\n--- Concatenating FastQ Files ---")
    pool_cat = Pool(6)
    for purity in PURITY_LIST:
        for dir_name in DIRS_LIST:
            pool_cat.apply_async(concatenate_fastq, args=(dir_name, purity))
    pool_cat.close()
    pool_cat.join()

    # PHASE 5: Run Minimap2 Alignment
    print("\n--- Running Minimap2 Alignment (Multi-process) ---")
    pool_minimap2 = Pool(6)
    
    for purity in PURITY_LIST:
        for dir_name in DIRS_LIST:
            PURITY_dir = os.path.join(SIM_REPEAT_DIR, purity)
            HACKBED_dir = os.path.join(PURITY_dir, 'HACKbed')
            work_dir = HACKBED_dir
            
            sim_fq = os.path.join(work_dir, dir_name, f'{dir_name}.fq.gz')
            sam = os.path.join(work_dir, dir_name, f'{dir_name}_minimap2_2.sam')
            bam = os.path.join(work_dir, dir_name, f'{dir_name}_minimap2_2.bam')
            bam_sort = os.path.join(work_dir, dir_name, f'{dir_name}_minimap2_2_sort.bam')
            
            pool_minimap2.apply_async(minimap2_run, args=(HG38_FA, sim_fq, sam, bam, bam_sort))

    pool_minimap2.close()
    pool_minimap2.join()

    print("\n--- Full Simulation Pipeline Execution Complete ---")


if __name__ == "__main__":
    if not os.path.isdir(WINDOW_DIR):
        print(f"Error: Window directory (WINDOW_DIR) not found: {WINDOW_DIR}")
        sys.exit(1)
        
    main()