#! /bin/bash
export PATH=/home/user/anaconda3/envs/savana_env/bin:$PATH 
savana --tumour {tumorBam} --normal {controlBam} --ref {reference} --outdir {saveDir} --length 50 --mapq 5 --buffer 10 --cn_binsize 10 --chunksize 1000000 --threads {threads} --sample {sampleID}
