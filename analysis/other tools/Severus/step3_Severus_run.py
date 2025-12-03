#! /bin/bash
export PATH=/home/user/anaconda3/envs/severus-1.6/bin:$PATH 
{Severus_path}/severus.py --target-bam {saveDir}/tumor_haplotagged.bam --control-bam {saveDir}/normal_haplotagged.bam --vntr-bed {vntr.bed} --out-dir {saveDir} -t {threads} --output-read-ids --phasing-vcf {saveDir}/phased_merge_output.vcf.gz --min-support 3  --min-mapq 10
