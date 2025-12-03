#! /bin/bash

{longphase} haplotag -r {reference} -s {saveDir}/phased_merge_output.vcf.gz -b {tumorBam} -t {threads} -o {saveDir}/tumor_haplotagged --tagSupplementary --qualityThreshold 10 
{samtools} index -@ {threads} {saveDir}/tumor_haplotagged.bam 
{longphase} haplotag -r {reference} -s {saveDir}/phased_merge_output.vcf.gz -b {controlBam} -t {threads} -o {saveDir}/normal_haplotagged --tagSupplementary --qualityThreshold 10 
{samtools} index -@ {threads} {saveDir}/normal_haplotagged.bam
