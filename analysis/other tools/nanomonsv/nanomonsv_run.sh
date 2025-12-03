### ont data:

#! /bin/bash
export PATH=/home/user/miniconda3/envs/nanomonsv/bin:$PATH 
nanomonsv parse {tumorBam} {tumorSaveDir}  
nanomonsv parse {controlBam} {controlSaveDir} 
nanomonsv get {tumorSaveDir} {tumorBam} {reference} --control_prefix {controlSaveDir} --control_bam {controlBam} --single_bnd --use_racon --processes {threads} --control_panel_prefix {nanomonsv_data_path}/hprc_year1_data_freeze_nanopore_guppy4_minimap2_2_24_merge_control_GRCh38


### PacBio data:

#! /bin/bash
export PATH=/home/user/miniconda3/envs/nanomonsv/bin:$PATH 
nanomonsv parse {tumorBam} {tumorSaveDir}  
nanomonsv parse {controlBam} {controlSaveDir} 
nanomonsv get {tumorSaveDir} {tumorBam} {reference} --control_prefix {controlSaveDir} --control_bam {controlBam} --single_bnd --use_racon --processes {threads} --control_panel_prefix {nanomonsv_data_path}/hprc_year1_data_freeze_PacBio_HiFi_minimap2_2_24_merge_control_GRCh38
