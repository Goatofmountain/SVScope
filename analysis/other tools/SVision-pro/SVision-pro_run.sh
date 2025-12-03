### ont data:

#! /bin/bash
export PATH=/home/user/anaconda3/envs/svision-pro-env/bin:$PATH 
{SVision-pro_path}/SVision-pro --target_path {tumorBam} --base_path {ControlBam} --genome_path {reference} --model_path {SVision-pro_path}/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {saveDir} --sample_name {sampleID} --preset error-prone --detect_mode somatic --min_supp 3  --min_mapq 5
python {SVision-pro_path}/extract_op.py --input_vcf {saveDir}/{sampleID}.svision_pro_v1.8.s3.vcf --extract somatic


### PacBio data:

#! /bin/bash
export PATH=/home/user/anaconda3/envs/svision-pro-env/bin:$PATH 
{SVision-pro_path}/SVision-pro --target_path {tumorBam} --base_path {ControlBam} --genome_path {reference} --model_path {SVision-pro_path}/src/pre_process/model_liteunet_1024_8_16_32_32_32.pth --out_path {saveDir} --sample_name {sampleID} --preset hifi --detect_mode somatic --min_supp 3  --min_mapq 5
python {SVision-pro_path}/extract_op.py --input_vcf {saveDir}/{sampleID}.svision_pro_v1.8.s3.vcf --extract somatic
