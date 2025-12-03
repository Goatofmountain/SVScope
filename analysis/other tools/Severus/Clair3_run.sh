### ont data:

#! /bin/bash
export PATH=/home/user/anaconda3/envs/clair3/bin:$PATH 
{Clair3_path}/run_clair3.sh --bam_fn={controlBam} --ref_fn={reference} --threads={threads} --platform="hifi" --model_path={Clair3_path}/hifi_revio --output={saveDir}  --enable_phasing --longphase_for_phasing 



### PacBio data:

#! /bin/bash
export PATH=/home/user/anaconda3/envs/clair3/bin:$PATH 
{Clair3_path}/run_clair3.sh --bam_fn={controlBam} --ref_fn={reference} --threads={threads} --platform="hifi" --model_path={Clair3_path}/hifi_revio --output={saveDir}  --enable_phasing --longphase_for_phasing 
