#!/bin/bash

TumorCell=$1
NormalCell=$2
MSA=$3
saveDir="/NAS/wg_tkl/SVScope_Data/SVScope_revise1017_CASTLEONT"
bamDir="/NAS/wg_tkl/SVScope_Data/bam"
ref="/NAS/wg_tkl/PanCancer_TKL/PanCancerRef/hg38_mainChr.fa"
mainScript="/NAS/wg_tkl/SVScope_Data/Code/SVScope/src/SVscope.py"

python ${mainScript} DataPrepare \
    -T ${bamDir}"/"${TumorCell}".bam" \
    -N ${bamDir}"/"${NormalCell}".bam" \
    -t ${TumorCell} \
    -n ${NormalCell} \
    -r ${ref} \
    -s ${saveDir}"/"${TumorCell}"_"${MSA} \
    -p 64 --selectwindows --FullProcess --cleanupDat \
    -M ${MSA} --platform ONT

# run with :
# /usr/bin/time -v ./Run_pipeline.sh HCC1937 HCC1937BL abPOA 2>&1|tee HCC1937_CASTLE_abPOA.log