import os,re
import argparse
# import pysam 
import numpy as np
import pandas as pd
import time
from statsmodels.stats import multitest
from Bio.Seq import Seq
from Bio import Align
from scipy import stats
# from matplotlib import pyplot as plt
from collections import Counter

############################################
#       Define  Function
############################################
##### Define Funtion
def AligmentScore(SomConsensus, GerConsensus, cutoff=0):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1         
    aligner.mismatch_score = 0      
    aligner.open_gap_score = -1     
    aligner.extend_gap_score = -1   
    
    alignments = aligner.align(SomConsensus, GerConsensus)
    best_alignment = alignments[0]
    
    alig = str(best_alignment).split('\n')[1]
    
    if cutoff > 0:
        TD_alig = alig[cutoff:len(alig)-cutoff]
    else:
        TD_alig = alig
        
    Match = Counter(TD_alig)["|"]
    AligLen = len(TD_alig)
    MisScore = AligLen - Match
    return(MisScore)

def smaller_absolute_value(a, b):
    if abs(a) < abs(b):
        return(a)
    else:
        return(b)

def Mismatch_abs(callLine):
    somSeqList = callLine['somSeqList'].split(';')
    germSeqList = callLine['germSeqList'].split(';')
    Res = []
    for Som in somSeqList:
        Abs = 1000000000000000000000
        for Ger in germSeqList:
            score = len(Som) -len(Ger)
            AbsScore =smaller_absolute_value(Abs, score)
        Res.append(AbsScore)
    if len(Res)>1:
        Score= ';'.join(list(map(str, Res)))
    else:
        Score = str(Res[0])
    return(Score)

def CalculateMisscore(callLine):
    somSeqList = callLine['somSeqList'].split(';')
    germSeqList = callLine['germSeqList'].split(';')
    MisScore = 1000000000000000000000
    for Som in somSeqList:
        for Ger in germSeqList:
            score = AligmentScore(Som, Ger)
            if len(Som) < len(Ger):
                score = (-1) * score
            MisScore = smaller_absolute_value(MisScore,score)
    return(MisScore)

def CallAlleleFreq(SomaticTD):
    # input somatic record, calculate allele frequency
    Raw_arr = SomaticTD[['somSupportReadID', 'germSupportReadID']].to_numpy()
    somReadCountList = np.array([len(x.split(",")) for x in Raw_arr[0].split(";")])
    germReadList = np.concatenate([x.split(",") for x in Raw_arr[1].split(";")])
    germTumorReads = [x for x in germReadList if re.search('_tumor|', x)]
    N = np.sum(somReadCountList) + len(germTumorReads)
    AlleleFreq = ";".join([str(x) for x in somReadCountList / N])
    return(AlleleFreq)

def MisScorePipe(filepath):
    df = pd.read_csv(filepath, sep="\t", header=None)
    df.columns = ['chrom', 'start', 'end', 'somSeqList', 'somSupportReadID', 'someventCount','germSeqList', 'germSupportReadID', 'germeventCount','flag']
    somDf = df.loc[df['flag']=='NormalOutput|EMOutput']
    SomaticRes = pd.DataFrame(columns=['chrom','start','end','window','somSupportReadID','germSupportReadID','MisScore', 'AF'])
    if somDf.shape[0] > 0:
        somDf['window'] = somDf['chrom']+'_'+somDf['start'].astype('str')+'-'+somDf['end'].astype('str')
        somDf['MisScore'] = somDf.apply(lambda x: CalculateMisscore(x), axis=1)
        somDf['AF'] = somDf.apply(lambda x: CallAlleleFreq(x), axis=1)
        SomaticRes = somDf[['chrom','start','end','window','somSupportReadID','germSupportReadID','MisScore','AF']]
    return(SomaticRes)

def main(args):
    filepath = os.path.join(args.workDir, args.sampleID, "%s.vs.%s.TandemRepeat.Raw.bed" % (args.sampleID, args.sampleID))
    output = os.path.join(args.outputDir, '%s.Somatic.bed' % args.sampleID)
    if not os.path.exists(args.outputDir):
        os.system('mkdir %s' % args.outputDir)
    SomaticRes = MisScorePipe(filepath)
    SomaticRes.to_csv(output, sep="\t", header=None, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-w", "--workDir", required=True, help="work Dir")
    parser.add_argument("-s", "--sampleID", required=True, help="sampleID")
    parser.add_argument("-o", "--outputDir", required=True, help="output dir")
    args = parser.parse_args()
    main(args)



