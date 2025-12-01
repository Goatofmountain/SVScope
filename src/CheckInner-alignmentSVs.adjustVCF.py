'''
Error Filter out 
Version 29 bug report:
- Continuous parameters could not filter out the none overlapped candidate windows 
- Make a vcf filter program to fillter out Report do not appeared in CandidateSpan.tumor.merged.somatic.bed 

'''
import os,re 
import numpy as np 
import pandas as pd 
import time
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def CheckFileLen(AimFile):
    with open(AimFile) as input:
        record = input.readlines()
    return(len(record))

def AdjustTandemRepeatINS(savedir):
    '''
    Tandem repeat result adjustment program
    input: savedir for TDScope Result 
    Output: 
        bedtools intersect file LC.somatic.bed vs CandidateSpan.tumor.merged.somatic.bed 
        dataframe recording the UUID relationship between LC region and CandidateSpanRegion
    '''
    # vcfFile = '{savedir}/{sampleID}_tumor.mergedSomatic.vcf'.format(savedir=savedir, sampleID=os.path.basename(savedir))
    FileName = [x for x in os.listdir(savedir) if re.search('mergedSomatic', x)][0]
    vcfFile = f'{savedir}/{FileName}'
    tmpFile = '{savedir}/somaticINS.region.bed'.format(savedir=savedir)
    output = open(tmpFile,'w')
    with open(vcfFile) as intput:
        WindowRecord = []
        for record in intput.readlines():
            if not re.search("#", record):
                UUID = record.strip().split("\t")[2]
                if UUID.split(".")[1] == 'INS':
                    region = "\t".join(re.split("[_-]", UUID.split(".")[-1])) + "\t10\t10\tINS\n"
                    output.write(region)
                    output.flush()
    output.close() 
    cmd = 'bedtools intersect -a {savedir}/CandidateSpan.tumorLC.merged.decision.somatic.bed -b {tmpFile} -wa -wb -F 1.0 > {savedir}/LC.vs.Candidate.bed'.format(
        savedir=savedir, tmpFile=tmpFile
    )
    os.system(cmd)
    Overlap = CheckFileLen('{savedir}/LC.vs.Candidate.bed'.format(savedir=savedir))
    if Overlap >0:
        df = pd.read_csv('{savedir}/LC.vs.Candidate.bed'.format(savedir=savedir), header=None, sep="\t")
        df.columns = ['chrom(LC)', 'start(LC)', 'end(LC)', 'Tread(LC)', 'Nread(LC)', 'Type(LC)', 
                      'chrom(Span)', 'start(Span)', 'end(Span)', 'Tread(Span)', 'Nread(Span)', 'Type(Span)']
        df['LCUUID'] = 'SVscope.'+df['Type(LC)']+"."+df['chrom(LC)']+"_"+df['start(LC)'].apply(str)+"-"+df['end(LC)'].apply(str)
        df['SpanUUID'] = 'SVscope.'+df['Type(Span)']+"."+df['chrom(Span)']+"_"+df['start(Span)'].apply(str)+"-"+df['end(Span)'].apply(str)
        return(df[['LCUUID', 'SpanUUID']])
    else:
        return(pd.DataFrame([]))

def VcfWindowLoading(savedir, excludechrom='chrM'):
    '''
    Load {savedir}/{savedir}_tumor.mergedSomatic.vcf file
    Record Head and 
    '''
    TRIDf = AdjustTandemRepeatINS(savedir)
    FileName = [x for x in os.listdir(savedir) if re.search('mergedSomatic', x)][0]
    vcfFile = f'{savedir}/{FileName}'
    # vcfFile = '{savedir}/{sampleID}_tumor.mergedSomatic.vcf'.format(savedir=savedir, sampleID=os.path.basename(savedir))
    if TRIDf.shape[0] > 0:
        TRIRegionUUID = np.unique(TRIDf['LCUUID'])
        TRIRegionUUID_small = np.unique(TRIDf['SpanUUID'])
    else:
        TRIRegionUUID = np.array([])
        TRIRegionUUID_small = np.array([])
    with open(vcfFile) as intput:
        Header = []
        WindowRecord = []
        WritedUUID = []
        for record in intput.readlines():
            if re.search('#', record):
                Header.append(record)
            else:
                info = record.strip().split("\t")
                UUID = info[2]
                if not re.search(excludechrom, UUID):
                    Type = UUID.split(".")[1]
                    if Type == 'INS':
                        if (UUID in TRIRegionUUID) and (UUID not in WritedUUID):
                            WindowRecord.append(record)
                            WritedUUID.append(UUID)
                        elif (UUID not in TRIRegionUUID_small) and (UUID not in WritedUUID):
                            WindowRecord.append(record)
                            WritedUUID.append(UUID)
                    elif Type == 'DEL':
                        if UUID not in WritedUUID:
                            WindowRecord.append(record)
                            WritedUUID.append(UUID)
                    else:
                        WindowRecord.append(record)
    vcf_adj = '{savedir}/{sampleID}_tumor.mergedSomatic.adjusted.vcf'.format(savedir=savedir, sampleID=os.path.basename(savedir))
    with open(vcf_adj, 'w') as output:
        for R in Header+WindowRecord:
            output.write(R)
            output.flush()
    return(vcf_adj)

def main(args):
    savedir = args.savedir
    VcfWindowLoading(savedir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    args = parser.parse_args()
    main(args)



