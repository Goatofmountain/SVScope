''''
__Author__: Kailing Tu
__Version__: v1.1.0
__ReleaseTime__: 2024-04-01

Release Note:
    Change default parameter mapQ = 5, offset=200
    Remove TD Number and time cutoff for program stop
    Merge pipeline DataMaker and DecisionMaker together for faster work 
    Using Cython script to accelerate datamaker and EM speed 

Requirement:
    os,re, argparse
    DecisionMaker : From our script 
    DataScanner : From Our script 
    UtilFunctions : From Our script 
'''
import os,re
import argparse
from DecisionMaker import *
from DataScanner import *
import functools
from multiprocessing import Pool
import time
import logging


def FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=200):
    # Input refFile, bamFile, bamLabel, TDRecord 
    # get TD associated sub sequence from each reads 
    TDchrom, TDStart, TDEnd = TDRecord.strip().split("\t")[0],int(TDRecord.strip().split("\t")[1]), int(TDRecord.strip().split("\t")[2])
    F5start,F5end,F3start,F3end = TDStart-offset, TDStart, TDEnd,TDEnd+offset
    readIDList, readTDSeq, FlankMQ = [],[],[]   # readID, readSequence, primary ALN mapQ
    for bamIDX in range(len(bamFileList)):
        tmpReadSeqArr = []                          # readID => readsequence extracted from primary alignment, mapQ
        F5_readsIDX, F3_readsIDX = [], []
        F5_readName, F3_readName = [], []
        for reads in pysam.AlignmentFile(bamFileList[bamIDX]).fetch(TDchrom, TDStart, TDEnd):
            # primary alignment reads get the whole reads sequence 
            if not(reads.is_secondary or reads.is_supplementary):
                tmpReadSeqArr.append([reads.query_name, reads.query_sequence, reads.mapq])
            # Work For F5 flank 
            if (reads.reference_start<F5start) and (reads.reference_end>F5end) and (not reads.is_secondary):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F5_readsIDX.append([reads.qname] +  ReadsLoci(reads, F5start, F5end, offset))
                F5_readName.append(reads.qname)
            # Work For F3 flank 
            if (reads.reference_start<F3start) and (reads.reference_end>F3end) and (not reads.is_secondary):
                offset = 0
                if reads.is_supplementary:
                    CIGAR = np.array(reads.cigartuples)
                    if CIGAR[0][0] == 5:
                        offset = CIGAR[0][1]
                F3_readsIDX.append([reads.qname] +  ReadsLoci(reads, F3start, F3end, offset))
                F3_readName.append(reads.qname)
        # remove reads 
        F5_N, F5_C = np.unique(F5_readName, return_counts=True)
        blackList_F5 = F5_N[np.where(F5_C>=2)[0]]
        F3_N, F3_C = np.unique(F3_readName, return_counts=True)
        blackList_F3 = F3_N[np.where(F3_C>=2)[0]]
        blackList = np.union1d(blackList_F3, blackList_F5)
        # Fetch 5Flank + TD + 3Flank sequence 
        if len(F5_readsIDX) * len(F3_readsIDX) * len(tmpReadSeqArr) > 0:
            spanReadIDs = np.intersect1d(np.intersect1d(np.array(tmpReadSeqArr)[:,0], np.array(F5_readsIDX)[:,0]), np.array(F3_readsIDX)[:,0])
            if blackList.shape[0] > 0:
                spanReadIDs = np.setdiff1d(spanReadIDs, blackList)
            if spanReadIDs.shape[0] >= 3:
                F5_df = pd.DataFrame(F5_readsIDX, columns=['readID', 'start', 'end'])
                F3_df = pd.DataFrame(F3_readsIDX, columns=['readID', 'start', 'end'])
                SeqDf = pd.DataFrame(tmpReadSeqArr, columns=['readID', 'qseq', 'mapQ'])
                SeqDf.index = SeqDf['readID']
                SummaryDf = pd.concat([F5_df.loc[F5_df['readID'].isin(spanReadIDs)].groupby(['readID'])['start'].apply(min), 
                                       F3_df.loc[F3_df['readID'].isin(spanReadIDs)].groupby(['readID'])['end'].apply(max), 
                                       SeqDf.loc[spanReadIDs, ['qseq', 'mapQ']]], axis=1)
                SummaryDf['SubSeq'] = SummaryDf.apply(lambda x: x['qseq'][x['start']:x['end']].replace("N",""), axis=1)
                readIDList += [LabelList[bamIDX] + "|" + x for x in SummaryDf.index]
                readTDSeq += [x for x in SummaryDf['SubSeq']]
                FlankMQ += [int(x) for x in SummaryDf['mapQ']]
    return(readTDSeq, readIDList, FlankMQ)

def DataMaker(TDRecord, refFile, bamFileList, LabelList, offset=200,mapQ=5):
    # input reference, bamfilelist, label, TDrecord
    # selecting reads and feature 
    # output seqdatamx for further clustering 
    # raw sequence generation 
    readTDSeq, readIDList, FlankMQ = FetchTDsubSeq(refFile, bamFileList, LabelList, TDRecord, offset=offset)
    CertainIDX = [I for I in range(len(FlankMQ)) if np.min(FlankMQ[I])>=mapQ]
    refFasta = pysam.FastaFile(refFile)
    chrom,start,end = TDRecord.strip().split("\t")[0:3]
    flank_5,flank_3 = refFasta.fetch(chrom, int(start)-offset, int(start)).upper(), refFasta.fetch(chrom, int(end),int(end)+offset).upper()
    exampleSequence = refFasta.fetch(chrom, int(start)-offset, int(end)+offset).upper()
    if re.search('N', flank_5) or re.search('N', flank_3) or re.search('N', exampleSequence) or (len(CertainIDX) <= 3):
        sequenceList = np.array([])
        ReadIDs = np.array([])
    else:
        # datafilter remove imcomplete flank sequence reads 
        readTDSeqNew = [readTDSeq[I] for I in CertainIDX]                                                # sub read selection uncertain flank aln reads would be removed
        ReadIDs = np.array([readIDList[I] for I in CertainIDX])
        sequenceList = [refFasta.fetch(chrom, int(start)-offset, int(end)+offset).upper()] + readTDSeqNew
    return(sequenceList, ReadIDs, flank_5,flank_3, TDRecord)

def TDscope(TDRecord, DataMaker, readcutoff=3, Tlabel="tumor"):
    '''
    Complete pipeline of TDscope 
    :param:
        TDRecord: bed record for Tandem Repeat window with \t split 
        DataMarker: data extractor from bam Files, preset by functools in the main pipeline 
        DecisionMaker: pipeline decide whether a TD region is somatic or not
    :return:
        Record: TD region record including chrom,start,end, somTDconsensus seq, somTD support reads, germTDconsensus seq, germTD support reads, pvalue
    '''
    start_time = time.time()
    logging.info(f"pipeline for region {TDRecord} start to work")
    sequenceList, ReadIDs, flank_5,flank_3, TDRecord = DataMaker(TDRecord)
    SaveRecord = np.array([sequenceList, ReadIDs, flank_5, flank_3, TDRecord], dtype=object)
    time_span = time.time() - start_time
    logging.info(f"pipeline for region {TDRecord} finished Take {time_span}s")
    return(SaveRecord)

def main(args):
    logging.info('Start working')
    start_time = time.time()
    # somatic repeat calling software main program 
    ## First check input parameters 
    tumorbamList = args.Tumorbam.split(",")
    normalbamList = args.Normalbam.split(",")
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    offset = int(args.offset)
    mapQ = int(args.mapQ)
    if len(TsampleID) != len(tumorbamList):
        print("SampleID not meet tumor bam file, exit !")
        exit(1)
    if len(NsampleID) != len(normalbamList):
        print("SampleID not meet normal bam file, exit !")
        exit(1)
    bamFileList = tumorbamList + normalbamList
    LabelList = [x+"_tumor" for x in TsampleID]  + [x+'_normal' for x in NsampleID]
    with open(args.windowBed) as bedin:
        TDRecordList = ["\t".join(x.strip().split("\t")[0:3]) for x in bedin.readlines()]
    DataMaker_exe = functools.partial(DataMaker, refFile=args.Reference, bamFileList=bamFileList, LabelList=LabelList, offset=offset,mapQ=mapQ)
    if not os.path.exists(args.savedir):
        os.mkdir(args.savedir)
    ## Run program
    if int(args.thread) > 4:
        P = Pool(processes=4)
    else:
        P = Pool(processes=int(args.thread))
    # P = Pool(processes=int(args.thread))
    results = []
    for TDRecord in TDRecordList:
        result = P.apply_async(TDscope, (TDRecord,DataMaker_exe,))
        results.append(result)
    batch = 0 
    dblockCount, dblockCountCutoff = 0, 8192
    batch_List = []
    max_runtime = 10 * 60 * 60
    start_time = time.time()
    while results:
        for result in results:
            if result.ready():
                output = result.get()
                if output.shape[0] > 0:
                    dblockCount += 1
                    batch_List.append(output)
                    if len(batch_List) >= dblockCountCutoff:
                        rawoutput = '%s.vs.%s.TandemRepeat.batch%s.npz' % ("-".join(TsampleID), '-'.join(NsampleID), batch)
                        np.savez("{outDir}/{rawoutput}".format(outDir=args.savedir, rawoutput=rawoutput), DatSet=np.array(batch_List))
                        batch_List = []
                        batch += 1
                        dblockCount = 0
                results.remove(result)
        if time.time() - start_time >= max_runtime:
            logging.info(f"Time Out, program will break and leave unfinished Data unfinished")
            break
    if len(batch_List) > 0:
        rawoutput = '%s.vs.%s.TandemRepeat.batch%s.npz' % ("-".join(TsampleID), '-'.join(NsampleID), batch)
        np.savez("{outDir}/{rawoutput}".format(outDir=args.savedir, rawoutput=rawoutput), DatSet=np.array(batch_List))
    time_span = (time.time() - start_time) / 3600
    logging.info(f'work finished with {time_span} hour')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-w", "--windowBed", required=True, help="pre made tandem repeat windows chrom,repeatStart,repeatEnd")
    parser.add_argument("-T", "--Tumorbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser.add_argument("-N", "--Normalbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser.add_argument("-r", "--Reference", required=True, help="reference file fasta path")
    parser.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser.add_argument("-o", "--offset", type=int, default=200, help="offset default value is 200")
    parser.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    args = parser.parse_args()
    main(args)


