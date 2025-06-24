import os,re
import pandas as pd 
import numpy as np 
import pysam 
import argparse
from DecisionMaker import *
from DataScanner import *
import PairwiseCompare
import SomTDDetector
from OutVCF import *
import functools
from multiprocessing import Pool
import time
import logging
import gzip
import sqlite3
from joblib import dump, load
from sklearn.model_selection import StratifiedShuffleSplit, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import subprocess
import multiprocessing
import WindowSelection_v8_bam
import SomTDDetector_AimDatFetch
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Support Functions
def check_arguments(args):
    """
    check args logic
    """
    if args.FullProcess and not args.genomeWindow:
        parser.error("--FullProcess requires --genomeWindow to be specified.")
        return(1)
    return(0)

# Major functions 
def DataPrepare(args):
    # Prepare data for TDScope 
    check_res = check_arguments(args)
    if check_res != 0:
        sys.exit(check_res)
    logging.info('Data prepare : Start working')
    args.faiFile = args.Reference + ".fai"
    if not os.path.exists(args.savedir):
        os.system('mkdir -p %s' % args.savedir)
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    # Update V31: remove bed.gz as tmp file to optimize time cost
    # TbedFile = ",".join([os.path.join(args.savedir, '%s.bed.gz' % T) for T in TsampleID])
    # NbedFile = ",".join([os.path.join(args.savedir, '%s.bed.gz' % N) for N in NsampleID])
    tumorbamList = args.Tumorbam.split(",")
    normalbamList = args.Normalbam.split(",")
    # # Prepare bed.gz files
    # bamTobed_Process = []
    # if not os.path.exists(TbedFile.split(",")[-1]):
    #     for Tbam,Tbed in zip(tumorbamList, TbedFile.split(",")):
    #         cmd = 'bedtools bamtobed -i {Tbam} -cigar | bgzip > {Tbed} && tabix {Tbed}'.format(Tbam=Tbam, Tbed=Tbed)
    #         process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #         bamTobed_Process.append((process, Tbed))
    # if not os.path.exists(NbedFile.split(",")[-1]):
    #     for Nbam, Nbed in zip(normalbamList, NbedFile.split(",")):
    #         cmd = 'bedtools bamtobed -i {Nbam} -cigar | bgzip > {Nbed} && tabix {Nbed}'.format(Nbam=Nbam, Nbed=Nbed)
    #         process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #         bamTobed_Process.append((process, Nbed))
    # while bamTobed_Process:
    #     for process, FileName in bamTobed_Process[:]:
    #         if process.poll() is not None:
    #             logging.info('Data prepare : bamTobed process for {FileName} finished  return code {returncode}'.format(
    #                 FileName = FileName, returncode = process.returncode
    #             ))
    #             bamTobed_Process.remove((process, FileName))
    #     if bamTobed_Process:
    #         time.sleep(30)
    # # bed.gz File finished 
    # args.bedFileTumor = TbedFile
    # args.bedFileNormal = NbedFile
    dbFile_tumor_name = os.path.join(args.savedir, "Tumor")
    dbFile_normal_name = os.path.join(args.savedir, "Normal")
    # Main process continues to do makeupDB
    logging.info('Data prepare : bed to sqlite database process config background processes')
    background_processes = []
    background_processes.append(multiprocessing.Process(target=makeupDB_bam, args=(args.Tumorbam, dbFile_tumor_name)))
    background_processes.append(multiprocessing.Process(target=makeupDB_bam, args=(args.Normalbam, dbFile_normal_name)))
    for process in background_processes:
        process.start()
    logging.info('Data prepare : bed to sqlite database process start in background processes')
    # Make Candidate windows 
    if args.selectwindows:
        args.windowBed = WindowSelection_v8_bam.main(args)
        if args.saveData:
            SomTDDetector_AimDatFetch.main(args)
            for process in background_processes:
                process.join()
        elif args.FullProcess:
            args.rawBedFile = localGraph(args)
            for process in background_processes:
                process.join()
            mergeVCF = AlnFeature(args)
            if args.cleanupDat and os.path.exists(mergeVCF):
                TsampleID = args.TSampleID.split(",")
                NsampleID = args.NSampleID.split(",")
                # Update 31: remove bed.gz file as tmp file to optimize time cost 
                # TbedFile = [os.path.join(args.savedir, '%s.bed.gz' % T) for T in TsampleID]
                # NbedFile = [os.path.join(args.savedir, '%s.bed.gz' % N) for N in NsampleID]
                # for B in TbedFile + NbedFile:
                #     os.system('rm %s' % B)
                os.system('rm %s' % os.path.join(args.savedir, 'Tumor.sqlite'))
                os.system('rm %s' % os.path.join(args.savedir, 'Normal.sqlite'))
        else:
            for process in background_processes:
                process.join()
    else:
        for process in background_processes:
            process.join()
    logging.info('Data prepare : All process finished')

def localGraph_npz(args):
    # Handle the logic for the localGraph module
    logging.info('Local Graph : Start working')
    start_time = time.time()
    # somatic repeat calling software main program 
    ## First check input parameters 
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    npzFileList = [os.path.join(args.savedir, x) for x in os.listdir(args.savedir) if re.search('npz', x)]
    offset = int(args.offset)
    mapQ = int(args.mapQ)
    rawoutput = '%s.vs.%s.TandemRepeat.Raw.bed' % ("-".join(TsampleID), '-'.join(NsampleID))
    file_path = os.path.join(args.savedir, rawoutput)
    # Update V29: Add work continuous parameter to avoid breakpoint 
    FinishedTDRecord = []
    if args.Continue and os.path.exists(file_path):
        with open(file_path) as input:
            FinishedTDRecord = ["\t".join(x.strip().split("\t")[0:3]) for x in input.readlines()]
    ## Run program 
    if int(args.thread) > 6:
        P = Pool(processes=6)
    else:
        P = Pool(processes=int(args.thread))
    results = []
    for npzFile in npzFileList:
        Dat = np.load(npzFile, allow_pickle=True)['DatSet']
        for I in range(Dat.shape[0]):
            sequenceList, ReadIDs, flank_5, flank_3, TDRecord = Dat[I]
            # Update V29: Add work continuous parameter to avoid breakpoint 
            if FinishedTDRecord:
                if TDRecord not in FinishedTDRecord:
                    result = P.apply_async(SomTDDetector.TDscope_npz, (TDRecord, sequenceList, ReadIDs, flank_5,flank_3))
                    results.append(result)
            else:
                result = P.apply_async(SomTDDetector.TDscope_npz, (TDRecord, sequenceList, ReadIDs, flank_5,flank_3))
                results.append(result)
    # Update V29: Add work continuous parameter to avoid breakpoint 
    if FinishedTDRecord:
        f = open(file_path, 'a')
    else:
        f = open(file_path, 'w')
    outRecord = 0
    while results:
        for result in results:
            if result.ready():
                output = result.get()
                f.write("\t".join([str(x) for x in output]) + '\n')
                f.flush()
                results.remove(result)
                outRecord += 1
    f.close()
    os.system('sort -k1,1 -k2,2n {outDir}/{rawoutput} -o {outDir}/{rawoutput}'.format(outDir=args.savedir, rawoutput=rawoutput))
    time_span = (time.time() - start_time) / 3600
    logging.info(f'Local Graph : work finished with {time_span} hour')
    return(file_path)

def localGraph(args):
    # Handle the logic for the localGraph module
    logging.info('Local Graph : Start working')
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
    rawoutput = '%s.vs.%s.TandemRepeat.Raw.bed' % ("-".join(TsampleID), '-'.join(NsampleID))
    file_path = os.path.join(args.savedir, rawoutput)
    with open(args.windowBed) as bedin:
        TDRecordList = ["\t".join(x.strip().split("\t")) for x in bedin.readlines()]
    Decision_exe = functools.partial(Decision, Tlabel='tumor', readcutoff=3, hcutoff=3, scutoff=0.05)
    DataMaker_exe = functools.partial(DataMaker, refFile=args.Reference, bamFileList=bamFileList, LabelList=LabelList, offset=offset,mapQ=mapQ)
    DataMaker2_exe = functools.partial(DataMaker2, refFile=args.Reference, bamFileList=bamFileList, LabelList=LabelList, offset=offset,mapQ=mapQ)
    if not os.path.exists(args.savedir):
        os.mkdir(args.savedir)
    ## Run program 
    if int(args.thread) > 6:
        P = Pool(processes=6)
    else:
        P = Pool(processes=int(args.thread))
    results = []
    for TDRecord in TDRecordList:
        result = P.apply_async(SomTDDetector.TDscope, (TDRecord,DataMaker_exe,DataMaker2_exe,Decision_exe,))
        results.append(result)
    with open(file_path, 'w') as f:
        outRecord = 0
        while results:
            for result in results:
                if result.ready():
                    output = result.get()
                    f.write("\t".join([str(x) for x in output]) + '\n')
                    f.flush()
                    results.remove(result)
                    outRecord += 1
    os.system('sort -k1,1 -k2,2n {outDir}/{rawoutput} -o {outDir}/{rawoutput}'.format(outDir=args.savedir, rawoutput=rawoutput))
    time_span = (time.time() - start_time) / 3600
    logging.info(f'Local Graph : work finished with {time_span} hour')
    return(file_path)

def AlnFeature(args):
    # Handel the logic for the AlnFeature module 
    logging.info('Alignment feature collection module : Start working')
    if not os.path.exists(args.savedir):
        os.system('mkdir -p %s' % args.savedir)
    TsampleID = args.TSampleID.split(",")
    NsampleID = args.NSampleID.split(",")
    # Update V31: remove bed.gz tmp file process to optimize time cost 
    # TbedFile = ",".join([os.path.join(args.savedir, '%s.bed.gz' % T) for T in TsampleID])
    # NbedFile = ",".join([os.path.join(args.savedir, '%s.bed.gz' % N) for N in NsampleID])
    tumorbamList = args.Tumorbam.split(",")
    normalbamList = args.Normalbam.split(",")
    # Update V31: remove bed.gz tmp file process to optimize time cost 
    # # Update Version 19: Consider multiple file status 
    # if not os.path.exists(TbedFile.split(",")[-1]):
    #     for Tbam,Tbed in zip(tumorbamList, TbedFile.split(",")):
    #         os.system('bedtools bamtobed -i {TbamFile} -cigar | bgzip > {TbedFile} && tabix {TbedFile}'.format(TbamFile=Tbam, TbedFile=Tbed))
    # if not os.path.exists(NbedFile.split(",")[-1]):
    #     for Nbam,Nbed in zip(normalbamList, NbedFile.split(",")):
    #         os.system('bedtools bamtobed -i {NbamFile} -cigar | bgzip > {NbedFile} && tabix {NbedFile}'.format(NbamFile=Nbam, NbedFile=Nbed))
    # get sqlite3 database for each bed.gz file # Update V19: Consider multiple bed status, if more than one data, merge all bed.gz record into one sqlite
    dbFile_tumor = os.path.join(args.savedir, "Tumor") + ".sqlite"
    dbFile_normal = os.path.join(args.savedir, "Normal") + '.sqlite'
    # if not os.path.exists(dbFile_tumor):
    #     dbFile_tumor = makeupDB(TbedFile, os.path.join(args.savedir, "Tumor"))
    #     dbFile_normal = makeupDB(NbedFile, os.path.join(args.savedir, "Normal"))
    # Update V31: remove bed.gz tmp file process to optimize time cost 
    if not os.path.exists(dbFile_tumor):
        dbFile_tumor = makeupDB(args.Tumorbam, os.path.join(args.savedir, "Tumor"))
        dbFile_normal = makeupDB(args.Normalbam, os.path.join(args.savedir, "Normal"))
    # Update V31: remove bed.gz tmp file process to optimize time cost
    # Load COV and MapQ 
    # bg_df_T = background(args.genomeWindow, TbedFile, dbFile_tumor, showchromSpan=False, workthread=int(args.thread))
    # bg_df_N = background(args.genomeWindow, NbedFile, dbFile_normal, showchromSpan=False, workthread=int(args.thread))
    # SV_df_T = background(args.rawBedFile, TbedFile, dbFile_tumor, showchromSpan=True, workthread=int(args.thread))
    # SV_df_N = background(args.rawBedFile, NbedFile, dbFile_normal, showchromSpan=True, workthread=int(args.thread))
    # Load COV and MapQ 
    bg_df_T = background_bam(args.genomeWindow, args.Tumorbam, dbFile_tumor, showchromSpan=False, workthread=int(args.thread))
    bg_df_N = background_bam(args.genomeWindow, args.Normalbam, dbFile_normal, showchromSpan=False, workthread=int(args.thread))
    SV_df_T = background_bam(args.rawBedFile, args.Tumorbam, dbFile_tumor, showchromSpan=True, workthread=int(args.thread))
    SV_df_N = background_bam(args.rawBedFile, args.Normalbam, dbFile_normal, showchromSpan=True, workthread=int(args.thread))
    # Normalization 
    bg_df_T_nafree = bg_df_T.dropna()
    bg_df_N_nafree = bg_df_N.dropna()
    SV_df_T['COV_Zscore'] = SV_df_T['COV'].apply(lambda x: (x-np.mean(bg_df_T_nafree['COV'])) / np.std(bg_df_T_nafree['COV']))
    SV_df_T['mapQ_Zscore'] = SV_df_T['mapQRate'].apply(lambda x: (x-np.mean(bg_df_T_nafree['mapQRate'])) / np.std(bg_df_T_nafree['mapQRate']))
    SV_df_N['COV_Zscore'] = SV_df_N['COV'].apply(lambda x: (x-np.mean(bg_df_N_nafree['COV'])) / np.std(bg_df_N_nafree['COV']))
    SV_df_N['mapQ_Zscore'] = SV_df_N['mapQRate'].apply(lambda x: (x-np.mean(bg_df_N_nafree['mapQRate'])) / np.std(bg_df_N_nafree['mapQRate']))
    SV_df_T = SV_df_T.drop_duplicates()
    SV_df_N = SV_df_N.drop_duplicates()
    SV_df_T.index = SV_df_T['window'].apply(lambda x: "_".join(x.split("_")[:2])+"-" + x.split("_")[-1])
    SV_df_N.index = SV_df_N['window'].apply(lambda x: "_".join(x.split("_")[:2])+"-" + x.split("_")[-1])
    # Get Mis Score 
    SeqCompareDf_Filter = PairwiseCompare.MisScorePipe(args.rawBedFile).drop_duplicates()
    SeqCompareDf_Filter.index = SeqCompareDf_Filter['chrom'] + "_" + SeqCompareDf_Filter['start'].apply(str) + "-" + SeqCompareDf_Filter['end'].apply(str)
    SeqCompareDf_Filter['ABSMisScore'] = SeqCompareDf_Filter['MisScore'].apply(lambda x: np.abs(x))
    SeqCompareOut = os.path.join(args.savedir, '%s.Somatic.bed' % "_".join(TsampleID))
    SeqCompareDf_Filter.to_csv(SeqCompareOut, sep="\t", index=False, header=None)
    df_SVwindow = pd.read_csv(args.rawBedFile, header=None, sep="\t")
    df_SVwindow.columns = ['chrom', 'start','end', 'SomSeq','SomReads','SomCount','GermSeq','GermReads','GermCount','Label']
    df_SVwindow_Filter = df_SVwindow.loc[df_SVwindow['Label']=='NormalOutput|EMOutput'].drop_duplicates()
    df_SVwindow_Filter['window'] = df_SVwindow_Filter['chrom'] + "_" +df_SVwindow_Filter['start'].apply(str) +"-"+df_SVwindow_Filter['end'].apply(str)
    df_SVwindow_Filter.index = df_SVwindow_Filter['window']
    # Merge to get all Features for each candidate SV windows 
    SVwindowList = np.intersect1d(SeqCompareDf_Filter.index, df_SVwindow_Filter.index)
    # Calculate adapted reads percentage from Normal and tumor 
    ReadPool = pd.concat([SV_df_T.loc[SVwindowList, ['window', 'COV_Zscore', 'mapQ_Zscore', 'chromSpan']],
                          SV_df_N.loc[SVwindowList, ['COV_Zscore', 'mapQ_Zscore', 'chromSpan']],
                          df_SVwindow_Filter.loc[SVwindowList, 'SomReads'], 
                          df_SVwindow_Filter.loc[SVwindowList, 'SomReads'].apply(lambda x: [a.split("|")[-1] for a in ",".join(x.split(";")).split(",")]) + df_SVwindow_Filter.loc[SVwindowList,'GermReads'].apply(lambda x: [a.split("|")[-1] for a in ",".join(x.split(";")).split(",")]),
                          SV_df_T.loc[SVwindowList, 'TotalReadID'].apply(lambda x: x.split(",")),
                          SV_df_T.loc[SVwindowList, 'mapQRate'],
                          SV_df_N.loc[SVwindowList, 'TotalReadID'].apply(lambda x: x.split(",")),
                          SV_df_N.loc[SVwindowList, 'mapQRate'], 
                          SeqCompareDf_Filter.loc[SVwindowList, 'ABSMisScore']], axis=1)
    ReadPool.columns = ['window', 'COV_Tumor', 'mapQ_Tumor', 'chromSpan_Tumor', 'COV_Normal', 'mapQ_Normal', 'chromSpan_Normal', 'SomReads', 'AdaptReads', 'TotalRead_T', 'mapQRate_T', 'TotalRead_N', 'mapQRate_N', 'ABSMisScore']
    ReadPool['AdaptRatio_T'] = ReadPool.apply(lambda x: np.intersect1d(x['AdaptReads'], x['TotalRead_T']).shape[0] / (len(x['TotalRead_T'])*(1-x['mapQRate_T'])) if len(x['TotalRead_T'])*(1-x['mapQRate_T'])>0 else 0, axis=1)
    ReadPool['AdaptRatio_N'] = ReadPool.apply(lambda x: np.intersect1d(x['AdaptReads'], x['TotalRead_N']).shape[0] / (len(x['TotalRead_N'])*(1-x['mapQRate_N'])) if len(x['TotalRead_N'])*(1-x['mapQRate_N'])>0 else 0, axis=1)
    ReadPool['SupportReadSpanRatio'] = ReadPool.loc[ReadPool.index, 'SomReads'].apply(lambda x: spanchrRatio([a.split("|")[-1] for a in ",".join(x.split(";")).split(",")], dbFile_tumor))
    # Load random forest model to get predicition 
    script_dir = os.path.dirname(os.path.abspath(__file__))
    model = load(os.path.join(script_dir,'RandomForest.1218.WholeData8-2.FinalModel.joblib'))
    X_scaled = ReadPool[['COV_Tumor', 'mapQ_Tumor', 'COV_Normal', 'mapQ_Normal', 'ABSMisScore', 'chromSpan_Tumor', 'chromSpan_Normal', 'AdaptRatio_T', 'AdaptRatio_N', 'SupportReadSpanRatio']]
    yprob = model.predict_proba(X_scaled)[:, 1]
    y_hat = model.predict(X_scaled)
    ReadPool['yprob'] = yprob
    ReadPool['y_hat'] = y_hat
    ReadPoolOut = os.path.join(args.savedir, 'RandomForestResult.tsv')
    ReadPool.to_csv(ReadPoolOut, sep="\t")
    out_vcf = os.path.join(args.savedir, '%s.vcf' % "_".join(args.TSampleID.split(",")))
    bed2vcf(args.rawBedFile, SeqCompareOut, ReadPoolOut,out_vcf,args.TSampleID, args.Reference)
    # MergeVCFs 
    if len([x for x in os.listdir(args.savedir) if x.split(".")[-1]=='vcf']) >= 2:
        VCFList = [os.path.join(args.savedir, x) for x in os.listdir(args.savedir) if x.split(".")[-1]=='vcf']
        mergeVCF = os.path.join(args.savedir, '%s.mergedSomatic.vcf' % "_".join(args.TSampleID.split(",")))
        mergeOutput = open(mergeVCF, 'w')
        # rebuild head 
        with open(out_vcf) as rawVCF:
            headRecord = [x for x in rawVCF.readlines() if re.search('#', x)]
            for record in headRecord:
                if re.search('##FORMAT',record):
                    mergeOutput.write('##ALT=<ID=INV,Description="Invasion">\n##ALT=<ID=BND,Description="Translocation">\n'+record)
                else:
                    mergeOutput.write(record)
        mergeOutput.close()
        tmpVCF = os.path.join(args.savedir, 'tmp.vcf')
        os.system('grep "True" %s >> %s' % (out_vcf, tmpVCF))
        os.system('grep -v "#" %s >> %s' % (os.path.join(args.savedir, "InterALNSVs.vcf"), tmpVCF))
        os.system('sort -k1,1 -k2,2n %s >> %s' % (tmpVCF, mergeVCF))
        os.system('rm %s' % tmpVCF)
    return(mergeVCF)

def callsomaticSV(args):
    # Handel the logic for the callsomaticSV module 
    args.rawBedFile = localGraph(args)
    mergeVCF = AlnFeature(args)
    if args.cleanupDat and os.path.exists(mergeVCF):
        TsampleID = args.TSampleID.split(",")
        NsampleID = args.NSampleID.split(",")
        TbedFile = [os.path.join(args.savedir, '%s.bed.gz' % T) for T in TsampleID]
        NbedFile = [os.path.join(args.savedir, '%s.bed.gz' % N) for N in NsampleID]
        # for B in TbedFile + NbedFile:
        #     os.system('rm %s' % B)
        os.system('rm %s' % os.path.join(args.savedir, 'Tumor.sqlite'))
        os.system('rm %s' % os.path.join(args.savedir, 'Normal.sqlite'))
    # out_vcf = os.path.join(args.savedir, '%s.vcf' % args.TSampleID)
    # bed2vcf(args.rawBedFile, SeqCompareOut, ReadPoolOut,out_vcf,args.TSampleID, args.Reference)

def main():
    # default file path setting 
    doc_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../doc')
    # Main parser 
    parser = argparse.ArgumentParser(
        prog='TDScope.py',
        description='TDScope: A computational system for somatic SV calling with local graph genome optimization and whole genome alignment feature adjustment',
        epilog='''
        Use python TDScope.py <command> -h to view help information for a specific command
        Available commands:
            -DataPrepare: Run dataprepare process user can also choose to run directly until TDScope finished
            -localGraph: Run Local graph optimization module on candidate somatic SV window 
            -AlnFeature: collect feature from spanned reads of each potential somatic SV selected by localGrpah 
            -callsomaticSV: Call somatic SV and output in vcf format
        '''
    )
    parser.add_argument('-W', '--genomeWindow', type=str, default=os.path.join(doc_dir, 'hg38_mainChr.10kb.window.bed'), 
                        help="gernomic window file, by default 10kb window bed could fetch by using bedtools makewindows command")
    parser.add_argument("-D", "--tandemRepeatFile", type=str, default=os.path.join(doc_dir, 'hg38.RepeatMasker.TD.Low.mainChr.sort.bed'), 
                        help="BedFile annotation for tandem repeat region masked by repeatmasker")
    # Create sub-command parsers
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands:',
        description='These are the commands supported by TDScope.py:'
    )
    # Update Version 26: Add data prepare sub-command 
    parser_DataPrepare = subparsers.add_parser(
        'DataPrepare',
        help='Run window selction for candidate somatic SV window selection',
        description='This command should be run before callsomaticSV module'
    )
    parser_DataPrepare.add_argument("-T", "--Tumorbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_DataPrepare.add_argument("-N", "--Normalbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_DataPrepare.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser_DataPrepare.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser_DataPrepare.add_argument("-r", "--Reference", required=True, help="reference file fasta path")
    parser_DataPrepare.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser_DataPrepare.add_argument("-p", "--thread", required=True, help="CPU use for program")
    # parser_DataPrepare.add_argument("-D", "--tandemRepeatFile", type=str, default=os.path.join(doc_dir, 'hg38.RepeatMasker.TD.Low.mainChr.sort.bed'), help="BedFile annotation for tandem repeat region masked by repeatmasker")
    # parser_DataPrepare.add_argument('-W', '--genomeWindow', type=str, default=os.path.join(doc_dir, 'hg38_mainChr.10kb.window.bed'), help="gernomic window file, by default 10kb window bed could fetch by using bedtools makewindows command. Required if --FullProcess is specified")
    parser_DataPrepare.add_argument("-o", "--offset", type=int, default=50, help="offset default value is 50")
    parser_DataPrepare.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    parser_DataPrepare.add_argument("--selectwindows", action='store_true', default=False, help='If specified, TDScope will run windowselection process. Once saveData set True, this parameter should be specified first. Default: False')
    parser_DataPrepare.add_argument("--saveData", action='store_true', default=False, help='If specified, save intermediate read sequence data within candidate window into .npzFile. Default: False')
    parser_DataPrepare.add_argument('--FullProcess', action='store_true', default=False, help='If specified, run TDScope process after window selection')
    parser_DataPrepare.add_argument("-c", '--cleanupDat', action='store_true', default=False, help="If set, clean up bed.gz and sqlite files for space, by default False")
    parser_DataPrepare.set_defaults(func=DataPrepare)
    # Add LocalGraph sub-command 
    parser_localGraph = subparsers.add_parser(
        'localGraph',
        help='Run Local graph optimization module on candidate somatic SV window',
        description='This command is used to perform the localGraph module'
    ) 
    parser_localGraph.add_argument("-w", "--windowBed", required=True, help="pre made tandem repeat windows chrom,repeatStart,repeatEnd")
    parser_localGraph.add_argument("-T", "--Tumorbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_localGraph.add_argument("-N", "--Normalbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_localGraph.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser_localGraph.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser_localGraph.add_argument("-r", "--Reference", required=True, help="reference file fasta path")
    parser_localGraph.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser_localGraph.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser_localGraph.add_argument("-o", "--offset", type=int, default=50, help="offset default value is 50")
    parser_localGraph.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    parser_localGraph.set_defaults(func=localGraph)
    # Add LocalGraph_npz sub-command 
    parser_localGraph_npz = subparsers.add_parser(
        'localGraph_npz',
        help='Run Local graph optimization module on candidate somatic SV window',
        description='This command is used to perform the localGraph module'
    )
    parser_localGraph_npz.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser_localGraph_npz.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser_localGraph_npz.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser_localGraph_npz.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser_localGraph_npz.add_argument("-o", "--offset", type=int, default=50, help="offset default value is 50")
    parser_localGraph_npz.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    parser_localGraph_npz.add_argument("-C", '--Continue', action='store_true', default=False, help="If set, continue local graph work according to current Raw.bed File")
    parser_localGraph_npz.set_defaults(func=localGraph_npz)
    # Add AlnFeature sub-command 
    parser_AlnFeature = subparsers.add_parser(
        'AlnFeature',
        help='Run alignment feature collection and random forest process for candidate somatic SV window',
        description='This command is used to perform the AlnFeature module'
    ) 
    parser_AlnFeature.add_argument('-B', '--rawBedFile', required=True, help="path to RawBedFile from LocalGraph Step")
    # parser_AlnFeature.add_argument('-W', '--genomeWindow', required=True, help="gernomic window file, by default 10kb window bed could fetch by using bedtools makewindows command")
    parser_AlnFeature.add_argument("-T", "--Tumorbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_AlnFeature.add_argument("-N", "--Normalbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_AlnFeature.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser_AlnFeature.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser_AlnFeature.add_argument("-r", "--Reference", required=True, help="reference file fasta path must have fai file in the same pathway")
    parser_AlnFeature.add_argument("-s", "--savedir", required=True, help="dir for result file save the bed dir ")
    parser_AlnFeature.add_argument("-p", "--thread", required=True, help="CPU use for program")
    parser_AlnFeature.set_defaults(func=AlnFeature)
    # Add callsomaticSV process, a merge of localGraph and AlnFeature
    parser_callsomaticSV = subparsers.add_parser(
        'callsomaticSV',
        help='Run Full process of TDscope finally output the vcf format file'
    )
    parser_callsomaticSV.add_argument("-w", "--windowBed", required=True, help="pre made candidate somatic SV windows chrom,repeatStart,repeatEnd")
    parser_callsomaticSV.add_argument("-T", "--Tumorbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_callsomaticSV.add_argument("-N", "--Normalbam", required=True, help="ONT read alignment bam file, multiple bam should seperated with ','")
    parser_callsomaticSV.add_argument("-t", "--TSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with tumor bam")
    parser_callsomaticSV.add_argument("-n", "--NSampleID", required=True, help="SampleID for this somatic project, if multiple sample should seperated with ',' and have same length with normal bam")
    parser_callsomaticSV.add_argument("-r", "--Reference", required=True, help="reference file fasta path")
    parser_callsomaticSV.add_argument("-s", "--savedir", required=True, help="dir for result file save")
    parser_callsomaticSV.add_argument("-p", "--thread", required=True, help="CPU use for program")
    # parser_DataPrepare.add_argument('-W', '--genomeWindow', type=str, default=os.path.join(doc_dir, 'hg38_mainChr.10kb.window.bed'), help="gernomic window file, by default 10kb window bed could fetch by using bedtools makewindows command. Required if --FullProcess is specified")
    parser_callsomaticSV.add_argument("-o", "--offset", type=int, default=50, help="offset default value is 50")
    parser_callsomaticSV.add_argument("-q", "--mapQ", type=int, default=5, help="mapQ default value is 5")
    parser_callsomaticSV.add_argument("-c", '--cleanupDat', action='store_true', default=False, help="If set, clean up bed.gz and sqlite files for space, by default False")
    parser_callsomaticSV.set_defaults(func=callsomaticSV)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
