# SVscope
Long-read-based somatic SV detector via full-length sequence model-based local graph-genome optimization
## Introduction
### A somatic structural variation caller based on long-read technology
The SVScope is a computational framework that leverages full-length sequence information and local graph genome optimization to accurately detect somatic SVs. The framework utilizes read alignment breakpoint information from the whole-genome scale to cluster and identify split-alignment somatic SVs and candidate inner-alignment somatic SVs. To mitigate the impact of alignment errors on inner-alignment somatic SV detection, SVScope re-analyzes the alignment relationships among all full-length sequences spanning the candidate somatic SV interval using a partial order alignment (sPOA) graph with multi-sequence alignment representation and accurately clusters reads with a sequence mixture model. To avoid read coordination errors affected by centromeres, telomeres, and segmental duplication sequences, SVScope also implements a random forest machine learning approach based on local alignment features to filter high-confidence somatic SV events. 

---
## Features
- Utilizes full-length sequence from long-read data
- Optimizes local graph genome with sequence mixture model 
- Provides detailed sequence and support read ID for each componenet of local genome including somatic conponents

## Installation
To get started with this software, follow these steps:

### Dependencies
Ensure you have the following Python packages installed:
- pysam version 0.19.1
- pyspoa version 0.2.1
- numpy version 1.21.5
- scipy version 1.7.3
- sklearn version 1.0.2
- Levenshtein version 0.23.0

### Installation Steps
- Clone the repository:
   ```bash
   git clone https://github.com/Goatofmountain/SVScope.git
   cd SVScope
   python src/SVScope.py -h 
   ```

## Usage
The SVScope algorithm consists of three main modules: the DataPrepare module initialize the detection process, generate candidate somatic SV intervals for further analysis, the local graph genome optimization module (`localGraph`), which optimizes the local graph genome, and the local graph confidence assessment module (`AlnFeature`), which evaluates the confidence of the local graph genome. We have designed the `callsomaticSV` module to link the above modules together to directly obtain the somatic SV calculation results in VCF format. The specific usage is as follows:
### command line
```bash
python src/SVscope.py} DataPrepare \ 
-D ${RepeatWindow} \                               # low complexity and tandem repeat window annotated by RepeatMasker, by default: SVScope/doc/hg38.RepeatMasker.TD.Low.mainChr.sort.bed
-W <Genome interval window> \                      # genomic window for normalization by default 10kb window, by default: SVScope/doc/hg38_mainChr.10kb.window.bed
-T <CaseBam> \                                     # Path of Case sample long-read data alignment data in bam format, we recommand to use minimap2.22 for reads alignment.
-N <ControlBam> \                                  # Path of Control sample long-read data alignment data in bam format, we recommand to use minimap2.22 for reads alignment.
-t <CaseID> \                                      # Case SampleID 
-n <ControlID> \                                   # Control SampleID
-r <Reference sequence> \                          # reference file in fasta format 
-s <SaveDir> \                                     # path for result output
--selectwindows \                                  # If specified, TDScope will run windowselection process. Once saveData set True, this parameter should be specified first. Default: False
--FullProcess \                                    # If specified, run TDScope process after window selection
--cleanupDat \                                     # If set, clean up bed.gz and sqlite files for space, by default False
-p <Thread>                                        # Number of CPU used for calculation

python src/CheckInner-alignmentSVs.adjustVCF.py \
-s <SaveDir> \                                     # path for result output
```
### output
- \<CaseID\>.vs.\<ControlID\>.Raw.bed             # Local genome component phasing result of TDScope in bed format, consist of 10 columns
- \<CaseID\>_tumor.vcf                            # result of raw inner-alignment somatic SV calling including INS and DEL without randomforest selection.
- \InterALNSVs.vcf                                # result of split-alignment somatic SV calling including BND, DUP and INV.
- \<CaseID>_tumor.mergedSomatic.vcf               # result of all somatic SV calling with randomforest selection.
- \<CaseID>_tumor.mergedSomatic.adjusted.vcf      # result of all somatic SV calling with randomforest selection, INS located within tandem repeat and low complexity regions are represented as regions annotated by Repeat Masker.


#### \<CaseID\>.vs.\<ControlID\>.Raw.bed Description
| Column | Name | Description |
|--------|------|-------------|
| 1      | Chromosome | Chromosome name |
| 2      | Start | Start position of the region |
| 3      | End | End position of the region |
| 4      | somatic Sequence | sequence of somatic genome, split by ";" if more than 1 components |
| 5      | support reads for somatic component | ID of reads supporting somatic genome components, split by ";" if more than 1 components |
| 6      | Number of somatic component | Number of somatic component |
| 7      | germline Sequence | sequence of germline genome, split by ";" if more than 1 components |
| 8      | support reads for germline component | ID of reads supporting germline genome components, split by ";" if more than 1 components |
| 9      | Number of germline component | Number of germline component |
| 10     | label | Label of interval |

