# SVScope
Resolving somatic SVs via full-length sequence model-based local graph-genome optimization 
## Introduction
### A somatic structural variation caller based on long-read technology
The SVScope is a computational framework that leverages full-length sequence information and local graph genome optimization to accurately detect somatic SVs. The framework utilizes read alignment breakpoint information from the whole-genome scale to cluster and identify split-alignment somatic SVs and candidate inner-alignment somatic SVs. To mitigate the impact of alignment errors on inner-alignment somatic SV detection, SVScope re-analyses the alignment relationships among all full-length sequences spanning the candidate somatic SV interval using a partial order alignment (sPOA) graph with multi-sequence alignment representation and accurately clusters reads with a sequence mixture model. To avoid read coordination errors affected by centromeres, telomeres, and segmental duplication sequences, SVScope also implements a random forest machine learning approach based on local alignment features to filter high-confidence somatic SV events. 

---
## Features
- Utilizes full-length sequence from long-read data
- Optimizes local graph genome with sequence mixture model 
- Provides a detailed sequence and support read ID for each component of the local genome, including somatic components

## Installation
To get started with this software, follow these steps:

### Dependencies
Ensure you have the following Python packages installed:
- python=3.8
- pysam version 0.22.1
- pyspoa version 0.2.1
- pyabpoa >= 1.5.3
- numpy version 1.21.5
- scipy version 1.7.3
- sklearn version 1.0.2
- pandas
- matplotlib
- bedtools 
- biopython
- statsmodels
- mafft

Conda users can install all dependencies with conda command:
```shell
conda env create -f environment.yml
```
By default, this command line will build up a new conda environment named Release-env. Once activated, all dependencies will be ready for work:
```shell
conda activate SVScope-env
```

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
python src/SVScope.py DataPrepare \
    -T <CaseBam> \                                     # Path of Case sample long-read data alignment data in bam format. Using "," to divide if there are multiple files, like -T <bam1>,<bam2>. We recommend using minimap2.22+ for reads alignment.
    -N <ControlBam> \                                  # Path of Control sample long-read data alignment data in bam format. Using "," to divide if there are multiple files, like -N <bam1>,<bam2>. We recommend using minimap2.22+ for reads alignment.
    -t <CaseID> \                                      # Case SampleID. Using "," to divide if there are multiple files, like -t sample1,sample2. The length of CaseID should be the same as the CaseBam file.
    -n <ControlID> \                                   # Control SampleID. Using "," to divide if there are multiple files, like -t sample1,sample2. The length of the CaseID should be the same as the ControlBam file.
    -r <Reference sequence> \                          # Reference file in fasta format. An index file in fai format is also required in the same path.
    -s <SaveDir> \                                     # Path for result output
    -p <Thread>                                        # Number of CPU used for calculation
    --selectwindows \                                  # If set, select the candidate window first.
    --FullProcess \                                    # Run full SVScope pipeline from window selection -> split-alignment SV calling -> Local graph and read clustering -> somatic confidence checking process, --selectwindows parameter should be set first.
    --cleanupDat \                                     # If set, remove tmp files of SVScope 
    -M <MSA> \                                         # MSA algorithm used in local graph process, user can choose 'spoa' for SIMD POA, abPOA for abPOA and mafft for mafft, by default abPOA.
    --platform <LRS platform>                          # Platform choose, now we support "ont" and "PacBio", by default "PacBio".

```
### output
```shell
- <CaseID>.vs.<ControlID>.Raw.bed                      # Local genome component phasing result of TDScope in bed format, consisting of 10 columns.
- <CaseID>.vcf                                         # Result of raw inner-alignment somatic SV calling, including INS and DEL without randomforest selection.
- InterALNSVs.vcf                                      # Result of split-alignment somatic SV calling including BND, DUP and INV.
- <CaseID>.mergedSomatic.vcf                           # Result of all somatic SV calling with randomforest selection.
```

#### \<CaseID\>.vs.\<ControlID\>.Raw.bed Description
| Column | Name | Description |
|--------|------|-------------|
| 1      | Chromosome | Chromosome name |
| 2      | Start | Start position of the region |
| 3      | End | End position of the region |
| 4      | Somatic sequence | Sequence of somatic genome, split by ";" if more than 1 component |
| 5      | Support reads for somatic component | ID of reads supporting somatic genome components, split by ";" if more than 1 component |
| 6      | Number of somatic component | Number of somatic component |
| 7      | Germline sequence | Sequence of germline genome, split by ";" if more than 1 components |
| 8      | Support reads for germline component | ID of reads supporting germline genome components, split by ";" if more than 1 component |
| 9      | Number of germline component | Number of germline component |
| 10     | Label | Label of interval |

### Option1: Change tandem repeat and low complex window 
By default, SVScope will use pre-defined low complexity and tandem repeat window annotated by RepeatMasker at SVScope/doc/hg38.RepeatMasker.TD.Low.mainChr.sort.bed;
By default, SVScope will use pre-defined background 10kb window at SVScope/doc/hg38_mainChr.10kb.window.bed for coverage normalization. 
Both window lists mentioned above are made for the human genome hg38. User can change window list with parameter -D for low complexity and tandem repeat window and -W for background window like:
```bash
python src/SVScope.py \
-D <TANDEMREPEATFILE> \                                # low complexity and tandem repeat window annotated by RepeatMasker 
-W <GENOMEWINDOW> \                                    # background 10kb window identified by bedtools makewindows 
   DataPrepare \
    -T <CaseBam> \                                     # Path of Case sample long-read data alignment data in bam format. Using "," to divide if there are multiple files, like -T <bam1>,<bam2>. We recommend using minimap2.22+ for reads alignment.
    -N <ControlBam> \                                  # Path of Control sample long-read data alignment data in bam format. Using "," to divide if there are multiple files, like -N <bam1>,<bam2>. We recommend using minimap2.22 for reads alignment.
    -t <CaseID> \                                      # Case SampleID. Using "," to divide if there are multiple files, like -t sample1,sample2. The length of CaseID should be the same as the CaseBam file.
    -n <ControlID> \                                   # Control SampleID. Using "," to divide if there are multiple files, like -t sample1,sample2. The length of the CaseID should be the same as the ControlBam file.
    -r <Reference sequence> \                          # Reference file in fasta format. A index file in fai format is also required in the same path.
    -s <SaveDir> \                                     # Path for result output
    -p <Thread>                                        # Number of CPU used for calculation
    --selectwindows \                                  # If set, select candidate window first.
    --FullProcess \                                    # Run full SVScope pipeline from window selection -> split-alignment SV calling -> Local graph and read clustering -> somatic confidence checking process, --selectwindows parameter should be set first.
    --cleanupDat \                                     # If set, remove tmp files of SVScope 
    -M <MSA> \                                         # MSA algorithm use in local graph process, user can choose from 'spoa' for SIMD POA, abPOA for abPOA and mafft for mafft, by default abPOA.
    --platform <LRS platform>                          # Platform choose, now we support "ont" and "PacBio", by default "PacBio".

```
### Option2: Interval Correction for Low Complexity and Tandem Repeat Regions
- Rationale: 
During our analysis, we observed that for somatic Insertions (INS) occurring within Low Complexity (LC) or Tandem Repeat (TR) regions, a single breakpoint coordinate (or even a median position) often fails to fully represent the genomic alteration due to alignment ambiguities and local sequence repetition.
We propose that the entire LC or TR genomic interval serves as a more robust and biologically meaningful representation for these specific somatic events.

- Script Function: 
To facilitate this representation for downstream analysis, we provide a correction tool: CheckInner-alignmentSVs.adjustVCF.py. This script performs a post-hoc analysis on the SVScope VCF output. It checks if a detected somatic INS falls within a targeted LC or TR region. If such an overlap is confirmed, the script merges the specific INS breakpoint record with the corresponding full LC/TR genomic interval, converting the call into a region-based format.

- Usage: 
```shell
# Optional: Convert specific INS breakpoints to full LC/TR intervals
python src/CheckInner-alignmentSVs.adjustVCF.py \
    -s <SaveDir>   # Path to the directory containing SVScope results
```
> Note: This step is optional and designed for users who prefer interval-based annotations for repetitive regions. The performance benchmarks in the SVScope manuscript were calculated based on the precise breakpoints output by the main pipeline.
#### Output:
```shell
- <CaseID>.mergedSomatic.adjusted.vcf                  # result of all somatic SV calling with randomforest selection, INS located within tandem repeat and low complexity regions are represented as regions annotated by Repeat Masker.
```
