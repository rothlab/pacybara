## Pre-requisites and installation:
1. Pacybara and barseqPro require a SLURM HPC cluster.
2. Install [`clusterutil`](https://github.com/jweile/clusterutil) 
3. Make sure R version 4 or higher and the following R packages are installed: `yogitools`, `yogiseq`, `hgvsParseR`,`argparser`,`pbmcapply`,`hash`,`bitops`
	* yogitools, yogiutil and hgvsParser can be installed via `remotes::install_github("jweile/yogitools")`, etc
4. Make sure python3.8 or higher is installed.
5. Make sure emboss, bwa, bowtie2, muscle and samtools are installed. (For example via conda).
6. Install [`pacbiotools`](https://github.com/jweile/pacbiotools)
7. Install [`barseqPro`](https://github.com/jweile/barseqPro) (this includes `pacybara`)


## Running pacybara
1. 	Prepare a reference fasta file for the amplicon (let's call it `amplicon.fasta` here) and note the start and end positions of the ORF within. Also create a directory that will contain the output, for example `outputDir/`. Also take note of the barcode degeneracy code(s) in the reference file. By default, pacybara assumes barcodes of length 25 with a strong-weak nucleotide pattern (`SWSWSWSWSWSWSWSWSWSWSWSWS`). If your barcoding strategy differs from this, you will need to indicate this when you run pacybara.
2. 	Let's assume your reference fasta file is called `amplicon.fasta` and your ORF within starts at position 207 and ends at 4604. Let's also assume you're using the default barcode pattern as above. The Execute pacybara on the sample: `pacybara.sh --orfStart 207 --orfEnd 4604  reads_demuxed.bc1001-bc1001_RQ998.fastq.gz amplicon.fasta outputDir/`. However, it is recommended that you run pacybara as a SLURM job. Using clusterutil, you can do so like this: `submitjob.sh -n myPacybaraJob -c 12 -m 24G -t 36:00:00 -l pacybara.log -e pacybara.log -- pacybara.sh --orfStart 207 --orfEnd 4604  reads_demuxed.bc1001-bc1001_RQ998.fastq.gz amplicon.fasta outputDir/`. As you can see we requested 12CPU cores, 24GB of RAM and 36hours runtime for the job. You can modify this as needed. 

Here's the full usage information on the pacybara executable: 
```
pacybara.sh [-b|--barcode <BARCODE>] [-s|--orfStart <ORFSTART>] 
   [-e|--orfEnd <ORFEND>] [--minQual <MINQUAL>] 
   [-m|--minMatches <MINMATCHES>] [--maxDiff <MAXDIFF>] 
   [-j|--minJaccard <MINJACCARD>] [-v|--virtualBC]
   [-c|--cpus <NUMBER>] [-q|--queue <QUEUE>] 
   [--blacklist {<NODE>,}]
   <INFASTQ> <FASTA> [<WORKSPACE>]

-b|--barcode   : The barcode degeneracy code sequence, defaults to
                 $BARCODE
-s|--orfStart  : The ORF start position, defaults to $ORFSTART
-e|--orfEnd    : The ORF end position, defaults to $ORFEND
--minQual      : The minimum PHRED quality for variant basecall to be
                 considered real. Defaults to $MINQUAL
-m|--minMatches: The minimum number of variant matches for a merge
                 to occur. Defaults to $MINMATCHES
--maxDiff      : The maxium allowed edit distance between two clusters
                 for a merge to occur. Defaults to $MAXDIFF
-j|--minJaccard: The minimum Jaccard coefficient between to clusters
                 for a merge to occur. Defaults to $MINJACCARD
-v|--virtualBC : Use virtual barcodes (fusion of up- and down-tags) 
                 for clustering. Otherwise only use uptags.
-d|--downTag   : Cluster based on second barcode (i.e. "down-tag") 
                 instead of first or virtual barcode
-c|--cpus      : Number of CPUs to use, defaults to $THREADS
-q|--queue     : The queue (slurm partition) to use
--blacklist    : A comma-separated list of compute nodes not to be used
<INFASTQ>      : The input fastq.gz file to process
<FASTA>        : The raw reference fasta file 
<WORKSPACE>    : The workspace directory. Defaults to 'workspace/'
```

## Pacybara Output
The output directory will contain multiple items: 
  * A bam file of all the reads aligned against the amplicon
  * A tarball containing the log files of all the parallel alignments
  * A directory ending in `_extract`. This contains fastq files of the extracted barcodes for each read and a file called `genotypes.csv.gz` which contains the extracted ORF genotypes for each read.
  * A directory ending in `_clustering`. This contains a number of intermediate files, as well as the final `clusters_transl.csv.gz` and `clusters_transl_filtered.csv.gz`
	  * `clusters_transl.csv.gz` contains an unfiltered table of all final clusters (i.e. clones) and barcodes and genotypes.
	  * `clusters_transl_filtered.csv.gz` is the same, but filtered for clones supported by >1 CCS read and no barcode collisions.
    * `clusters_transl_softfilter.csv.gz` is filtered less strictly. It still requires >1 CCS read, but allows clones from barcode collisions as long as they dominate that barcode with at least a 2/3 majority.
  * Within the `*_clustering` directory you will also find a `qc/` subdirectory. This contains a number of QC plots.

### Converting the Pacybara output for use with barseqPro.
The barseqPro software requires a library of barcode associations. This table can be made using any of the `clusters_transl*` tables, depending on the desired filtering level. I recommend using the `softfilter` version for best results. To convert it to the required format, extract it using gunzip (or similar software) and open it in a spreadsheet software such as `Calc` or `Excel`. Then:
1. Delete the following columns: `virtualBarcode`, `reads`, `collisions`, `upTagCollisions`, and `geno`.
2. Rename the `upBarcode` column to `barcode`.
3. Move the `size` column to the end (so that is now column `I`)

You should now have the following 9 columns, in the following order:
`barcode`, `hgvsc`, `hgvsp`, `codonChanges`, `codonHGVS`, `aaChanges`, `aaChangeHGVS`, `offTarget`, `size`. Save it as a CSV file, for example `myBarcodeLibrary.csv`.

## Running barseqPro
1. Prepare a parameter sheet for your run. An example can be found below.
2. Prepare a folder containing the FASTQ files for your samples. Make sure the names of the FASTQ files match the sample names in the parameter sheet.
3. The results will be written to the current working directory, so I recommend creating a new directory for this purpose first; e.g. `mkdir barseq_MFG_2022-08-15 && cd $_`.
4. BarseqPro just takes two arguments, the path to the FASTQ folder and the parameter sheet, e.g. `barseq.sh fastqFolder/ parameters.txt`. Again, it is recommended to submit this as a slurm job: `submitjob.sh -n myBarseqJob -c 12 -m 24G -t 36:00:00 -l barseq.log -e barseq.log -- barseq.sh fastqFolder/ parameters.txt`. 

Runtime is very dependent on the maxError parameter in the parameter sheet. Allowing for only 1 error is much faster than 2 or 3!

### BarseqPro Output:

BarseqPro will create the following folders: `counts/`, `logs/`, and `scores/`
If some of these are missing or empty, you can check the log files for errors. The most important results will be in the `scores/` folder:
  * joint_scores_*.csv will contain the map scores for each amino acid change in a given assay in MaveDB format.
  * all_aa_scores_*.csv is similar, but contains a more detailed breakdown of scores based on single mutants only, multi-mutant averaging, and inverse multiplicative model inference.
  * allScores.csv contains fitness scores for each individual barcoded clone.

### Running barseq QC
You can run `barseq_qc.R` on the output of barseqPro. It takes the following arguments: The path to the `allLRs.csv` file, the path to the `allCounts.csv` file, and the desired output directory, e.g.: `barseq_qc.R scores/allLRs.csv counts/allCounts.csv qc/`

It will generate a number of QC plots in the specified output folder.

### Example parameter sheet
A valid parameter sheet is a simple plain-text .TXT file with four sections:
`ARGUMENTS`, `FLANKING SEQUENCES`, `CODING SEQUENCE`, and `SAMPLE TABLE`. Each of them are demarkated using `#BEGIN` and `#END` statements. 

The `ARGUMENTS` section is a bash script that needs to define the following variables: `TITLE`, `LIBRARY`, `BCLEN`, `BCMAXERR`, `REVCOMP`, and `PAIREDEND`. See the example below for individual explanations. 

The `FLANKING SEQUENCES` section should follow FASTA format and contain two entries: `upstream` and `downstream`. They denote the upstream and downstream sequences of the barcode.

The `CODING SEQUENCE` section should also follow FASTA format and contain the subject coding sequence / open reading frame.

Finally, the `SAMPLE TABLE` should contain a tab-separated table with the following columns: `sample`, `assay`, `condition`, and `replicate`. Here the `sample` should correspond to the label of the corresponding FASTQ file.

```
#########################
#BarseqPro parameter file
#########################

#FASTQ location: $HOME/projects/barseqPro/LDLR_R5_downtag_fastq

#BEGIN ARGUMENTS
#Experiment title
TITLE="LDLR-pacybara-downtag-R05"
#Barcode library CSV table (from Pacbio pipeline)
LIBRARY="$HOME/projects/barseqPro/libraries/LDLR_R05_downtag.csv"
#Barcode length:
BCLEN=25
#Maximum number of errors (mismatches) allowed in barcode
BCMAXERR=1
#Whether reads are in reverse complement relative to
# the barcodes defined in the library table. 1=yes, 0=no
REVCOMP=0
#Whether the sequencing run uses paired-end mode. 1=yes, 0=no
PAIREDEND=0
#Frequency filter cutoff. The minimum relative frequency required in the nonselect (i.e. "all") condition.
FREQFILTER=5e-7
#Bottleneck filter cutoff. The minimum read count required in any select condition. 0 = off
BNFILTER=0
#END ARGUMENTS

#BEGIN FLANKING SEQUENCES
>upstream
TGTGAAGG
>downstream
CCTCAGTC
#END FLANKING SEQUENCES

#BEGIN CODING SEQUENCE
>CDS
ATGGGGCCCTGGGGCTGGAAATTGCGCTGGACCGTCGCCTTGCTCCTCGCCGCGGCGGGGACTGCAGTG
GGCGACAGATGTGAAAGAAACGAGTTCCAGTGCCAAGACGGGAAATGCATCTCCTACAAGTGGGTCTGC
GATGGCAGCGCTGAGTGCCAGGATGGCTCTGATGAGTCCCAGGAGACGTGCTTGTCTGTCACCTGCAAA
TCCGGGGACTTCAGCTGTGGGGGCCGTGTCAACCGCTGCATTCCTCAGTTCTGGAGGTGCGATGGCCAA
GTGGACTGCGACAACGGCTCAGACGAGCAAGGCTGTCCCCCCAAGACGTGCTCCCAGGACGAGTTTCGC
TGCCACGATGGGAAGTGCATCTCTCGGCAGTTCGTCTGTGACTCAGACCGGGACTGCTTGGACGGCTCA
GACGAGGCCTCCTGCCCGGTGCTCACCTGTGGTCCCGCCAGCTTCCAGTGCAACAGCTCCACCTGCATC
CCCCAGCTGTGGGCCTGCGACAACGACCCCGACTGCGAAGATGGCTCGGATGAGTGGCCGCAGCGCTGT
AGGGGTCTTTACGTGTTCCAAGGGGACAGTAGCCCCTGCTCGGCCTTCGAGTTCCACTGCCTAAGTGGC
GAGTGCATCCACTCCAGCTGGCGCTGTGATGGTGGCCCCGACTGCAAGGACAAATCTGACGAGGAAAAC
TGCGCTGTGGCCACCTGTCGCCCTGACGAATTCCAGTGCTCTGATGGAAACTGCATCCATGGCAGCCGG
CAGTGTGACCGGGAATATGACTGCAAGGACATGAGCGATGAAGTTGGCTGCGTTAATGTGACACTCTGC
GAGGGACCCAACAAGTTCAAGTGTCACAGCGGCGAATGCATCACCCTGGACAAAGTCTGCAACATGGCT
AGAGACTGCCGGGACTGGTCAGATGAACCCATCAAAGAGTGCGGGACCAACGAATGCTTGGACAACAAC
GGCGGCTGTTCCCACGTCTGCAATGACCTTAAGATCGGCTACGAGTGCCTGTGCCCCGACGGCTTCCAG
CTGGTGGCCCAGCGAAGATGCGAAGATATCGATGAGTGTCAGGATCCCGACACCTGCAGCCAGCTCTGC
GTGAACCTGGAGGGTGGCTACAAGTGCCAGTGTGAGGAAGGCTTCCAGCTGGACCCCCACACGAAGGCC
TGCAAGGCTGTGGGCTCCATCGCCTACCTCTTCTTCACCAACCGGCACGAGGTCAGGAAGATGACGCTG
GACCGGAGCGAGTACACCAGCCTCATCCCCAACCTGAGGAACGTGGTCGCTCTGGACACGGAGGTGGCC
AGCAATAGAATCTACTGGTCTGACCTGTCCCAGAGAATGATCTGCAGCACCCAGCTTGACAGAGCCCAC
GGCGTCTCTTCCTATGACACCGTCATCAGCAGGGACATCCAGGCCCCCGACGGGCTGGCTGTGGACTGG
ATCCACAGCAACATCTACTGGACCGACTCTGTCCTGGGCACTGTCTCTGTTGCGGATACCAAGGGCGTG
AAGAGGAAAACGTTATTCAGGGAGAACGGCTCCAAGCCAAGGGCCATCGTGGTGGATCCTGTTCATGGC
TTCATGTACTGGACTGACTGGGGAACTCCCGCCAAGATCAAGAAAGGGGGCCTGAATGGTGTGGACATC
TACTCGCTGGTGACTGAAAACATTCAGTGGCCCAATGGCATCACCCTAGATCTCCTCAGTGGCCGCCTC
TACTGGGTTGACTCCAAACTTCACTCCATCTCAAGCATCGATGTCAACGGGGGCAACCGGAAGACCATC
TTGGAGGATGAAAAGAGGCTGGCCCACCCCTTCTCCTTGGCCGTCTTTGAGGACAAAGTATTTTGGACA
GATATCATCAACGAAGCCATTTTCAGTGCCAACCGCCTCACAGGTTCCGATGTCAACTTGTTGGCTGAA
AACCTACTGTCCCCAGAGGATATGGTCCTCTTCCACAACCTCACCCAGCCAAGAGGAGTGAACTGGTGT
GAGAGGACCACCCTGAGCAATGGCGGCTGCCAGTATCTGTGCCTCCCTGCCCCGCAGATCAACCCCCAC
TCGCCCAAGTTTACCTGCGCCTGCCCGGACGGCATGCTGCTGGCCAGGGACATGAGGAGCTGCCTCACA
GAGGCTGAGGCTGCAGTGGCCACCCAGGAGACATCCACCGTCAGGCTAAAGGTCAGCTCCACAGCCGTA
AGGACACAGCACACAACCACCCGGCCTGTTCCCGACACCTCCCGGCTGCCTGGGGCCACCCCTGGGCTC
ACCACGGTGGAGATAGTGACAATGTCTCACCAAGCTCTGGGCGACGTTGCTGGCAGAGGAAATGAGAAG
AAGCCCAGTAGCGTGAGGGCTCTGTCCATTGTCCTCCCCATCGTGCTCCTCGTCTTCCTTTGCCTGGGG
GTCTTCCTTCTATGGAAGAACTGGCGGCTTAAGAACATCAACAGCATCAACTTTGACAACCCCGTCTAT
CAGAAGACCACAGAGGATGAGGTCCACATTTGCCACAACCAGGACGGCTACAGCTACCCCTCGAGACAG
ATGGTCAGTCTGGAGGATGACGTGGCGTGA
#END CODING SEQUENCE

#BEGIN SAMPLE TABLE     
sample  assay condition replicate
LDLR_Reg5_Rep1_All  Uptake  All 1
LDLR_Reg5_Rep2_All  Uptake  All 2
LDLR_Reg5_Rep3_All  Uptake  All 3
LDLR_Reg5_Rep4_All  Uptake  All 4
LDLR_Reg5_Rep1_F5 Uptake  F5  1
LDLR_Reg5_Rep2_F5 Uptake  F5  2
LDLR_Reg5_Rep3_F5 Uptake  F5  3
LDLR_Reg5_Rep4_F5 Uptake  F5  4
#END SAMPLE TABLE
```
