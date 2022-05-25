## Pre-requisites and installation:
1. Pacybara and barseqPro require a SLURM HPC cluster.
2. Install [`clusterutil`](https://github.com/jweile/clusterutil) 
3. Make sure the following R version 4 or higher and the following R packages are installed: `yogitools`, `yogiseq`, `hgvsParseR`,`argparser`,`pbmcapply`,`hash`,`bitops`
	* yogitools, yogiutil and hgvsParser can be installed via `remotes::install_github("jweile/yogitools")`, etc
4. Make sure emboss, bwa, bowtie2 and samtools are installed.
4. Install [`pacbiotools`](https://github.com/jweile/pacbiotools)
5. Install [`barseqPro`](https://github.com/jweile/barseqPro) (this includes `pacybara`)


## Running pacybara
1. 	Prepare a reference fasta file for the amplicon (let's call it `amplicon.fasta` here) and note the start and end positions of the ORF within. Also create a directory that will contain the output, for example `outputDir/`
2. 	Execute pacybara on the sample: `pacybara.sh --orfStart 207 --orfEnd 4604  reads_demuxed.bc1001-bc1001_RQ998.fastq.gz amplicon.fasta outputDir/`
  * Use SLURM to submit pacybara as its own job. You can use the submitjob.sh script to do so.
	* Pacybara can take a very long time to run. Runtime roughly increases with the of the number of reads ***squared***! That means 300K reads may take 3-4 hours; 1M reads make take 1.5 days; and 2M reads may take over a week! Make sure to set the time cutoff accordingly during job submission.
3. 	Run pacybara_qc.R to generate a QC report.

## Pacybara Output
The output folder will contain multiple items: 
  * A bam file of all the reads aligned against the amplicon
  * A tarball containing the log files of all the parallel alignments
  * A folder ending in `extract`. This contains fastq files of the extracted barcodes for each read and a file called `genotypes.csv.gz` which contains the extracted ORF genotypes for each read.
  * A folder ending in `clustering`. This contains a number of intermediate files, as well as the final `clusters_transl.csv.gz` and `clusters_transl_filtered.csv.gz`
	  * `clusters_transl.csv.gz` contains an unfiltered table of all final clusters (i.e. clones) and barcodes and genotypes.
	  * `clusters_transl.csv.gz` is the same, but filtered for clones supported by >1 CCS read and no barcode collisions.

## Running barseqPro
1. Prepare a parameter sheet for your run. Examples can be found in the `templates` folder.
2. Prepare a folder containing the FASTQ files for your samples. Make sure the names of the FASTQ files match the sample names in the parameter sheet.
3. Execute barseqPro: `barseqPro.sh fastqFolder/ parameters.txt`
  * Use SLURM to submit barseqPro as its own job. You can use the submitjob.sh script to do so.
  * Runtime is very dependent on the maxError parameter in the parameter sheet. Allowing for only 1 error is much faster!
