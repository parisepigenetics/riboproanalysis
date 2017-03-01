# RiboProAnalysis

**RiboProAnalysis** is a pipeline for Ribosome Profiling analysis of any eukaryotic genome from Ensembl 75+ or later.
It performs all the necessary pre-processing steps (quality control, filtering, trimming and size selection), filter reads mapped to rRNA, map to reference genome, counting on CDS for each gene and differential analysis from raw Ribosome Profiling data.

## Usage
RiboProAnalysis can be used either via a Docker image (_URL_HERE!!!_) or a standard Bash script with several options. It can be used to perform demultiplexing on multiplexed FASTQ (reads MUST begin with the index sequence). Availability of RNA-seq count data allows comparison of transcript abundance with riboseq counts, providing a proxy to analyse the regulation level during differential analysis.
If you use FASTQ files (no demultiplexing), the extension has to be .fastq

A configuration file .conf is mandatory to launch the pipeline.
Without RNA-seq count data, a tabulated design file -named target.txt- is needed.
The user have to build rRNA and genome index files before start running the pipeline.
If available, RNA-seq count files must be named: SAMPLENAME_mRNAcounts.txt for counts file and SAMPLENAME_mRNA.transcriptome.mapping.bam for mapping to transcriptome BAM files.

### Steps to run the pipeline
* Build a bowtie index for rRNA sequences (use Bowtie1):
```
bowtie-build rRNA.fasta rRNA
```
* Build a STAR index for reference genome:
```
STAR --runMode genomeGenerate --genomeDir /path/to/genome/index --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ... --sjdbGTFfile /path/to/gtf/annotations \
--sjdbOverhang 28
```
<!--- (TODO I think that the script is already doing that)
* Create a tmp/ directory in your working directory with the command:
```
mkdir tmp/
```--->
* Run the pipeline (there are two options):
  * Run the RiboProAnalysis docker container (docker mode) with the following command in the working directory:
  ```
  docker run --rm --privileged --name ribopro -v /var/run/docker.sock:/var/run/docker.sock -v $(pwd):/home -w /home \
  -v /etc/passwd:/etc/passwd
  -v /path/to/rRNA/index:/rRNAindexdirectory \
  -v /path/to/genome/index:/genomeindexdirectory \
  -v /path/to/directory/containig/genome/fasta/file:/genomefastafiledirectory \
  -v /path/to/directory/containing/transcriptome/fasta/file:/transcriptomedirectory \
  -v /path/to/directory/containig/GTF/Ensembl/annotations:/root \
  -v $(pwd)/tmp:/tmp \
  parisepigenetics/riboproanalysis bash -c "riboproanalysisDocker.sh My_configuration_file.conf"
  ```
  * Run the RiboProAnalysis bash script (local mode) with the following command in the working directory :
  ```
  riboproanalysis.sh MyConfigurationFile.conf
  ```

### Variables to set in the configuration file

Variables | Explanation |  Choices-Examples  | Default
----------|-------------|--------------------|--------
PATH_TO_GENOME_INDEX | Absolute path to genome index previously built with STAR | /absolute/path/to/genome/index | Mandatory for the local mode
PATH_TO_rRNA_INDEX | Absolute path to rRNA index previously built with Bowtie1 | /absolute/path/to/rRNA/index | Mandatory for the local mode
PATH_TO_ANNOTATION_FILE | Absolute path to GTF annotations file (Ensembl 75+) | /absolute/path/to/gtf/annotations | Mandatory
PATH_TO_REFERENCE_TRANSCRIPTOME_FILE | Absolute path to transcriptome FASTA file (based on genome fasta file from Ensembl 75+) | /absolute/path/to/transcriptome/fasta/file | Mandatory with RNA-seq data
PATH_TO_REFERENCE_GENOME_FILE | Absolute file to genome FASTA file from Ensembl 75+ | /absolute/path/to/genome/fasta/file | Mandatory
USER_IDS | Result of the bash command : $(id -u):$(id -g) | UserId:GroupId | Mandatory for Docker mode
SAMPLE_ARRAY | Array containing sample names (if demultiplexing) or FASTQ file for each sample | (Sample1 Sample 2 Sample 3) OR (Samp1.fastq Samp2.fastq Samp3.fastq) | Mandatory
CONDITION_ARRAY | Array containig condition name of each sample respecting the same order | (Cond_Samp1 Cond_Samp2 Cond_Samp3) | Mandatory
ADAPTER_SEQUENCE_THREE_PRIME | Adapter sequence for 3' trimming | AAAAAAAGGTCCTAA | Mandatory
STRANDED | Answer for stranded option of HTSeq-Count | yes/no/reverse | Mandatory
PATH_TO_RAW_UNDEMULTIPLEXED_FILE | Absolute path to multiplexed FASTQ file | /absolute/path/to/multiplexed/fastq | Mandatory for demultiplexing
SAMPLE_INDEX_ARRAY | Array containing 5' index used for demultiplexing. Must be in the same order as in SAMPLE_ARRAY, index match with respective sample name | (IndexSamp1 IndexSamp2 IndexSamp3) | Mandatory for demultiplexing
ANSWER_REMOVE_POLYN_READS | Parameter to remove reads containing more than 2 N bases | YES / NO | NO
ANSWER_DEMULTIPLEXING | Parameter to launch demultiplexing step | YES / NO | NO
ANSWER_REMOVE_PCR_DUPLICATES | Parameter to launch PCR duplicates removing | YES / NO | NO
ANSWER_RNASEQ_DATA | Parameter to launch differential analysis with or without RNA-seq data | YES / NO | NO
ANSWER_PE_RNASEQ | Parameter to know if RNA-seq data (if provided) are paired-end or not | YES / NO | NO
ANSWER_PSITE_CORRECTION | Parameter to know if the script of p-site offset auto-correction will ba used or not | YES / NO | YES
RNASEQ_LIBTYPE | Parameter to know the orientation of RNA-seq data (if provided) | reversestrand/forwardstrand/unstranded | Mandatory for RNA-seq data
STOP_EXEC_PSITE_CORRECTION | Parameter to stop the pipeline to correct p-site offset file by hand | YES / NO | NO
DIFFERENTIAL_ANALYSIS_PACKAGE | Choice of the R package used by SARTools | DESEQ2 / EDGER | EDGER
AUTHOR | Author's name | Alexandra | Mandatory for SARTools
REFERENCE_CONDITION | Reference condition for the statistical analysis of SARTools | WT | Mandatory for SARTools
CHECK_DOCKER_IMAGES | Check the tags of Docker images | YES / NO | NO


## Installation
This software could be launched from a docker container launching docker containers itself (docker mode), or from a bash script launching docker containers (bash mode).

You should:
* Install Docker on your machine.
* Pull the following docker images from the [Genomic Paris Centre](https://github.com/GenomicParisCentre) dockerfiles public repository on github:
	* genomicpariscentre/fastqc:0.11.5
	* genomicpariscentre/cutadapt:1.8.3
	* genomicpariscentre/bowtie1:1.1.1
	* genomicpariscentre/star:2.5.1b
	* genomicpariscentre/gff3-ptools:0.4.0
	* genomicpariscentre/samtools:0.1.19
	* genomicpariscentre/htseq:0.6.1p1
	* genomicpariscentre/babel:0.3-0
	* genomicpariscentre/sartools:1.3.2
	* genomicpariscentre/bcbio-nextgen:1.0.0a0
	* parisepigenetics/plastid:0.4.6
	* genomicpariscentre/ribomap:1.2

* Pull this Github RiboProAnalysis package.


## Input files

### Configuration file
You have to create your configuration file .conf in the working directory. It is a little Bash script which is imported in the main Bash script.
You put mandatory and interesting variables presented in **Available variables to set in the configuration file**.

The syntax to declare a variable is :
```
export VARIABLE_NAME=MyVariable
```

#### Example configuration file to run with the Bash script
```
export PATH_TO_RAW_UNDEMULTIPLEXED_FILE=/import/disir01/bioinfo/RiboPro/Riboprotma_project/2015_240_NoIndex_L008_R1_001.fastq
export PATH_TO_GENOME_INDEX=/import/disir01/bioinfo/RiboPro/IndexAlignement/STAR/yeastGenomeEnsembl
export PATH_TO_rRNA_INDEX=/import/disir01/bioinfo/RiboPro/IndexAlignement/Bowtie1/rRNALevureNCBI
export PATH_TO_ANNOTATION_FILE=/import/rhodos01/shares-net/bioinfo/RiboPro/FichiersLevure/GenomeAnnotations_Ensembl/Saccharomyces_cerevisiae.R64-1-1.75.gtf
export PATH_TO_REFERENCE_TRANSCRIPTOME_FILE=/absolute/path/to/transcriptome/fasta/file.fasta
export PATH_TO_REFERENCE_GENOME_FILE=/absolute/path/to/genome/fasta/file.fasta
export ANSWER_DEMULTIPLEXING=YES
export ANSWER_REMOVE_PCR_DUPLICATES=YES
export ANSWER_RNASEQ_DATA=NO
export DIFFERENTIAL_ANALYSIS_PACKAGE=EDGER
export SAMPLE_ARRAY=(RT1 RT11 RT7 RT13)
export CONDITION_ARRAY=(RT RT RD RD)
export SAMPLE_INDEX_ARRAY=(NNNGGTTNN NNNAACCNN NNNTTAGNN NNNCGGANN)
export ADAPTER_SEQUENCE_THREE_PRIME=AGATCGGAAGAGCGGTTCAG
export STRANDED=yes
export AUTHOR=User
export REFERENCE_CONDITION=WildType
```

#### Example configuration file to run with Docker container
```
export USER_IDS=2747:100
export PATH_TO_RAW_UNDEMULTIPLEXED_FILE=/home/2015_240_NoIndex_L008_R1_001.fastq
export PATH_TO_ANNOTATION_FILE=/root/Saccharomyces_cerevisiae.R64-1-1.75.gtf
export PATH_TO_REFERENCE_TRANSCRIPTOME_FILE=/transcriptomedirectory/transcrioptome.fasta
export PATH_TO_REFERENCE_GENOME_FILE=/genomefastafiledirectory/genome.fasta
export ANSWER_DEMULTIPLEXING=NO
export ANSWER_REMOVE_PCR_DUPLICATES=YES
export ANSWER_RNASEQ_DATA=YES
export SAMPLE_ARRAY=(RT1.fastq RT2.fastq RD1.fastq RD2.fastq)
export CONDITION_ARRAY=(RT RT RD RD)
export SAMPLE_INDEX_ARRAY=(NA)
export ADAPTER_SEQUENCE_THREE_PRIME=AGATCGGAAGAGCGGTTCAG
export STRANDED=yes
```
<!---TODO: perhaps we can upload the yeast data so that people can run the exemplary analysis--->
