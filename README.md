# RiboProAnalysis

**RiboProAnalysis** is a pipeline for Ribosome Profiling analysis af any eukaryotic genome from Ensembl 75+.
It performs all the neccessry pre-processing steps (quality control, filtering, trimming and size selection), reads mapping to rRNA and reference genome, counting on CDS for each gene and differential analysis from raw Ribosome Profiling data.

## Usage
RiboProAnalysis can be used either via a Docker image (_URL!_) or a standard Bash script with several cases: it can performs demultiplexing on multiplexed FASTQ (reads MUST begin with the index sequence) and use of RNA-seq counts to give a study of the mode of regulation of the translation.
If you use FASTQ files (no demultiplexing), the extension have to be .fastq

A configuration file .conf is mandatory to launch the pipeline.
If there is no use of RNA-seq counts, a tabulated design file -named target.txt- is needed.
The user have to build rRNA and genome index files before start running the pipeline.
If you have RNA-seq counts files they must be named: SAMPLENAME_mRNAcounts.txt for counts file and SAMPLENAME_mRNA.transcriptome.mapping.bam for mapping to transcriptome BAM files.

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
  * Run RiboProAnalysis docker container with the following command in the working directory:
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
  * Run RiboProAnalysis bash program with following command in the working directory :
  ```
  riboproanalysis.sh MyConfigurationFile.conf
  ```

### Variables to set in the configuration file

Variables | Explanation | Choices/Examples | Default
----------|-------------|------------------|--------
PATH_TO_GENOME_INDEX | Absolute path to genome index previously built with STAR | /absolute/path/to/genome/index | Mandatory if not Docker mode
PATH_TO_rRNA_INDEX | Absolute path to rRNA index previously built with Bowtie1 | /absolute/path/to/rRNA/index | Mandatory if not Docker mode

## Installation
This software could be launched from a Docker container launching Docker containers itself, or from a Bash script launching Docker containers.

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

* Pull RiboProAnalysis image.

## Input files

### Configuration file
You have to create your configuration file .conf in the working directory. It is a little Bash script which is imported in the main Bash script.
You put mandatory and interesting variables presented in **Available variables to set in the configuration file**.

The syntax to declare a variable is :
```
export VARIABLE_NAME=MyVariable
```

#### Example of configuration file for run with the Bash script
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

#### Example of configuration file for run with Docker container
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
