#!/bin/bash -ex

#######################################################################
## The main bash script to run the ribosome profile pipeline.        ##
##                                                                   ##
## Version 0.9                                                       ##
## Developer: Alexandra Bomane                                       ##
##            <alexandra.bomane@laposte.fr>                          ##
## Maintainer: Costas Bouyioukos                                     ##
##      email: costas.bouyioukos@univ-paris-diderot.fr               ##
## Developed in the PTER team at UMR7216, 2016-18                    ##
## This program is a free software released under the GNU GPL v. 3   ##
#######################################################################

######################### Variables section ###########################
## Environment

# For debugging
#set -xv
# Allow to stop the program after an error, BUT doesn't display the error
# TODO fix this, display the error. Look at the shebang line.

#### TODO
# Integrate parameter CASAVA_VERSION --> if 1.8 --> Run remove bad IQF (Illumina Qiality Filter) ; if not 1.8 --> DON'T run remove bad IQF (Illumina Qiality Filter)
# Integrate parameter ANSWER_SPLICING_JUNCTIONS --> if YES let the STAR command line for mapping has it is ; if NO --> DON'T do splcing using --alignIntronMax 1 parameter of STAR

# Working directory
export WORKDIR=$(pwd)

# Default variables
export SAMPLE_INDEX_ARRAY=(NONE)
export ANSWER_REMOVE_POLYN_READS=NO
export ANSWER_DEMULTIPLEXING=NO
export ANSWER_REMOVE_PCR_DUPLICATES=NO
export ANSWER_RNASEQ_DATA=NO
#ANSWER_KEEP_MULTIREAD=NO
export DIFFERENTIAL_ANALYSIS_PACKAGE=EDGER
export CHECK_DOCKER_IMAGES=NO
export ANSWER_PSITE_CORRECTION=YES
export STOP_EXEC_PSITE_CORRECTION=NO

# Import configuration (.conf) file specified by the user: it overides default variables
source $1

# Error exit function
function error_exit {
  PROGNAME=$(basename $0)
  echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}
export -f error_exit

# Check ANSWER_* variables
export WORKING_ANSWER_REMOVE_POLYN_READS=${ANSWER_REMOVE_POLYN_READS^^}
if [ ! $WORKING_ANSWER_REMOVE_POLYN_READS = NO ] && [ ! $WORKING_ANSWER_REMOVE_POLYN_READS = YES ]
then
  echo "Check your ANSWER_REMOVE_POLYN_READS parameter. It must be YES or NO."
  error_exit "$LINENO: An error has occurred."
fi

export WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
if [ ! $WORKING_ANSWER_DEMULTIPLEXING = NO ] && [ ! $WORKING_ANSWER_DEMULTIPLEXING = YES ]
then
  echo "Check your ANSWER_DEMULTIPLEXING parameter. It must be YES or NO."
  error_exit "$LINENO: An error has occurred."
fi

export WORKING_ANSWER_REMOVE_PCR_DUPLICATES=${ANSWER_REMOVE_PCR_DUPLICATES^^}
if [ ! $WORKING_ANSWER_REMOVE_PCR_DUPLICATES = NO ] && [ ! $WORKING_ANSWER_REMOVE_PCR_DUPLICATES = YES ]
then
  echo "Check your ANSWER_REMOVE_PCR_DUPLICATES parameter. It must be YES or NO."
  error_exit "$LINENO: An error has occurred."
fi
export WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
if [ ! $WORKING_ANSWER_RNASEQ_DATA = NO ] && [ ! $WORKING_ANSWER_RNASEQ_DATA = YES ]
then
  echo "Check your ANSWER_RNASEQ_DATA parameter. It must be YES or NO."
  error_exit "$LINENO: An error has occurred."
fi

if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
then
  WORKING_ANSWER_PE_RNASEQ=${ANSWER_PE_RNASEQ^^}
  if [ ! $WORKING_ANSWER_PE_RNASEQ = NO ] && [ ! $WORKING_ANSWER_PE_RNASEQ = YES ]
  then
    echo "Check your ANSWER_PE_RNASEQ parameter. It must be YES or NO."
    error_exit "$LINENO: An error has occurred."

  elif [ -z $ANSWER_PE_RNASEQ ]
  then
    echo "ANSWER_PE_RNASEQ parameter is mandatory if you have RNA-seq data. It must be YES or NO."
    error_exit "$LINENO: An error has occurred."
  fi

  if [ ! $RNASEQ_LIBTYPE = "forwardstrand" ] && [ ! $RNASEQ_LIBTYPE = "unstranded" ] && [ ! $RNASEQ_LIBTYPE = "reversestrand" ]
  then
    echo "Check your RNASEQ_LIBTYPE parameter. It must be 'forwardstrand' (Read 1/single-end reads belongs to the forward strand) or 'unstranded' or 'reversestrand' (Read 1/single-end reads belongs to reverse strand)"
    error_exit "$LINENO: An error has occurred."

  elif [ -z $RNASEQ_LIBTYPE ]
  then
    echo "RNASEQ_LIBTYPE parameter is mandatory if you have RNA-seq data. It must be 'forwardstrand' (Read 1/single-end reads belongs to the forward strand) or 'unstranded' or 'reversestrand' (Read 1/single-end reads belongs to the reverse strand)"
    error_exit "$LINENO: An error has occurred."
  fi
fi

if [ "$WORKING_ANSWER_PE_RNASEQ" = YES ] && [ $RNASEQ_LIBTYPE = "reversestrand" ]
then
  SALMON_LIBTYPE="ISR"

elif [ "$WORKING_ANSWER_PE_RNASEQ" = YES ] && [ $RNASEQ_LIBTYPE = "unstranded" ]
then
  SALMON_LIBTYPE="IU"

elif [ "$WORKING_ANSWER_PE_RNASEQ" = YES ] && [ $RNASEQ_LIBTYPE = "forwardstrand" ]
then
  SALMON_LIBTYPE="ISF"

elif [ "$WORKING_ANSWER_PE_RNASEQ" = NO ] && [ $RNASEQ_LIBTYPE = "reversestrand" ]
then
  SALMON_LIBTYPE="SR"

elif [ "$WORKING_ANSWER_PE_RNASEQ" = NO ] && [ $RNASEQ_LIBTYPE = "unstranded" ]
then
  SALMON_LIBTYPE="U"

elif [ "$WORKING_ANSWER_PE_RNASEQ" = NO ] && [ $RNASEQ_LIBTYPE = "forwardstrand" ]
then
  SALMON_LIBTYPE="SF"
fi

#WORKING_ANSWER_KEEP_MULTIREAD=${ANSWER_KEEP_MULTIREAD^^}
#if [ ! $WORKING_ANSWER_KEEP_MULTIREAD = NO ] && [ ! $WORKING_ANSWER_KEEP_MULTIREAD = YES ]
#then
#  echo "Check your ANSWER_KEEP_MULTIREAD parameter. It must be YES or NO."
#  error_exit "$LINENO: An error has occurred."
#fi

if [ ! "$DIFFERENTIAL_ANALYSIS_PACKAGE" = EDGER ] && [ ! $DIFFERENTIAL_ANALYSIS_PACKAGE = DESEQ2 ]
then
  echo "Unavailable R package. Choose : EDGER or DESEQ2 (case sensitive)"
  error_exit "$LINENO: An error has occurred."
fi

export WORKING_STOP_EXEC_PSITE_CORRECTION=${STOP_EXEC_PSITE_CORRECTION^^}
if [ ! "$WORKING_STOP_EXEC_PSITE_CORRECTION" = NO ] && [ ! $WORKING_STOP_EXEC_PSITE_CORRECTION = YES ]
then
  echo "Check your STOP_EXEC_PSITE_CORRECTION parameter. It must be YES or NO."
  error_exit "$LINENO: An error has occurred."
fi


# Tmp directory #TODO check if it is nessecary.
if [ ! -e tmp/ ]
then
  mkdir -p tmp/
fi
export TMPDIR=$(readlink -f tmp/)


# Muti-mapped reads analysis directory
# TODO we will probably not include that in the final script.
if [ ! -e multi_reads_analysis/ ]
then
  mkdir -p multi_reads_analysis/
fi


## Scripts

# Main Bash script
export MAIN_SCRIPT_CANONICAL_PATH=$(readlink -f $0) ## basename $0
export CANONICAL_PATH=$(dirname $MAIN_SCRIPT_CANONICAL_PATH)

# Python and R scripts paths
export PYTHON_SCRIPTS_PATH="${CANONICAL_PATH}/PythonScripts/"
export R_SCRIPTS_PATH="${CANONICAL_PATH}/RScripts/"

#TODO probably there is NO NEED to have all these parameters as individual params. A call of the name of the script from the containing directory will be enough!
# Python scripts
export PYTHON_SCRIPT_DEMULTIPLEXING="run_demultiplexing.py"
export PYTHON_SCRIPT_REMOVE_PCR_DUP="rmDupPCR.py"
export PYTHON_SCRIPT_REMOVE_BAD_IQF="remove_bad_reads_Illumina_passing_filter.py"
export PYTHON_SCRIPT_READ_LENGTH_DISTRIBUTION="read_length_distribution.py"
export PYTHON_SCRIPT_SAM_FILTERING="sam_file_filter.py"
export PYTHON_SCRIPT_LONGEST_TRANSCRIPT="get_longest_transcripts_from_ensembl_gtf.py"
export PYTHON_SCRIPT_COUNT_CLUSTER="clusterCountings.py"
export PYTHON_SCRIPT_ADD_START_GTF="add_cds_start_to_gtf.py"
export PYTHON_SCRIPT_PSITE_AUTOCORRECTION="plastid_psite_offsets_autocorrection.py"
export PYTHON_SCRIPT_CDS_RANGE_GENERATOR="cds_range_generator.py"

# R scripts
export R_SCRIPT_BUILD_COUNTING_TABLE_RNASEQ="RNAseqCountDataMatrix.R"
export R_SCRIPT_BUILD_COUNTING_TABLE_RP="RPCountDataMatrix.R"
export R_SCRIPT_ANADIFF_BABEL="babel_RP_differentialAnalysis.R"
export R_SCRIPT_PERMT_TEST_BABEL="babel_RP_permutationTest.R"
export R_SCRIPT_ANADIFF_SARTOOLS_DESEQ2="script_DESeq2.R"
export R_SCRIPT_ANADIFF_SARTOOLS_EDGER="script_edgeR.R"


# Check mandatory parameters
if [ -z $PATH_TO_GENOME_INDEX ]
then
  echo "Give your genome index path."
  error_exit "$LINENO: An error has occurred."
fi

if [ -z "$PATH_TO_rRNA_INDEX" ]
# TODO check if we really need to do that also.
then
  echo "Give your rRNA index path."
  error_exit "$LINENO: An error has occurred."
fi

if [ -z $SAMPLE_ARRAY ]
then
  echo "Give the samples array in $1 (eg : SAMPLE_ARRAY=(sample1 sample2 ...))."
  error_exit "$LINENO: An error has occurred."
fi

if [ -z $CONDITION_ARRAY ]
then
  echo "Give the biological conditions array in $1 (eg : CONDITION_ARRAY=(sample1condition sample2condition ...))."
  error_exit "$LINENO: An error has occurred."
fi

if [ -z $ADAPTER_SEQUENCE_THREE_PRIME ]
then
  echo "Give the 3' adapter sequence in $1 (eg : ADAPTER_SEQUENCE_THREE_PRIME=MYADAPTERSEQUENCETHREEPRIME)"
  error_exit "$LINENO: An error has occurred."
fi

export WORKING_SAMPLE_ARRAY=$(echo ${SAMPLE_ARRAY[*]})

export WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
then
  if [ -z $SAMPLE_INDEX_ARRAY ]
  then
    echo "Give your sample index array in $1 (eg : SAMPLE_INDEX_ARRAY=(sample1index sample2index ...))."
    error_exit "$LINENO: An error has occurred."
  fi
fi

export WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
if [ $WORKING_ANSWER_RNASEQ_DATA = NO ]
then
  if [ -z $AUTHOR ]
  then
    $AUTHOR=UserName
  fi
fi

if [ $WORKING_ANSWER_RNASEQ_DATA = NO ]
then
  if [ -z $REFERENCE_CONDITION ]
  then
    echo "Give your reference (biological) condition in $1 (eg: REFERENCE_CONDITION=MYREFERENCECONDITION)"
    error_exit "$LINENO: An error has occurred."
  fi
fi

if [ -z $PATH_TO_REFERENCE_GENOME_FILE ]
then
  echo "Give the absolute path to your reference genome FASTA. in $1 (eg: PATH_TO_REFERENCE_GENOME_FILE=/path/to/my/ref/genome/mygenome.fasta)"
  error_exit "$LINENO: An error has occurred."
fi

if [ -z $PATH_TO_REFERENCE_TRANSCRIPTOME_FILE ]
then
  echo "Give the absolute path to your reference transcriptome FASTA. in $1 (eg: PATH_TO_REFERENCE_TRANSCRIPTOME_FILE=/path/to/my/ref/transcriptome/mytranscriptome.fasta)"
  error_exit "$LINENO: An error has occurred."
fi

if [ -z $STRANDED ]
then
  echo "Set the --stranded option of HTSeq-Count (For help : http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)"
  error_exit "$LINENO: An error has occurred."
fi

if [ -z $PATH_TO_ANNOTATION_FILE ]
then
  echo "Give the path to your GTF annotations."
  error_exit "$LINENO: An error has occurred."
fi

export WORKING_SAMPLE_INDEX_ARRAY=$(echo ${SAMPLE_INDEX_ARRAY[*]})

export SAMPLES=($(echo ${SAMPLE_ARRAY[@]%.fastq})) ## TODO No need because handling by GNU parallel (FIXME NO prtobably I made a wrong comment here.)

#WORKING_CONDITION_ARRAY=$(echo ${CONDITION_ARRAY[*]})
export WORKING_CONDITION_ARRAY=$(echo ${CONDITION_ARRAY[*]})

export SHELL=$(type -p bash)

export CONDITION_ARRAY_UNIQ=$(echo "${CONDITION_ARRAY[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')


# Check arrays length
export NB_SAMPLE=$(echo ${#SAMPLE_ARRAY[@]})

if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
then
  NB_SAMPLE_INDEX=$(echo ${#SAMPLE_INDEX_ARRAY[@]})
  if [ $NB_SAMPLE_INDEX -ne $NB_SAMPLE ]
  then
    echo "SAMPLE_INDEX_ARRAY and SAMPLE_ARRAY have different lengths. Check them."
    error_exit "$LINENO: An error has occurred."
  fi
fi

if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
then
  NB_CONDITION=$(echo ${#CONDITION_ARRAY[@]})
  if [ $NB_CONDITION -ne $NB_SAMPLE ]
  then
    echo "CONDITION_ARRAY and SAMPLE_ARRAY have different lengths. Check them."
    error_exit "$LINENO: An error has occurred."
  fi
fi


# Check Docker images
# TODO UPDATE ALL docker images and make them all available from PAris Epigenetics.
export WORKING_CHECK_DOCKER_IMAGES=${CHECK_DOCKER_IMAGES^^}
if [ ! $WORKING_CHECK_DOCKER_IMAGES = NO ]
then
  if [ $WORKING_CHECK_DOCKER_IMAGES = YES ]
  then
    docker pull genomicpariscentre/fastqc
    docker pull genomicpariscentre/cutadapt
    docker pull genomicpariscentre/bowtie
    docker pull genomicpariscentre/star
    docker pull genomicpariscentre/samtools
    docker pull genomicpariscentre/gff3-ptools
    docker pull genomicpariscentre/htseq
    docker pull genomicpariscentre/babel
    docker pull genomicpariscentre/sartools
    docker pull genomicpariscentre/ribomap
#    docker pull genomicpariscentre/bcbio-nextgen
    docker pull parisepigenetics/bcbio-nextgen
    docker pull parisepigenetics/plastid
#    docker pull parisepigenetics/riboseqr  TODO Integrate RiboSeqR
#    docker pull parisepigenetics/umitools  TODO Integrate Umitools
  else
    echo "Check your CHECK_DOCKER_IMAGES parameter. It must be YES or NO."
    error_exit "$LINENO: An error has occurred."
  fi
fi


### Tools parameters

## 3' trimming : Cutadapt
export MIN_READ_LENGTH="25"
export MAX_READ_LENGTH="45"
export FILTER_MAX_N="2"

## Prepare samples step : Seqcluster
export MIN_COUNT="2"
export MIN_SIZE="25"
export MAX_SIZE="35"
export MIN_SHARED="2"

## Align to rRNA sequences : Bowtie 1
# TODO I do not know if we need to use yet another aligner to filter the rDNA. We can use just one (either STAR or hisat2) and then remove the rDNA from the annotaions.
# Bowtie 1 Options details : -q --> Fastq file as input ; --un --> write unaligned reads to another file (.fastq) ; -S --> write hits in SAM format
export BOWTIE_OPTIONS="-q -S --un"

## Align to reference genome : STAR
export MAX_ALLOWED_MISMATCHES="2"  # alignment will be output only if it has no more mismatches than this value
export SEED_SEARCH_POINT="16"  # defines the search start point through the read - the read is split into pieces no longer than this value
export FILTER_SCORE_MIN="0"  # alignment will be output if its ratio of score to *read* length is higher than this value
export FILTER_MATCH_MIN="0.85"  # alignment will be output if its ratio of number of matched bases to *read* length is higher than this value
export MAX_LOCI_ALLOWED="1000"  # max number of loci anchors are allowed to map to
export MULTIMAP_SCORE_RANGE="0"  # the score range below the maximum score for multimapping alignments

## Multi-mapped analysis
#TODO Multi-mapped analysis will be removed from the published version.
export MULTIMAP_N_MAX="500"

## P-site offsets analysis
export MIN_LENGTH_POFFSET="25"
export MAX_LENGTH_POFFSET="35"
export DEFAULT_POFFSET="0"

## Ribomap
export MIN_FPLEN="25"
export MAX_FPLEN="35"

## HTSeq-Count
export MODE_FOR_MULTIPLE_FEATURES_READS="union"
export FEATURE_TYPE="CDS"
export IDATTR="gene_id" ## Ok for Secluster-clustering step
export FILETYPE="bam"


#### FUNCTIONS ########################################################
# TODO ALL I mean ALL the scripts we use will be dockerised (inclusing Python and R scripts.)

# Run demultiplexing to get Fastq files
# $1 = SAMPLE $2 = ADAPTER
demultiplexing() {
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
  then
    if [ -z $PATH_TO_RAW_UNDEMULTIPLEXED_FILE ]
    then
      echo "Give the path to your multiplexed FASTQ file."
      error_exit "$LINENO: An error has occurred."
    fi
    LOGFILE="$1_demultiplexing.log"
    OUTFILE="$1_demultiplex.fastq"
    if [ -s $LOGFILE ]
    then
      return 0
    else
      echo "Starting of demultiplexing :"
      $PYTHON_SCRIPT_DEMULTIPLEXING -i $PATH_TO_RAW_UNDEMULTIPLEXED_FILE -o $OUTFILE -a $2 > $LOGFILE
      if [ $? -ne 0 ]
      then
        echo "run_demultiplexing cannot run correctly ! Check your mutliplexed FASTQ path and your index adapter sequence."
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Log file : $LOGFILE generated."
      echo "End of demultiplexing."
    fi
  else
    return 0
  fi
}


# Run FastQC to check input
# $1 = directory output ; $2 = input
fastqc_quality_control() {
  if [ "$(ls -1 $1)" ]
  then
    return 0
  else
    mkdir -p $1
    if [ $? -ne 0 ]
    then
      echo "$1 image cannot be created !"
      error_exit "$LINENO: An error has occurred."
    fi
      echo "Running FastQC..."
        echo docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/fastqc -o $1 $2
        docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/fastqc -o $1 $2
      if [ $? -ne 0 ]
      then
        echo "FastQC cannot run correctly ! $@"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Finish FastQC."
  fi
}


# Run FastQC to check our demultiplexing
#TODO This function can be renamed raw_quality_control_report() ALSO this step can be merged with the previous one.
raw_quality_report() {
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
  then
    INPUT_RAW_FASTQ="$1_demultiplex.fastq"
  else
    INPUT_RAW_FASTQ="${1}.fastq"
  fi
  DIR_RAW_FASTQ_REPORT="$1_raw_fastqc_report"
  if [ -s $INPUT_RAW_FASTQ ]
  then
    fastqc_quality_control $DIR_RAW_FASTQ_REPORT $INPUT_RAW_FASTQ
  else
    echo "$INPUT_RAW_FASTQ doesn't exist ! Check your SAMPLE_ARRAY."
    error_exit "$LINENO: An error has occurred."
  fi
}


# Remove bad passing filter reads
removeBadIQF() {
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
  then
    INPUT_FASTQ="$1_demultiplex.fastq"
  else
    INPUT_FASTQ="$1.fastq"
  fi
  LOGFILE="$1_rmIQF.log"
  RM_BADIQF_OUTPUT="$1_rmIQF.fastq"
  if [ -s $LOGFILE ] && [ -s $RM_BADIQF_OUTPUT ]
  then
    echo "Filtering on Illumina Quality Filter already done."
    return 0
  else
    echo "Removing bad IQF :"
    $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_REMOVE_BAD_IQF -i $INPUT_FASTQ -o $RM_BADIQF_OUTPUT > $LOGFILE
    if [ $? -ne 0 ]
    then
      echo "Removing bad IQF cannot run correctly !"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "Log file : $LOGFILE generated"
    echo "End of removing bad IQF"
  fi
}


# Check remove bad passing filter
removeBadIQF_report() {
  RM_BADIQF_DIR="$1_rmIQF_report"
  RM_IQF_INPUT="$1_rmIQF.fastq"
  if [ -s $RM_IQF_INPUT ]
  then
    fastqc_quality_control $RM_BADIQF_DIR $RM_IQF_INPUT
  else
    echo "$RM_IQF_INPUT doesn't exist"
    error_exit "$LINENO: An error has occurred."
  fi
}


# Remove PCR duplicates --> % Amplification in log file
# TODO Can be replaced by UMI-tools
removePCRduplicates() {
  WORKING_ANSWER_REMOVE_PCR_DUPLICATES=${ANSWER_REMOVE_PCR_DUPLICATES^^}
  if [ $WORKING_ANSWER_REMOVE_PCR_DUPLICATES = "YES" ]
  then
    LOGFILE="$1_rmPCR.log"
    RM_PCRDUP_OUTPUT="$1_rmPCR.fastq"
    RM_PCRDUP_INPUT="$1_rmIQF.fastq"
    if [ -s $RM_PCRDUP_INPUT ]
    then
      if [ -s $RM_PCRDUP_OUTPUT ] && [ -s $LOGFILE ]
      then
        return 0
      else
        echo "Removing PCR duplicates :"
        awk '{ i=(NR-1) % 4; tab[i]=$0 ; if (i==3) { print tab[1]"\t"tab[0]"\t"tab[3]"\t"tab[2]} }' $RM_PCRDUP_INPUT | sort | $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_REMOVE_PCR_DUP $RM_PCRDUP_OUTPUT > $LOGFILE
        if [ ! -s $RM_PCRDUP_OUTPUT ]
        then
          echo "Cannot run rmDupPCR correctly !"
          error_exit "$LINENO: An error has occurred."
        fi
        echo "Log file : $LOGFILE generated."
        echo "End of PCR duplicates removing."
      fi
    else
      echo "You need a file which was filtered on bad Illumina Qualitiy Filter (_rmIQF.fastq) !"
      error_exit "$LINENO: An error has occurred."
    fi
  else
    return 0
  fi
}


# Run the 5' trimming
index_Adapter_trimming() {
  WORKING_ANSWER_REMOVE_PCR_DUPLICATES=${ANSWER_REMOVE_PCR_DUPLICATES^^}
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
  then
    INDEX_TRIM_OUTPUT="$1_TrimIndex.fastq"
    if [ $WORKING_ANSWER_REMOVE_PCR_DUPLICATES = YES ]
    then
      INDEX_TRIM_INPUT="$1_rmPCR.fastq"
    else
      INDEX_TRIM_INPUT="$1_rmIQF.fastq"
    fi
    INDEX_LENGTH=$(expr length $2)
    LOGFILE="$1_TrimIndex.log"
    if [ -s $INDEX_TRIM_OUTPUT ]
    then
      return 0
    else
      echo "Index adapter trimming :"
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/cutadapt bash -c "cutadapt -u $INDEX_LENGTH -o $INDEX_TRIM_OUTPUT $INDEX_TRIM_INPUT" > $LOGFILE
      if [ ! -s $INDEX_TRIM_OUTPUT ]
      then
        echo "Index adapter trimming cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Log file : $LOGFILE generated."
      echo "End of index adapter trimming."
    fi
  else
    return 0
  fi
}


# Check the 5' trimming
index_Adapter_trimming_report() {
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
  then
    DIR_INDEX_TRIM_FASTQC="$1_TrimIndex_report"
    INDEX_TRIM_INPUT="$1_TrimIndex.fastq"
    if [ -s $INDEX_TRIM_INPUT ]
    then
      fastqc_quality_control $DIR_INDEX_TRIM_FASTQC $INDEX_TRIM_INPUT
    else
      echo "$INDEX_TRIM_INPUT doesn't exist"
      error_exit "$LINENO: An error has occurred."
    fi
  fi
}


# Run Cutadapt for the 3' trimming
threePrime_trimming() {
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  WORKING_ANSWER_REMOVE_PCR_DUPLICATES=${ANSWER_REMOVE_PCR_DUPLICATES^^}
  WORKING_ANSWER_REMOVE_POLYN_READS=${ANSWER_REMOVE_POLYN_READS^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
  then
    THREEPRIME_TRIM_INPUT="$1_TrimIndex.fastq"
  else
    if [ $WORKING_ANSWER_REMOVE_PCR_DUPLICATES = "YES" ]
    then
      THREEPRIME_TRIM_INPUT="$1_rmPCR.fastq"
    else
      THREEPRIME_TRIM_INPUT="$1_rmIQF.fastq"
    fi
  fi
  THREEPRIME_TRIM_OUTPUT="$1_ThreePrime_Trim.fastq"
  LOGFILE="$1_ThreePrimeTrim.log"
  if [ -s $THREEPRIME_TRIM_OUTPUT ] && [ -s $LOGFILE ]
  then
    return 0
  else
    echo "3' trimming :"
    if [ $WORKING_ANSWER_REMOVE_POLYN_READS = "YES" ]
    then
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/cutadapt bash -c "cutadapt -a $2 --discard-untrimmed --max-n $FILTER_MAX_N -o $THREEPRIME_TRIM_OUTPUT $THREEPRIME_TRIM_INPUT > $LOGFILE"
    else
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/cutadapt bash -c "cutadapt -a $2 --discard-untrimmed -o $THREEPRIME_TRIM_OUTPUT $THREEPRIME_TRIM_INPUT > $LOGFILE"
    fi
    if [ $? -ne 0 ]
    then
      echo "Cutadapt cannot run correctly !"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "Log file : $LOGFILE generated."
    echo "End of Cutadapt."
  fi
}


# Check the 3' trimming
threePrime_trimming_report() {
  DIR_THREEPRIME_TRIM_FASTQC="$1_ThreePrime_Trim_report"
  THREEPRIME_TRIM_INPUT="$1_ThreePrime_Trim.fastq"
  if [ -s $THREEPRIME_TRIM_INPUT ]
  then
    fastqc_quality_control $DIR_THREEPRIME_TRIM_FASTQC $THREEPRIME_TRIM_INPUT
  else
    echo "$THREEPRIME_TRIM_INPUT doesn't exist"
    error_exit "$LINENO: An error has occurred."
  fi
}


# Perform size selection.
size_Selection() {
  WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
  WORKING_ANSWER_REMOVE_PCR_DUPLICATES=${ANSWER_REMOVE_PCR_DUPLICATES^^}
  if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
  then
    THREEPRIME_TRIM_INPUT="$1_ThreePrime_Trim.fastq"
    SIZE_SELECT_OUTPUT="$1_SizeSelection.fastq"
    LOGFILE="$1_SizeSelection.log"
  else
    BASENAME=$(basename $1 .f*q)
    THREEPRIME_TRIM_INPUT="${BASENAME}_ThreePrime_Trim.fastq"
    SIZE_SELECT_OUTPUT="${BASENAME}_SizeSelection.fastq"
    LOGFILE="${BASENAME}_SizeSelection.log"
  fi
  if [ -s $THREEPRIME_TRIM_OUTPUT ] && [ -s $LOGFILE ]
  then
    return 0
  else
    echo "Size selection :"
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/cutadapt bash -c "cutadapt -m $MIN_READ_LENGTH -M $MAX_READ_LENGTH -o $SIZE_SELECT_OUTPUT $THREEPRIME_TRIM_INPUT > $LOGFILE"
    if [ $? -ne 0 ]
    then
      echo "Cutadapt cannot run correctly !"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "Log file : $LOGFILE generated."
    echo "End of Cutadapt."
  fi
}


# We shake (TODO, what shake means here...??) the size selection
size_Selection_report() {
  DIR_SIZE_SELECT_FASTQC="$1_Size_Selection_report"
  SIZE_SELECT_INPUT="$1_SizeSelection.fastq"
  if [ -s $SIZE_SELECT_INPUT ]
  then
    fastqc_quality_control $DIR_SIZE_SELECT_FASTQC $SIZE_SELECT_INPUT
  else
    echo "$SIZE_SELECT_INPUT doesn't exist"
    error_exit "$LINENO: An error has occurred."
  fi
}


## Prepare Seqcluster.
#TODO look at previous comments we do not know if we need to keep seqcluster.
collapse_step_seqcluster() {
  echo "Start of multi-mapped reads analysis (by Seqcluster) :"
  echo "Collapse reads step"
  DIR_COLLAPSE_OUTPUT="multi_reads_analysis/$1_collapseStep"
  FASTQ_CLEAN_INPUT="$1_SizeSelection.fastq"
  if [ "$(ls -1 $DIR_COLLAPSE_OUTPUT)" ]
  then
    echo "Collapse reads step already done for $1"
    return 0
  else
    if [ -s $FASTQ_CLEAN_INPUT ]
    then
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v /etc/passwd:/etc/passwd -v ${WORKDIR}:/home -w /home genomicpariscentre/bcbio-nextgen bash -c "/usr/local/bin/seqcluster collapse -f $FASTQ_CLEAN_INPUT -o $DIR_COLLAPSE_OUTPUT"
      if [ $? -ne 0 ]
      then
        echo "Error during collapse reads step."
        error_exit "$LINENO: An error has occurred."
      fi
      mv ${DIR_COLLAPSE_OUTPUT}/$(basename $FASTQ_CLEAN_INPUT .fastq)_trimmed.fastq ${DIR_COLLAPSE_OUTPUT}/$(basename $FASTQ_CLEAN_INPUT .fastq)_collapsed.fastq
      echo "$DIR_COLLAPSE_OUTPUT directory was generated."
      echo "End of collapse reads step."
    else
      echo "$FASTQ_CLEAN_INPUT doesn't exist : this step cannot run correctly."
      error_exit "$LINENO: An error has occurred."
    fi
  fi
}

prepareSamples_step_seqcluster() {
  echo "Continuation of multi-mapped reads analysis (by Seqcluster) :"
  echo "Prepare samples step"
  DIR_PREPARE_SAMPLES_OUTPUT="multi_reads_analysis/prepareSamplesStep"
  if [ "$(ls -1 $DIR_PREPARE_SAMPLES_OUTPUT)" ]
  then
    echo "Prepare samples step alredy done."
    return 0
  else
    if [ -s seqclusterPrepareConfiguration.txt ]
    then
      rm seqclusterPrepareConfiguration.txt
    fi
# Building of the configuration table for Prepare samples step
    for index in ${!SAMPLE_ARRAY[*]}
    do
      WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
      if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
      then
        SAMPLE=${SAMPLE_ARRAY[$index]}
      else
        SAMPLE=$(basename ${SAMPLE_ARRAY[$index]} .fastq)
      fi
      echo -e "multi_reads_analysis/${SAMPLE}_collapseStep/${SAMPLE}_SizeSelection_collapsed.fastq"\\t$SAMPLE >> seqclusterPrepareConfiguration.txt
    done
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v /etc/passwd:/etc/passwd -v ${WORKDIR}:/home -w /home genomicpariscentre/bcbio-nextgen bash -c "/usr/local/bin/seqcluster prepare -c seqclusterPrepareConfiguration.txt -o $DIR_PREPARE_SAMPLES_OUTPUT --minc $MIN_COUNT --minl $MIN_SIZE --maxl $MAX_SIZE --min-shared $MIN_SHARED"
    if [ $? -ne 0 ]
    then
      echo "Error during prepare samples step."
      error_exit "$LINENO: An error has occurred."
    fi
    echo "$DIR_PREPARE_SAMPLES_OUTPUT directory was generated."
    echo "End of prepare samples step."
  fi
}


#TODO there is no need to employ Bowtie.
# Run Bowtie 1 to align reads to rRNA sequences : we get unmapped reads for next steps and mapped reads to have length distribution (Python script using matplotlib) # TODO Try with STAR
align_To_R_RNA() {
  for sample in ${SAMPLE_ARRAY[*]}
  do
    echo "Starting of mapping to rRNA :"
    WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
    if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
    then
      UNMAPPED_RNA_FASTQ_FILE="${sample}_no_rRNA.fastq"
      MAPPED_RNA_SAM_FILE="${sample}_rRNA_mapped.sam"
      LOGFILE_BOWTIE="${sample}_rRNA_mapping.log"
      INPUT_RNA_MAPPING="${sample}_SizeSelection.fastq"
    else
      BASENAME=$(basename $sample .fastq)
      UNMAPPED_RNA_FASTQ_FILE="${BASENAME}_no_rRNA.fastq"
      MAPPED_RNA_SAM_FILE="${BASENAME}_rRNA_mapped.sam"
      LOGFILE_BOWTIE="${BASENAME}_rRNA_mapping.log"
      INPUT_RNA_MAPPING="${BASENAME}_SizeSelection.fastq"
    fi
    rRNA_INDEX_BASENAME=$(echo $(basename ${PATH_TO_rRNA_INDEX}/*.1.ebwt | cut -f1 -d'.'))
    if [ -s $UNMAPPED_RNA_FASTQ_FILE ] && [ -s $MAPPED_RNA_SAM_FILE ]
    then
      echo "Mapping to rRNA already done for $sample"
      #return 0
    else
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home -v $PATH_TO_rRNA_INDEX:/root genomicpariscentre/bowtie1 bash -c "bowtie -p $(nproc) $BOWTIE_OPTIONS $UNMAPPED_RNA_FASTQ_FILE /root/$rRNA_INDEX_BASENAME $INPUT_RNA_MAPPING $MAPPED_RNA_SAM_FILE 2> $LOGFILE_BOWTIE"
      if [ $? -ne 0 ]
      then
        echo "Bowtie1 cannot run correctly ! Check the rRNA index path."
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Log file : $LOGFILE_BOWTIE generated."
      echo "End of mapping to rRNA."
    fi
  done
}


# Run FASTQC on unmapped fastq
unmmaped_to_rRNA_report() {
  DIR_UNMAPPED_RRNA_FASTQC="$1_no_rRNA_report"
  UNMAPPED_RRNA_INPUT="$1_no_rRNA.fastq"
  if [ -s $UNMAPPED_RRNA_INPUT ]
  then
    fastqc_quality_control $DIR_UNMAPPED_RRNA_FASTQC $UNMAPPED_RRNA_INPUT
  else
    echo "$UNMAPPED_RRNA_INPUT doesn't exist ! Check your rRNA index path."
    error_exit "$LINENO: An error has occurred."
  fi
}


# Filter rRNA from seqs.fastq generated by Seqcluster
align_To_R_RNA_seqcluster() {
  echo "Continuation of multi-mapped reads analysis (by Seqcluster) :"
  echo "Starting of mapping to rRNA for multi-mapped reads analysis :"
  DIR_RNA_MAPPING="multi_reads_analysis/mappingStep/"
  UNMAPPED_RNA_FASTQ_FILE="multi_reads_analysis/mappingStep/seqs_no_rRNA.fastq"
              MAPPED_RNA_SAM_FILE="multi_reads_analysis/mappingStep/seqs_rRNA_mapped.sam"
  UNMAPPED_RNA_LIST="multi_reads_analysis/mappingStep/seqs_no_rRNA.list"
  UNMAPPED_RNA_MATRIX="multi_reads_analysis/prepareSamplesStep/seqs_no_rRNA.ma"
  PREPARE_STEP_MATRIX="multi_reads_analysis/prepareSamplesStep/seqs.ma"
              LOGFILE="multimapped_analysis_rRNA_mapping.log"
              INPUT_RNA_MAPPING="multi_reads_analysis/prepareSamplesStep/seqs.fastq"
  if [ -e $DIR_RNA_MAPPING ]
  then
    if [ -s $UNMAPPED_RNA_FASTQ_FILE ] && [ -s $MAPPED_RNA_SAM_FILE ]
    then
      echo "Mapping to rRNA for multi-mapped reads analysis already done."
      return 0
    fi
  else
    rRNA_INDEX_BASENAME=$(echo $(basename ${PATH_TO_rRNA_INDEX}/*.1.ebwt | cut -f1 -d'.'))
    mkdir -p $DIR_RNA_MAPPING
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home -v $PATH_TO_rRNA_INDEX:/root genomicpariscentre/bowtie1 bash -c "bowtie -p $(nproc) $BOWTIE_OPTIONS $UNMAPPED_RNA_FASTQ_FILE /root/$rRNA_INDEX_BASENAME $INPUT_RNA_MAPPING $MAPPED_RNA_SAM_FILE 2> $LOGFILE"
    if [ $? -ne 0 ]
    then
      echo "Bowtie1 cannot run correctly ! Check the rRNA index path."
      error_exit "$LINENO: An error has occurred."
    fi
    if [ ! -s $UNMAPPED_RNA_MATRIX ]
    then
      grep "@seq" $UNMAPPED_RNA_FASTQ_FILE | cut -f2 -d "@" | sort -T $TMPDIR > $UNMAPPED_RNA_LIST
                        grep "^id" $PREPARE_STEP_MATRIX >> $UNMAPPED_RNA_MATRIX
      nawk 'FNR==NR{a[$0];next}($1 in a)' $UNMAPPED_RNA_LIST $PREPARE_STEP_MATRIX >> $UNMAPPED_RNA_MATRIX
    fi
    if [ $? -ne 0 ]
    then
      echo "Cannot generate $UNMAPPED_RNA_MATRIX"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "$UNMAPPED_RNA_MATRIX generated."
    echo "Log file : $LOGFILE generated."
    echo "End of mapping to rRNA (for multi-mapped analysis)."

  fi
}


# We run the Python library matplotlib TODO ...no it does not seem to do this...
mapped_to_R_RNA_distrib_length() {
  INPUT_SAM_MAPPED_RNA="$1_rRNA_mapped.sam"
  DISTR_LGT_PDF="$1_mapped2rRNA_read_length_distribution.pdf"
  if [ -s $DISTR_LGT_PDF ]
  then
    return 0
  else
    echo "Computing mapped to rRNA reads length distribution..."
    grep -v '^@' $INPUT_SAM_MAPPED_RNA | awk '{print length($10)}' | $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_READ_LENGTH_DISTRIBUTION $DISTR_LGT_PDF -e $INPUT_SAM_MAPPED_RNA
    if [ ! -s $DISTR_LGT_PDF ]
    then
      echo "Cannot compute mapped to rRNA reads length distribution !"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "PDF file : $DISTR_LGT_PDF generated"
    echo "End of mapping to rRNA reads length distribution."
  fi
}

#TODO, after ALL these calls to fastqc integrate all the reporst with MUltiQC! FIXME

# Run STAR to align reads to the reference genome
align_to_ref_genome() {
  for sample in ${SAMPLE_ARRAY[*]}
  do
    echo "Start mapping to reference genome :"
    WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
    if [ $WORKING_ANSWER_DEMULTIPLEXING = YES ]
    then
      DIR_ALIGN_STAR="${sample}_align_star/"
      INPUT_ALIGN_GENOME="${sample}_no_rRNA.fastq"
    else
      BASENAME=$(basename $sample .fastq)
      DIR_ALIGN_STAR="${BASENAME}_align_star/"
      INPUT_ALIGN_GENOME="${BASENAME}_no_rRNA.fastq"
    fi
    if [ -e $DIR_ALIGN_STAR ]
    then
      if [ -s "${DIR_ALIGN_STAR}Log.final.out" ] && [ -s "${DIR_ALIGN_STAR}Aligned.out.sam" ]
      then
        echo "Mapping already done for $sample."
        #return 0
      fi
    else
      mkdir -p $DIR_ALIGN_STAR
      if [ $? -ne 0 ]
      then
        echo "Cannot create the directory !"
        error_exit "$LINENO: An error has occurred."
      fi
      #FIXME check optimal STAR alignment parameters.
      docker run --rm -u $(id -u):$(id -g) -v $WORKDIR:/home -v $PATH_TO_GENOME_INDEX:/root -w /home genomicpariscentre/star bash -c "STAR --runThreadN $(nproc) --genomeDir /root --readFilesIn $INPUT_ALIGN_GENOME --outFileNamePrefix $DIR_ALIGN_STAR --outSAMunmapped Within --outFilterMismatchNmax $MAX_ALLOWED_MISMATCHES --quantMode TranscriptomeSAM --seedSearchStartLmax $SEED_SEARCH_POINT --outFilterScoreMinOverLread $FILTER_SCORE_MIN --outFilterMatchNminOverLread $FILTER_MATCH_MIN --winAnchorMultimapNmax $MAX_LOCI_ALLOWED --outFilterMultimapScoreRange $MULTIMAP_SCORE_RANGE"
      if [ ! -s "${DIR_ALIGN_STAR}Aligned.out.sam" ]
      then
        echo "STAR cannot run correctly ! Check your genome index path."
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Directory $DIR_ALIGN_STAR generated"
      echo "End of mapping to reference genome."
    fi
  done
}


# Performe the metagene analysis
#FIXME, we do not do this anymore, we use the longgest transcript.
metagene_analysis() {
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    echo "Metagene analysis :"
    BASENAME_METAGENE="metagene_orfs"
    LANDMARK="cds_start"
    INPUT_ANNOTATION=$(basename $PATH_TO_ANNOTATION_FILE)
    ANNOTATION_PREFIX=${INPUT_ANNOTATION:0:-4}
    ANNOTATION_CDSLANDMARK="${ANNOTATION_PREFIX}_cdslandmark.gtf"
    SORTED_ANNOTATION_CDSLANDMARK="${ANNOTATION_PREFIX}_cdslandmark_sorted.gtf"
    if [ -s metagene_orfs_rois.txt ] && [ -s metagene_orfs_rois.bed ]
    then
      echo "Metagene analysis already done."
      return 0
    else
      echo "Add CoDing Sequences landmark to annotations :"
      $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_ADD_START_GTF -i $PATH_TO_ANNOTATION_FILE -o $ANNOTATION_CDSLANDMARK
      if [ $? -ne 0 ]
      then
        echo "Cannot add CDS landmarks correctly."
        error_exit "$LINENO: An error has occurred."
      fi
      sort -T $TMPDIR -k1,1 -k4,4n $ANNOTATION_CDSLANDMARK > $SORTED_ANNOTATION_CDSLANDMARK
      if [ $? -ne 0 ]
      then
        echo "Annotations file cannot be sorted."
        error_exit "$LINENO: An error has occurred."
      fi
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home -v "${WORKDIR}/${SORTED_ANNOTATION_CDSLANDMARK}:/run/${SORTED_ANNOTATION_CDSLANDMARK}" parisepigenetics/plastid:0.4.6 bash -c "metagene generate $BASENAME_METAGENE --landmark $LANDMARK --sorted --annotation_files /run/${SORTED_ANNOTATION_CDSLANDMARK}"
      if [ $? -ne 0 ]
      then
        echo "Metagene analysis cannot run correctly."
        error_exit "$LINENO: An error has occurred."
      fi
      gzip -9 $SORTED_ANNOTATION_CDSLANDMARK
      if [ $? -ne 0 ]
      then
        echo "Cannot compress ${SORTED_ANNOTATION_CDSLANDMARK} correctly."
        error_exit "$LINENO: An error has occurred."
      fi
      echo "End of metagene analysis."
#      rm $ANNOTATION_CDS_EXON_START
      rm $ANNOTATION_CDSLANDMARK
    fi
  else
    return 0
  fi
}


# Align to reference genome for the seqcluster reads. (TODO check it... validate it)
align_to_ref_genome_seqcluster() {
  echo "Continuation of multi-mapped reads analysis (by Seqcluster) :"
  echo "Starting of mapping to reference genome for multi-mapped reads analysis :"
  DIR_MAPPING_GENOME="multi_reads_analysis/mappingStep/"
  INPUT_MAPPING_GENOME="multi_reads_analysis/mappingStep/seqs_no_rRNA.fastq"
  if [ -s "multi_reads_analysis/mappingStep/Aligned.out.sam" ]
  then
    echo "Mapping to reference genome for multi-mapped reads analysis already done."
    return 0
  else
    docker run --rm -u $(id -u):$(id -g) -v $WORKDIR:/home -v $PATH_TO_GENOME_INDEX:/root -w /home genomicpariscentre/star bash -c "STAR --runThreadN $(nproc) --genomeDir /root --readFilesIn $INPUT_MAPPING_GENOME --outSAMunmapped Within --outFilterMismatchNmax $MAX_ALLOWED_MISMATCHES  --seedSearchStartLmax $SEED_SEARCH_POINT --outFilterScoreMinOverLread $FILTER_SCORE_MIN --outFilterMatchNminOverLread $FILTER_MATCH_MIN --winAnchorMultimapNmax $MAX_LOCI_ALLOWED --outFilterMultimapScoreRange $MULTIMAP_SCORE_RANGE --outFilterMultimapNmax $MULTIMAP_N_MAX --outFileNamePrefix $DIR_MAPPING_GENOME"
    if [ $? -ne 0 ]
    then
      echo "STAR cannot run correcly ! Check your genome index path."
      error_exit "$LINENO: An error has occurred."
    fi
    echo "Outputs generated in $DIR_MAPPING_GENOME"
    echo "End of mapping to reference genome for multi-mapped reads analysis."
  fi
}


# Filter the SAM file to get conserve uniq reads
samFiltering() {
#  WORKING_ANSWER_KEEP_MULTIREAD=${ANSWER_KEEP_MULTIREAD^^}
  SAM_INPUT="$1_align_star/Aligned.out.sam"
  FILTERED_SAM_UNIQUE_OUTPUT="$1_align_filtered.sam"
  FILTERED_SAM_MULTI_OUTPUT="$1_align_multi.sam"
  LOGFILE="$1_align_filtering.log"
  echo "Starting of SAM file filtering :"
  if [ -s $SAM_INPUT ]
  then
    if [ -s $FILTERED_SAM_UNIQUE_OUTPUT ]
    then
      echo "SAM file filering already done."
      return 0
    else
#        if [ $WORKING_ANSWER_KEEP_MULTIREAD = YES ]
#        then
      grep -v '^@' $SAM_INPUT | awk '$2 != 4 {print $0}' | sort -k 1,1 -T $TMPDIR | $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_SAM_FILTERING -i $SAM_INPUT -o $FILTERED_SAM_UNIQUE_OUTPUT -m $FILTERED_SAM_MULTI_OUTPUT > $LOGFILE
#        else
#          grep -v '^@' $SAM_INPUT | awk '$2 != 4 {print $0}' | sort -k 1,1 -T $TMPDIR | $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_SAM_FILTERING -i $SAM_INPUT -o $FILTERED_SAM_UNIQUE_OUTPUT > $LOGFILE
#        fi
      if [ $? -ne 0 ]
      then
        echo "SAM file filtering cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Log file : $LOGFILE generated."
      echo "End of SAM file filtering."
    fi
  else
    echo "You need a SAM file to launch this step !"
    error_exit "$LINENO: An error has occurred."
  fi
}


# Compute the reads length distribution after alignment to the reference genome and the SAM file filtering (uniquely mapped only)
mapped_to_genome_distrib_length() {
  SAM_FILTERED_INPUT="$1_align_filtered.sam"
  DISTR_LGT_PDF="$1_uniquely_mapped2genome_read_length_distribution.pdf"
  if [ -s $DISTR_LGT_PDF ]
  then
    echo "Computing uniquely mapped to genome reads length distribution already done."
    return 0
  else
    echo "Computing uniquely mapped to genome reads length distribution :"
    grep -v '^@' $SAM_FILTERED_INPUT | awk '{print length($10)}' | $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_READ_LENGTH_DISTRIBUTION $DISTR_LGT_PDF -e $SAM_FILTERED_INPUT
    if [ ! -s $DISTR_LGT_PDF ]
    then
      echo "Cannot compute uniquely mapped to genome reads length distribution !"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "PDF file : $DISTR_LGT_PDF generated."
    echo "End of computing uniquely mapped to genome reads length distribution."
  fi
}


# Compute the multi-reads length distribution after alignment to reference genome
# TODO This can be written as a generic function for all the length distribution analysis
multimapped_to_genome_distrib_length() {
  WORKING_ANSWER_KEEP_MULTIREAD=${ANSWER_KEEP_MULTIREAD^^}
  if [ $WORKING_ANSWER_KEEP_MULTIREAD = YES ]
  then
    SAM_MULTIREAD_INPUT="$1_align_multi.sam"
    DISTR_LGT_PDF="$1_multimapped_to_genome_read_length_distribution.pdf"
    if [ -s $DISTR_LGT_PDF ]
    then
      echo "Computing multi-mapped to genome reads length distribution already done."
      return 0
    else
      echo "Computing multi-mapped to genome reads length distribution :"
      grep -v '^@' $SAM_MULTIREAD_INPUT | awk '{print length($10)}' | $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_READ_LENGTH_DISTRIBUTION $DISTR_LGT_PDF -e $SAM_MULTIREAD_INPUT
      if [ ! -s $DISTR_LGT_PDF ]
      then
        echo "Cannot compute multi-mapped to genome reads length distribution !"
        error_exit "$LINENO: An error has occurred."
      fi
    fi
    echo "PDF file : $DISTR_LGT_PDF generated."
    echo "End of computing multi-mapped to genome reads length distribution."
  else
    return 0
  fi
}


# Convert the filtered SAM file into a sorted BAM file which is indexed in BAI
sam_to_bam() {
  FILTERED_SORTED_ALIGNMENT="$1_align_filtered.sorted"  # Sorted alignment basename for BAM & BAI files
  FILTERED_SAM="$1_align_filtered.sam"
  if [ -s $FILTERED_SAM ]
  then
    if [ -s "${FILTERED_SORTED_ALIGNMENT}.bam" ]
    then
      return 0
    else
      echo "Starting of Samtools"
      # SAM to BAM conversion + sorting of BAM file
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/samtools bash -c "samtools view -Sb $FILTERED_SAM | samtools sort - $FILTERED_SORTED_ALIGNMENT"
      # Index BAI
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/samtools bash -c "samtools index "${FILTERED_SORTED_ALIGNMENT}.bam" "${FILTERED_SORTED_ALIGNMENT}.bai""
      if [ ! -s "${FILTERED_SORTED_ALIGNMENT}.bam" ]
      then
        echo "Samtools cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Sorted-indexed alignment : '${FILTERED_SORTED_ALIGNMENT}.bam' generated. You can use it in a genome browser (e.g IGV)"
      echo "End of Samtools."
    fi
  else
    echo "You need a filtered SAM file to launch this step !"
    error_exit "$LINENO: An error has occurred."
  fi
}


p_offset_analysis() {
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    echo "P-site offsets analysis :"
    INPUT_METAGENE_ORF="metagene_orfs_rois.txt"
    #BASENAME_RIBOPROFILE="$1_riboprofile"
#    INPUT_COUNT_FILES=()
    #METAGENE_PROFILE="${BASENAME_RIBOPROFILE}_metagene_profiles.txt"
    #P_OFFSETS="${BASENAME_RIBOPROFILE}_p_offsets.txt"
    #P_OFFSETS_CORRECTED="${BASENAME_RIBOPROFILE}_p_offsets_corrected.txt"
    WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
    if [ -s $INPUT_METAGENE_ORF ]
    then
      for cond in $CONDITION_ARRAY_UNIQ ## TODO Check if we can do it wuth GNU parallel. Problem is the passing of the arrays SAMPLE & CONDITIONS
      do
        echo  "P-site offsets analysis for $cond condition :"
        BASENAME_RIBOPROFILE="${cond}_riboprofile"
        METAGENE_PROFILE="${BASENAME_RIBOPROFILE}_metagene_profiles.txt"
        P_OFFSETS="${BASENAME_RIBOPROFILE}_p_offsets.txt"
        WORKING_P_OFFSETS="${cond}_p_offsets.txt"
        P_OFFSETS_CORRECTED="${cond}_p_offsets_corrected.txt"
        INPUT_COUNT_FILES=()
        if [ -s $WORKING_P_OFFSETS ] || [ -s $P_OFFSETS_CORRECTED ]
        then
          echo "P-site offsets analysis already done."
          return 0
        else
#            CONDITION_ARRAY+=($2)
#            SAMPLE_ARRAY+=($3)
#            echo ${CONDITION_ARRAY[*]}
#            echo ${SAMPLE_ARRAY[*]}
          for index in ${!CONDITION_ARRAY[*]}
          do
            CURRENT_CONDITION=${CONDITION_ARRAY[$index]}
#             echo $CURRENT_CONDITION
            if [ $CURRENT_CONDITION = $cond ]
            then
              if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
              then
                SAMPLE=${SAMPLE_ARRAY[$index]}
              else
                SAMPLE=$(basename ${SAMPLE_ARRAY[$index]} .fastq)
              fi
              INPUT_COUNT_FILES+=("${SAMPLE}_align_filtered.sorted.bam")
#                echo $SAMPLE
            fi
          done
          #echo ${INPUT_COUNT_FILES[*]}
          docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home parisepigenetics/plastid:0.4.6 bash -c "psite $INPUT_METAGENE_ORF $BASENAME_RIBOPROFILE --min_length $MIN_LENGTH_POFFSET --max_length $MAX_LENGTH_POFFSET --require_upstream --aggregate --count_files ${INPUT_COUNT_FILES[*]} --keep --default $DEFAULT_POFFSET"
          if [ $? -ne 0 ]
          then
            echo "P-site offsets analysis cannot run correctly."
            error_exit "$LINENO: An error has occurred."
          fi
          grep -v "##" $P_OFFSETS | grep -v "length" | grep -v "default" > $WORKING_P_OFFSETS
          WORKING_ANSWER_PSITE_CORRECTION=${ANSWER_PSITE_CORRECTION^^}
          if [ $WORKING_ANSWER_PSITE_CORRECTION = YES ]
          then
            $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_PSITE_AUTOCORRECTION -i $WORKING_P_OFFSETS -m $METAGENE_PROFILE -o $P_OFFSETS_CORRECTED
            if [ $? -ne 0 ]
            then
              echo "P-site offsets auto-correction cannot run correctly."
              error_exit "$LINENO: An error has occurred."
            fi
          fi
          echo "End of P-site offsets analysis."
          echo "If you see a P-site offset of 0 or some odd p-offsets in $WORKING_P_OFFSETS and $P_OFFSETS_CORRECTED, you should correct them manually."
        fi
      done
      WORKING_STOP_EXEC_PSITE_CORRECTION=${STOP_EXEC_PSITE_CORRECTION^^}
      if [ $WORKING_STOP_EXEC_PSITE_CORRECTION = YES ]
      then
        echo "Execution stoped by STOP_EXEC_PSITE_CORRECTION parameter."
        error_exit "$LINENO: An error has occurred."
      fi
    else
      echo "You need $INPUT_METAGENE_ORF file to launch P-site offsets analysis."
      error_exit "$LINENO: An error has occurred."
    fi
  else
    return 0
  fi
}


rna_seq_quantification() {
# The user MUST do transcriptome mapping of their RNA-seq data before (TODO, why it needs to be done like this..)
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    echo "mRNA quantification :"
    QUANTIFICATION_OUTPUT_DIR="$1_mRNA_quantification"
    TRANSCRIPTOME_FASTA_FILE=$(basename $PATH_TO_REFERENCE_TRANSCRIPTOME_FILE)
    RNA_BAM="$1_mRNA.transcriptome.mapping.bam"
    #TRANSCRIPTOME_FASTA_DIR=$(dirname $PATH_TO_REFERENCE_TRANSCRIPTOME_FILE)
    if [ -s "${QUANTIFICATION_OUTPUT_DIR}/quant.sf" ]
    then
      echo "mRNA quantification already done for $1."
      return 0
    else
      if [ -s $RNA_BAM ]
      then
        mkdir -p $QUANTIFICATION_OUTPUT_DIR
        if [ $? -ne 0 ]
        then
          echo "Cannot make $QUANTIFICATION_OUTPUT_DIR."
          error_exit "$LINENO: An error has occurred."
        fi
        docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home -v "${PATH_TO_REFERENCE_TRANSCRIPTOME_FILE}:/run/${TRANSCRIPTOME_FASTA_FILE}" genomicpariscentre/ribomap bash -c "salmon quant -t /run/${TRANSCRIPTOME_FASTA_FILE} -l $2 -a $RNA_BAM -o $QUANTIFICATION_OUTPUT_DIR --biasCorrect 2> /dev/null"
        if [ $? -ne 0 ]
        then
          echo "mRNA quantification cannot run correctly."
          error_exit "$LINENO: An error has occurred."
        fi
        echo "End of mRNA quantification."
      else
        echo "$RNA_BAM file doesn't exist."
        error_exit "$LINENO: An error has occurred."
      fi
    fi
  else
    return 0
  fi
}


cds_range_building() {
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    echo "CDS range building :"
    ANNOTATIONS_FILE=$(basename $PATH_TO_ANNOTATION_FILE)
    ANNOTATION_PREFIX="${ANNOTATIONS_FILE:0:-4}"
    CDS_RANGE="${ANNOTATION_PREFIX}_cds_range.txt"
    if [ -s $CDS_RANGE ]
    then
      echo "CDS range building already done."
      return 0
    else
      $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_CDS_RANGE_GENERATOR -i $PATH_TO_ANNOTATION_FILE -o $CDS_RANGE
      if [ $? -ne 0 ]
      then
        echo "Cannot build CDS range in $CDS_RANGE file correctly."
        error_exit "$LINENO: An error has occurred."
      fi
      echo "End of CDS range buiding"
    fi
  else
    return 0
  fi
}


isoform_level_estimation() {
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    echo "Starting of Isoform-level ribosome occupancy guided by transcript abundance :"
    WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
#      RIBOMAP_OUTPUT="$1_isoform_level_estimation/$1_"
    ANNOTATIONS_FILE=$(basename $PATH_TO_ANNOTATION_FILE)
    ANNOTATION_PREFIX="${ANNOTATIONS_FILE:0:-4}"
    CDS_RANGE="${ANNOTATION_PREFIX}_cds_range.txt"
    TRANSCRIPTOME_FASTA_FILE=$(basename $PATH_TO_REFERENCE_TRANSCRIPTOME_FILE)
#      RNASEQ_BAM="$1_mRNA.transcriptome.mapping.bam"
#      RIBO_BAM="$1_align_star/Aligned.toTranscriptome.out.bam"
    RIBOMAP_OUTDIR=isoform_level_estimation_analyis
    mkdir -p $RIBOMAP_OUTDIR
    if [ $? -ne 0 ]
    then
      echo "Cannot build $RIBOMAP_OUTDIR directory."
      error_exit "$LINENO: An error has occurred."
    fi
    for sample in ${SAMPLE_ARRAY[*]}
    do
      if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
      then
        SAMPLE=$sample
      else
        SAMPLE=$(basename $sample .fastq)
      fi
      RIBOMAP_OUTPUT="${SAMPLE}_isoform_level_estimation"
      RIBOMAP_OUTPUT_DIR_SAMPLE="${RIBOMAP_OUTDIR}/${SAMPLE}"
      RNASEQ_BAM="${SAMPLE}_mRNA.transcriptome.mapping.bam"
      RIBO_BAM="${SAMPLE}_align_star/Aligned.toTranscriptome.out.bam"
      RNA_QUANTIFICATION="${SAMPLE}_mRNA_quantification/quant.sf"
      LOGFILE="${SAMPLE}_isoform_level_assignment.log"
      if [ -e $RIBOMAP_OUTPUT_DIR_SAMPLE ]
      then
        if [ -s "${RIBOMAP_OUTPUT_DIR_SAMPLE}/${RIBOMAP_OUTPUT}_abundant.list" ]
        then
          echo "Isoform-level estimation already done for $SAMPLE"
        fi
      else
        for index in ${!SAMPLE_ARRAY[*]}
        do
          if [ $SAMPLE = ${SAMPLE_ARRAY[$index]} ]
          then
            SAMPLE_COND=${CONDITION_ARRAY[$index]}
          fi
        done
#        echo $SAMPLE
#        echo $SAMPLE_COND
        if [ -s "${SAMPLE_COND}_p_offsets_corrected.txt" ]
        then
          POFFSET="${SAMPLE_COND}_p_offsets_corrected.txt"
        elif [ -s "${SAMPLE_COND}_p_offsets.txt" ]
        then
          POFFSET="${SAMPLE_COND}_p_offsets.txt"
        else
          echo "P-site offsets file is needed for this step."
          error_exit "$LINENO: An error has occurred."
        fi
        if [ ! $RNASEQ_LIBTYPE = "unstranded" ]
        then
          docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -v "$PATH_TO_REFERENCE_TRANSCRIPTOME_FILE:/run/${TRANSCRIPTOME_FASTA_FILE}" -w /home genomicpariscentre/ribomap bash -c "riboprof --fasta /run/${TRANSCRIPTOME_FASTA_FILE} --mrnabam $RNASEQ_BAM --ribobam $RIBO_BAM --min_fplen $MIN_FPLEN --max_fplen $MAX_FPLEN --offset $POFFSET --cds_range $CDS_RANGE --sf $RNA_QUANTIFICATION --out $RIBOMAP_OUTPUT --useSecondary --tabd_cutoff 0 > $LOGFILE"
        else
          docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -v "$PATH_TO_REFERENCE_TRANSCRIPTOME_FILE:/run/${TRANSCRIPTOME_FASTA_FILE}" -w /home genomicpariscentre/ribomap bash -c "riboprof --fasta /run/${TRANSCRIPTOME_FASTA_FILE} --mrnabam $RNASEQ_BAM --ribobam $RIBO_BAM --min_fplen $MIN_FPLEN --max_fplen $MAX_FPLEN --offset $POFFSET --cds_range $CDS_RANGE --sf $RNA_QUANTIFICATION --out $RIBOMAP_OUTPUT --useSecondary --useRC --tabd_cutoff 0 > $LOGFILE"
        fi
        if [ $? -ne 0 ]
        then
          echo "Isoform-level estimation cannot run corrctly."
          error_exit "$LINENO: An error has occurred."
        fi
        echo "End of Isoform-level estimation for ${SAMPLE}."
        mkdir -p $RIBOMAP_OUTPUT_DIR_SAMPLE
        if [ $? -ne 0 ]
        then
          echo "Cannot build $RIBOMAP_OUTPUT_DIR_SAMPLE directory."
          error_exit "$LINENO: An error has occurred."
        fi
        for file in $(ls ${RIBOMAP_OUTPUT}*)
        do
          mv $file $RIBOMAP_OUTPUT_DIR_SAMPLE
        done
        if [ $? -ne 0 ]
        then
          echo "Cannot move $SAMPLE results in $RIBOMAP_OUTPUT_DIR_SAMPLE directory"
          error_exit "$LINENO: An error has occurred."
        fi
      fi
    done
  else
    return 0
  fi
}


sam_to_bam_seqcluster() {
#TODO Generic function of this similar to sam_to_bam
  echo "Continuation of multi-mapped reads analysis (by Seqcluster) :"
  echo "Conversation of genome mapping file for multi-mapped reads analysis from SAM into BAM :"
  SORTED_ALIGNMENT="multi_reads_analysis/mappingStep/Aligned.out.sorted"
  INPUT_ALIGNMENT="multi_reads_analysis/mappingStep/Aligned.out.sam"
  if [ -s $INPUT_ALIGNMENT ]
  then
    if [ -s "${SORTED_ALIGNMENT}.bam" ]
    then
      return 0
    else
      echo "Starting of Samtools"
      # SAM to BAM conversion + sorting of BAM file
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/samtools bash -c "samtools view -Sb $INPUT_ALIGNMENT | samtools sort - $SORTED_ALIGNMENT"
      # Index BAI
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/samtools bash -c "samtools index "${SORTED_ALIGNMENT}.bam" "${SORTED_ALIGNMENT}.bai""
      if [ ! -s "${SORTED_ALIGNMENT}.bam" ]
      then
        echo "Samtools cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "Sorted-indexed alignment : '${SORTED_ALIGNMENT}.bam' generated."
      echo "End of Samtools."
    fi
  else
    echo "You need $INPUT_ALIGNMENT to launch this step !"
    error_exit "$LINENO: An error has occurred."
  fi
}


clustering_step_seqcluster() {
  echo "Continuation of multi-mapped reads analysis (by Seqcluster) :"
  echo "Clustering step :"
  # Inputs files
  INPUT_ALIGNMENT_CLUSTERING="multi_reads_analysis/mappingStep/Aligned.out.sorted.bam"
  INPUT_SEQCLUSTER_MATRIX="multi_reads_analysis/prepareSamplesStep/seqs_no_rRNA.ma"
  INPUT_COUNTS_CLUSTER="multi_reads_analysis/clusteringStep/counts.tsv"
  # Outputs
  DIR_OUTPUT_CLUSTERING="multi_reads_analysis/clusteringStep"
  DIR_OUTPUT_ANADIFFCLUSTER="multi_reads_analysis/differential_analysis"
  # Input data (genome + annotations)
  DIRNAME_GENOME=$(dirname $PATH_TO_REFERENCE_GENOME_FILE)
  BASENAME_GENOME=$(basename $PATH_TO_REFERENCE_GENOME_FILE)
  BASENAME_ANNOTATION=$(basename $PATH_TO_ANNOTATION_FILE)
  if [ -s $INPUT_ALIGNMENT_CLUSTERING ]
  then
    if [ -e $DIR_OUTPUT_CLUSTERING ]
    then
      NB_FILE_DIR_OUTPUT_CLUSTERING=$(ls -R $DIR_OUTPUT_CLUSTERING | wc -l)
      let NB_FILE_DIR_OUTPUT_CLUSTERING=$NB_FILE_DIR_OUTPUT_CLUSTERING-1
      if [ $NB_FILE_DIR_OUTPUT_CLUSTERING -gt 1 ]
      then
        echo "Clustering step for multi-mapped analysis already done."
        return 0
      fi
    else
      # We generate FASTA genome index (.fai)
      if [ ! -s "${DIRNAME_GENOME}/${BASENAME_GENOME}.fai" ]
      then
        echo "Bulding of FASTA index (.fai) :"
        docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $DIRNAME_GENOME:/home -w /home genomicpariscentre/samtools bash -c "samtools faidx $BASENAME_GENOME"
        if [ $? -ne 0 ]
        then
          echo "Samtools cannot run correctly (for generate FASTA index '${BASENAME_GENOME}.fai')"
          error_exit "$LINENO: An error has occurred."
        fi
        echo "${BASENAME_GENOME%.*}.fai was generated in $DIRNAME_GENOME"
      fi
      # We run seqcluster clustering
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v /etc/passwd:/etc/passwd -v "${PATH_TO_REFERENCE_GENOME_FILE}:/run/${BASENAME_GENOME}" -v "${DIRNAME_GENOME}/${BASENAME_GENOME}.fai:/run/${BASENAME_GENOME}.fai" -v "${PATH_TO_ANNOTATION_FILE}:/run/${BASENAME_ANNOTATION}" -v ${WORKDIR}:/home -w /home genomicpariscentre/bcbio-nextgen bash -c "/usr/local/bin/seqcluster cluster -a $INPUT_ALIGNMENT_CLUSTERING -m $INPUT_SEQCLUSTER_MATRIX --feature_id $IDATTR -g /run/${BASENAME_ANNOTATION} -o $DIR_OUTPUT_CLUSTERING -ref /run/${BASENAME_GENOME} --db multi_mapped_analysis_database"
      if [ $? -ne 0 ]
      then
        echo "Clustering step cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      # We generate clusters countings
      echo "Building clusters counting files in $DIR_OUTPUT_CLUSTERING :"
      $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_COUNT_CLUSTER -i $INPUT_COUNTS_CLUSTER
      # We generate files to do differential analysis of clusters
      echo "Building files to do differential analysis of clusters in $DIR_OUTPUT_ANADIFFCLUSTER"
      mkdir -p $DIR_OUTPUT_ANADIFFCLUSTER
      cp *_Clustercount.tsv $DIR_OUTPUT_ANADIFFCLUSTER
      echo -e "label\\tfiles\\tgroup" >> "${DIR_OUTPUT_ANADIFFCLUSTER}/target.txt"
      for index in ${!SAMPLE_ARRAY[*]}
      do
        WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
        if [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
        then
          SAMPLE=${SAMPLE_ARRAY[$index]}
        else
          SAMPLE=$(basename ${SAMPLE_ARRAY[$index]} .fastq)
        fi
        echo -e "$SAMPLE\\t${SAMPLE}_Clustercount.tsv\\t${CONDITION_ARRAY[$index]}" >> "${DIR_OUTPUT_ANADIFFCLUSTER}/target.txt"
      done
      if [ $? -ne 0 ]
      then
        echo "Cannot building clusters counting files correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "End of building clusters counting files."
      echo "$DIR_OUTPUT_CLUSTERING was generated."
      echo "End of clustering step."
    fi
  else
    echo "$INPUT_ALIGNMENT_CLUSTERING doesn't exist !"
    error_exit "$LINENO: An error has occurred."
  fi
}


report_step_seqcluster() {
  echo "Continuation of multi-mapped reads analysis (by Seqcluster) :"
  echo "Report step :"
  INPUT_JSON="multi_reads_analysis/clusteringStep/seqcluster.json"
  DIROUT_REPORT="multi_reads_analysis/reportStep"
  BASENAME_GENOME=$(basename $PATH_TO_REFERENCE_GENOME_FILE)
  DIRNAME_GENOME=$(dirname $PATH_TO_REFERENCE_GENOME_FILE)
  if [ -s $INPUT_JSON ]
  then
    if [ -e $DIROUT_REPORT ]
    then
      NB_FILE_DIROUT_REPORT=$(ls -R $DIROUT_REPORT | wc -l)
      let NB_FILE_DIROUT_REPORT=$NB_FILE_DIROUT_REPORT-1
      if [ $NB_FILE_DIROUT_REPORT -gt 1 ]
      then
        echo "Report step already done."
        return 0
      fi
    else
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v /etc/passwd:/etc/passwd -v "${PATH_TO_REFERENCE_GENOME_FILE}:/run/${BASENAME_GENOME}" -v "${DIRNAME_GENOME}/${BASENAME_GENOME}.fai:/run/${BASENAME_GENOME}.fai" -v ${WORKDIR}:/home -w /home genomicpariscentre/bcbio-nextgen bash -c "/usr/local/bin/seqcluster report -j $INPUT_JSON -o $DIROUT_REPORT -r /run/${BASENAME_GENOME}"
      if [ $? -ne 0 ]
      then
        echo "Report step cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
    fi
  else
    echo "$INPUT_JSON doesn't exist !"
    error_exit "$LINENO: An error has occurred."
  fi
}


# Get longest transcript of each gene for CDS annotations from Ensembl 75 GTF
get_longest_transcripts_from_annotations() {
  INPUT_ANNOTATION=$(basename $PATH_TO_ANNOTATION_FILE)
  DIRNAME_ANNOTATIONS=$(dirname $PATH_TO_ANNOTATION_FILE)
  ANNOTATION_PREFIX=${INPUT_ANNOTATION:0:-4}
  CDS_ANNOTATIONS="${ANNOTATION_PREFIX}_only_cds.gtf"
  LONGEST_TRANSCRIPTS="${ANNOTATION_PREFIX}_longest_transcripts.txt"
  CDS_LONGEST_TRANSCRIPTS_LIST="${ANNOTATION_PREFIX}_only_cds_longest_transcripts.txt"
  CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS="${ANNOTATION_PREFIX}_only_cds_longest_transcripts.gtf"
  if [ ! -s $CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS ]
  then
    echo "Building annotations containing CoDing Sequences from longest transcripts :"
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -v $DIRNAME_ANNOTATIONS:/root -w /home genomicpariscentre/gff3-ptools bash -c "gtf-filter --keep-comments -o $CDS_ANNOTATIONS \"field feature == CDS\" /root/$INPUT_ANNOTATION"
    $PYTHON_SCRIPTS_PATH$PYTHON_SCRIPT_LONGEST_TRANSCRIPT -i $PATH_TO_ANNOTATION_FILE -o $CDS_LONGEST_TRANSCRIPTS_LIST
    grep -Ff $CDS_LONGEST_TRANSCRIPTS_LIST $CDS_ANNOTATIONS > $CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS
    echo "GTF annotations $CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS generated."
    echo "End of building annotations."
  else
    return 0
  fi
}


# We compute the number of reads in CDS (HTSeq-count)
htseq_count() {
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  FILTERED_SORTED_BAM="$1_align_filtered.sorted.bam"
  HTSEQCOUNT_FILE="$1_htseq.txt"
  HTSEQCOUNT_FILE_ANADIF_BABEL="$1_RPcounts.txt"
  ANNOTATIONS_FILE=$(basename $PATH_TO_ANNOTATION_FILE)
  ANNOTATION_PREFIX=${ANNOTATIONS_FILE:0:-4}
  CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS="${ANNOTATION_PREFIX}_only_cds_longest_transcripts.gtf"
  if [ -s $FILTERED_SORTED_BAM ] && [ -s $CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS ]
  then
    if [ -s $HTSEQCOUNT_FILE ] || [ -s DifferentialAnalysis/$HTSEQCOUNT_FILE ]
    then
      return 0
    else
      echo "Starting of expression estimation (counted reads/gene) :"
      docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR:/home -w /home genomicpariscentre/htseq bash -c "htseq-count --mode $MODE_FOR_MULTIPLE_FEATURES_READS --type $FEATURE_TYPE --idattr $IDATTR --stranded $STRANDED --format $FILETYPE $FILTERED_SORTED_BAM $CDS_LONGEST_TRANSCRIPTS_ANNOTATIONS > $HTSEQCOUNT_FILE"
      if [ $? -ne 0 ]
      then
        echo "HTSeq-Count cannot run correctly ! Check your --stranded option on http://www-huber.embl.de/users/anders/HTSeq/doc/count.html"
        error_exit "$LINENO: An error has occurred."
      fi
      mkdir -p DifferentialAnalysis
      chown $(id -u):$(id -g) -R DifferentialAnalysis
      # We copy counting files
      if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
      then
        grep -v __.* $HTSEQCOUNT_FILE > $HTSEQCOUNT_FILE_ANADIF_BABEL
        cp $HTSEQCOUNT_FILE_ANADIF_BABEL DifferentialAnalysis
        cp $1_mRNAcounts.txt DifferentialAnalysis
      else
        cp $HTSEQCOUNT_FILE DifferentialAnalysis
      fi
      # We generate (later) & copy target.txt --> we can't do this here because we are in GNU parallel --> done in anadiffSARTools fonction
      #echo "Building of the design file target.txt to do differential analysis :"
      #if [ -s target.txt ]
      #then
      #  cp target.txt DifferentialAnalysis
      #fi
      if [ $? -ne 0 ]
      then
        echo "HTSeq-Count cannot run correctly !"
        error_exit "$LINENO: An error has occurred."
      fi
      echo "DifferentialAnalysis generated. It contains :"
      ls DifferentialAnalysis
      echo "End of expression estimation."
    fi
  else
    echo "You need a filtered-sorted BAM file to launch this step !"
    error_exit "$LINENO: An error has occurred."
  fi
}


# Run Babel if the user has provided RNA-seq data: here, build of counts matrix
build_rnaseq_ribopro_counting_tables() {
  if [ -e DifferentialAnalysis ]
  then
    WORKDIR_ANADIFF=$(readlink -f DifferentialAnalysis)
  fi
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  # If user has RNA-seq countings, we build RNAseq and Ribosome Profiling counting tables
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    echo "Building matrix expression for Babel :"
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR_ANADIFF:/home -v $R_SCRIPTS_PATH:/root -w /home genomicpariscentre/babel Rscript "/root/${R_SCRIPT_BUILD_COUNTING_TABLE_RNASEQ}" ${SAMPLES[@]}
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR_ANADIFF:/home -v $R_SCRIPTS_PATH:/root -w /home genomicpariscentre/babel Rscript "/root/${R_SCRIPT_BUILD_COUNTING_TABLE_RP}" ${SAMPLES[@]}
    echo "End of building matrix expression."
  else
    return 0
  fi
}


# Run Babel if the user has provided RNA-seq data: here differntial analysis and its permutation test
#TODO revisit that one of the most important parts.
anadif_babel() {
  WORKDIR_ANADIFF=$(readlink -f DifferentialAnalysis)
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  # If user has RNA-seq counting, we use Babel R package
  if [ $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    #if [ -z $CONDITION_ARRAY ]
    #then
    #  echo "Give your (biological) condition array."
    #  error_exit "$LINENO: An error has occurred."
    #fi
    echo "Differential analysis :"
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR_ANADIFF:/home -w /home -v $R_SCRIPTS_PATH:/root genomicpariscentre/babel Rscript "/root/${R_SCRIPT_ANADIFF_BABEL}" ${CONDITION_ARRAY[@]}
    echo "Permutation test :"
    docker run --rm -u $(id -u):$(id -g) -v $TMPDIR:/tmp -v $WORKDIR_ANADIFF:/home -w /home -v $R_SCRIPTS_PATH:/root genomicpariscentre/babel Rscript "/root/${R_SCRIPT_PERMT_TEST_BABEL}" ${CONDITION_ARRAY[@]}
    echo "End of statistical analysis."
  else
    echo "Differential analysis is already done."
    return 0
  fi
}


anadif_sartools() {
#TODO also revisit that!
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  WORKING_DIFFERENTIAL_ANALYSIS_PACKAGE=${DIFFERENTIAL_ANALYSIS_PACKAGE^^}
  # If user hasn't RNA-seq counting, we use SARTools R package
  if [ ! $WORKING_ANSWER_RNASEQ_DATA = YES ]
  then
    if [ -e DifferentialAnalysis ]
    then
      WORKDIR_ANADIFF=$(readlink -f DifferentialAnalysis)
    else
      echo "DifferentialAnalysis doesn't exist !"
      error_exit "$LINENO: An error has occurred."
    fi
    echo "Building design file target.txt to do differential analysis :"
    echo -e "label\\tfiles\group" >> ${WORKDIR_ANADIFF}/target.txt
    for index in ${!SAMPLE_ARRAY[*]}
    do
      WORKING_ANSWER_DEMULTIPLEXING=${ANSWER_DEMULTIPLEXING^^}
      if  [ $WORKING_ANSWER_DEMULTIPLEXING = "YES" ]
      then
        SAMPLE=${SAMPLE_ARRAY[$index]}
      else
        SAMPLE=$(basename ${SAMPLE_ARRAY[$index]} .fastq)
      fi
      echo -e "${SAMPLE}_htseq.txt\\t${CONDITION_ARRAY[$index]}" >> ${WORKDIR_ANADIFF}/target.txt
    done
    PARAM=(/home $1 $2 target.txt /home $3) # 1 Working directory, 2 project name, 3 author name, 4 design (target.txt), 5 directory with counts files, 6 biological condition of reference
    WORK_PARAM=$(echo ${PARAM[*]})
    PARAMETERS=$(echo $WORK_PARAM)
    # EdgeR is launch by default if not specified (because Babel uses edgeR)
    if [ $WORKING_DIFFERENTIAL_ANALYSIS_PACKAGE = DESEQ2 ]
    then
      echo "Diffrential analysis :"
      docker run --rm -u $(id -u):$(id -g) -v $WORKDIR_ANADIFF:/home -w /home -v $R_SCRIPTS_PATH:/root genomicpariscentre/sartools Rscript "/root/${R_SCRIPT_ANADIFF_SARTOOLS_DESEQ2}" $PARAMETERS
    else
      echo "Diffrential analysis :"
      docker run --rm -u $(id -u):$(id -g) -v $WORKDIR_ANADIFF:/home -w /home -v $R_SCRIPTS_PATH:/root genomicpariscentre/sartools Rscript "/root/${R_SCRIPT_ANADIFF_SARTOOLS_EDGER}" $PARAMETERS
    fi
    echo "End of differential analysis."
  else
    return 0
  fi
}


# Launching of SARTools with the counts from SeqCluster
#TODO we do not weed the SeqCluster analysis anymore.
anadif_clusters() {
  WORKING_ANSWER_RNASEQ_DATA=${ANSWER_RNASEQ_DATA^^}
  WORKING_DIFFERENTIAL_ANALYSIS_PACKAGE=${DIFFERENTIAL_ANALYSIS_PACKAGE^^}
  DIR_ANADIFCLUSTER="multi_reads_analysis/differential_analysis"
  if [ -e $DIR_ANADIFCLUSTER ]
  then
    WORKDIR_ANADIFCLUSTER=$(readlink -f $DIR_ANADIFCLUSTER)
  else
    echo "$DIR_ANADIFCLUSTER doesn't exist !"
    error_exit "$LINENO: An error has occurred."
  fi
  PARAM=(/home $1 $2 target.txt /home $3)
  WORK_PARAM=$(echo ${PARAM[*]})
  PARAMETERS=$(echo $WORK_PARAM)
  echo "Diffrential analysis for cluster :"
  if [ $WORKING_DIFFERENTIAL_ANALYSIS_PACKAGE=DESEQ2 ]
  then
    docker run --rm -u $(id -u):$(id -g) -v $WORKDIR_ANADIFCLUSTER:/home -w /home -v $R_SCRIPTS_PATH:/root genomicpariscentre/sartools Rscript "/root/${R_SCRIPT_ANADIFF_SARTOOLS_DESEQ2}" $PARAMETERS
  else
    docker run --rm -u $(id -u):$(id -g) -v $WORKDIR_ANADIFCLUSTER:/home -w /home -v $R_SCRIPTS_PATH:/root genomicpariscentre/sartools Rscript "/root/${R_SCRIPT_ANADIFF_SARTOOLS_EDGER}" $PARAMETERS
  fi
  echo "End of differential analysis for clusters."
}

# Export functions so that they will be available to parallel.
export -f fastqc_quality_control
export -f demultiplexing
export -f raw_quality_report
export -f removeBadIQF
export -f removeBadIQF_report
export -f removePCRduplicates
export -f index_Adapter_trimming
export -f index_Adapter_trimming_report
export -f threePrime_trimming
export -f threePrime_trimming_report
export -f size_Selection
export -f size_Selection_report
export -f collapse_step_seqcluster
export -f align_To_R_RNA
export -f unmmaped_to_rRNA_report
export -f align_To_R_RNA_seqcluster
export -f mapped_to_R_RNA_distrib_length
export -f align_to_ref_genome
export -f metagene_analysis
export -f align_to_ref_genome_seqcluster
export -f samFiltering
export -f mapped_to_genome_distrib_length
export -f multimapped_to_genome_distrib_length
export -f sam_to_bam
export -f p_offset_analysis
export -f rna_seq_quantification
export -f cds_range_building
export -f isoform_level_estimation
export -f sam_to_bam_seqcluster
export -f clustering_step_seqcluster
export -f report_step_seqcluster
export -f get_longest_transcripts_from_annotations
export -f htseq_count
export -f build_rnaseq_ribopro_counting_tables
export -f anadif_babel
export -f anadif_sartools
export -f anadif_clusters


### MAIN ###
parallel --no-notice --xapply demultiplexing ::: $WORKING_SAMPLE_ARRAY ::: $WORKING_SAMPLE_INDEX_ARRAY
if [ $? -ne 0 ]
then
  # echo $?
  error_exit "$LINENO: demultiplexing An error has occurred." ## TODO Check if we can get the exit status of the last job
fi
wait

parallel --no-notice raw_quality_report {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: raw_quality_report An error has occurred."
fi
wait

parallel --no-notice removeBadIQF {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: removeBadIQF An error has occurred."
fi
wait

parallel --no-notice removeBadIQF_report {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: removeBadIQF_report An error has occurred."
fi
wait

parallel --no-notice removePCRduplicates {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: removePCRduplicates An error has occurred."
fi
wait

parallel --no-notice --xapply index_Adapter_trimming {.} ::: $WORKING_SAMPLE_ARRAY ::: $WORKING_SAMPLE_INDEX_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: index_Adapter_trimming An error has occurred."
fi
wait

parallel --no-notice index_Adapter_trimming_report {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: index_Adapter_trimming_report An error has occurred."
fi
wait

parallel --no-notice --xapply threePrime_trimming {.} ::: $WORKING_SAMPLE_ARRAY ::: $ADAPTER_SEQUENCE_THREE_PRIME
if [ $? -ne 0 ]
then
  error_exit "$LINENO: threePrime_trimming An error has occurred."
fi
wait

parallel --no-notice threePrime_trimming_report {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: threePrime_trimming_report An error has occurred."
fi
wait

parallel --no-notice size_Selection {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: size_Selection An error has occurred."
fi
wait

parallel --no-notice size_Selection_report {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: size_Selection_report An error has occurred."
fi
wait

parallel --no-notice collapse_step_seqcluster {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: collapse_step_seqcluster An error has occurred."
fi
wait

prepareSamples_step_seqcluster
if [ $? -ne 0 ]
then
  error_exit "$LINENO: prepareSamples_step_seqcluster An error has occurred."
fi
wait

align_To_R_RNA
if [ $? -ne 0 ]
then
  error_exit "$LINENO: align_To_R_RNA An error has occurred."
fi
wait

parallel --no-notice unmmaped_to_rRNA_report {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: align_To_R_RNA An error has occurred."
fi
wait

align_To_R_RNA_seqcluster
if [ $? -ne 0 ]
then
  error_exit "$LINENO: align_To_R_RNA_seqcluster An error has occurred."
fi
wait

parallel --no-notice mapped_to_R_RNA_distrib_length {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: mapped_to_R_RNA_distrib_length An error has occurred."
fi
wait

align_to_ref_genome
if [ $? -ne 0 ]
then
  error_exit "$LINENO: align_to_ref_genome An error has occurred."
fi
wait

metagene_analysis
if [ $? -ne 0 ]
then
  error_exit "$LINENO: metagene_analysis An error has occurred."
fi
wait

align_to_ref_genome_seqcluster
if [ $? -ne 0 ]
then
  error_exit "$LINENO: align_to_ref_genome_seqcluster An error has occurred."
fi
wait

parallel --no-notice samFiltering {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: samFiltering An error has occurred."
fi
wait

parallel --no-notice mapped_to_genome_distrib_length {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: mapped_to_genome_distrib_length An error has occurred."
fi
wait

parallel --no-notice multimapped_to_genome_distrib_length {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: multimapped_to_genome_distrib_length An error has occurred."
fi
wait

parallel --no-notice sam_to_bam {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: sam_to_bam An error has occurred."
fi
wait

#parallel --no-notice p_offset_analysis ::: $CONDITION_ARRAY_UNIQ ::: $WORKING_CONDITION_ARRAY ::: $WORKING_SAMPLE_ARRAY
p_offset_analysis
if [ $? -ne 0 ]
then
  error_exit "$LINENO: p_offset_analysis An error has occurred."
fi
wait

parallel --no-notice --xapply rna_seq_quantification {.} ::: $WORKING_SAMPLE_ARRAY ::: $SALMON_LIBTYPE
if [ $? -ne 0 ]
then
  error_exit "$LINENO: rna_seq_quantification An error has occurred."
fi
wait

cds_range_building
if [ $? -ne 0 ]
then
  error_exit "$LINENO: cds_range_building An error has occurred."
fi
wait

isoform_level_estimation
if [ $? -ne 0 ]
then
  error_exit "$LINENO: isoform_level_estimation An error has occurred."
fi
wait

sam_to_bam_seqcluster
if [ $? -ne 0 ]
then
  error_exit "$LINENO: sam_to_bam_seqcluster An error has occurred."
fi
wait

clustering_step_seqcluster
if [ $? -ne 0 ]
then
  error_exit "$LINENO: clustering_step_seqcluster An error has occurred."
fi
wait

report_step_seqcluster
if [ $? -ne 0 ]
then
  error_exit "$LINENO: report_step_seqcluster An error has occurred."
fi
wait

get_longest_transcripts_from_annotations
if [ $? -ne 0 ]
then
  error_exit "$LINENO: get_longest_transcripts_from_annotations An error has occurred."
fi
wait

parallel --no-notice htseq_count {.} ::: $WORKING_SAMPLE_ARRAY
if [ $? -ne 0 ]
then
  error_exit "$LINENO: htseq_count An error has occurred."
fi
wait

build_rnaseq_ribopro_counting_tables
if [ $? -ne 0 ]
then
  error_exit "$LINENO: build_rnaseq_ribopro_counting_tables An error has occurred."
fi
wait

anadif_babel
if [ $? -ne 0 ]
then
  error_exit "$LINENO: anadif_babel An error has occurred."
fi
wait

anadif_sartools $PROJECT_NAME $AUTHOR $REFERENCE_CONDITION
if [ $? -ne 0 ]
then
  error_exit "$LINENO: anadif_sartools An error has occurred."
fi
wait

anadif_clusters "${PROJECT_NAME}_multimapped_analysis" $AUTHOR $REFERENCE_CONDITION
if [ $? -ne 0 ]
then
  error_exit "$LINENO: anadif_clusters An error has occurred."
fi
wait

# Write final report
#TODO DEFINITELY ,ake something about that!
#FINALLOGFILE="${PROJECT_NAME}.final.report"
#for file in $(ls -c *log); do stat -c '%y' $file >> $FINALLOGFILE; printf "\n" >> $FINALLOGFILE; cat $file >> $FINALLOGFILE; done

# Put log files in log directory
mkdir -p log
chown $(id -u):$(id -g) -R log
mv *.log log
echo "End of the analysis, log files in log/ directory."
