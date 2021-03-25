#!/usr/bin/env nextflow

params.genome = "$baseDir/data/genome/*.fa" //TODO perhaps compile a regexp to capture the name of the organism.
params.tarnscriptome = "$baseDir/data/transcriptome/*.fa" //TODO same as above
params.reads = "$baseDir/data/reads/*.fq"
params.annotation = "$baseDir/data/*.{gtf, gff}" //TODO perhaps a more robust regexp
params.resultsDir = "$baseDir/results"

println = """\
R I B O P I P E an integrated pipeline for ribo/rna-seq best practices analyses
===============================================================================

Genome        : ${params.genome}
Transcriptome : ${params.transcriptome}
Reads         : ${params.reads}
Annotation    : ${params.annotation}
Results       : ${params.resultsDir}
"""
.stripIndent()


sequences = file(params.reads)

/* split a fasta file in multiple files */
process genomeAlignment {
    input:
    file 'input.fa' from sequences
    output:
    file 'seq_*' into records
    """
    STAR ...
    """
}

/* RSEM counts */
process reverse {
    input:
    file x from records
    output:
    stdout result
    """
    RSEM ...
    """
}
