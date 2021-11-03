---
title: "Nextflow coding practices"
teaching: 30
exercises: 15
questions:
- "How do I make my code readable?"
- "How do I make my code portable?"
- "How do I make my code maintainable?"
objectives:
- "Learn how to use whitespace and comments to improve code readability."
- "Understand coding pitfalls that reduce portability."
- "Understand coding pitfalls that reduce maintainability."
keypoints:
- "Nextflow is not sensitive to whitespace. Use it to layout code for readability."
- "Use comments and whitespace to group chunks of code to describe big picture functionality."
- "Avoid `params.parameter` in a process. Pass all parameters using input channels."
- "Input files should be passed using input channels."
---

## Nextflow coding practices

Nextflow is a powerful flexible language that one can code in a variety of ways.
This can lead to poor practices in coding. For example, this can lead
to the workflow only working under certain configurations or execution platforms.
Alternatively, it can make it harder for someone to contribute to a codebase,
or for you to amend two years later for article submission.
These are some useful coding tips that make maintaining and porting your
workflow easier.

### Use whitespace to improve readability.

Nextflow is generally not sensitive to whitespace in code. This allows you
to use indentation, vertical spacing, new-lines, and increased spacing to
improve code readability.

~~~
#! /usr/bin/env nextflow

// Tip: Allow spaces around assignments ( = )
nextflow.enable.dsl = 2     

// Tip: Separate blocks of code into groups with common purpose
//      e.g., parameter blocks, include statements, workflow blocks, process blocks
// Tip: Align assignment operators vertically in a block
params.reads     = ''
params.gene_list = ''
params.gene_db   = 'ftp://path/to/database'

// Tip: Align braces or instruction parts vertically
include { BAR        } from 'modules/bar'
include { TAN as BAZ } from 'modules/tan'

workflow {

    // Tip: Indent process calls
    // Tip: Use spaces around process/function parameters
    FOO( Channel.fromPath( params.reads, checkIfExists: true ) )
    BAR( FOO.out )
    // Tip: Use vertical spacing and indentation for many parameters.
    BAZ(
        Channel.fromPath( params.gene_list, checkIfExists:true ),
        FOO.out,
        BAR.out,
        file( params.gene_db, checkIfExists: true )
    )

}

process FOO {

    // Tip: Separate process parts into distinct blocks
    input:
    path fastq

    output:
    path

    script:
    prefix = fastq.baseName
    """
    tofasta $fastq > $prefix.fasta
    """
}
~~~
{: .language-groovy}

> ## Challenge
> Use whitespace to improve the readability of the following code.
~~~
#! /usr/bin/env nextflow

nextflow.enable.dsl=2
params.reads = ''
workflow {
foo(Channel.fromPath(params.reads))
bar(foo.out)
}
process foo {
input:
path fastq
output:
path
script:
prefix=fastq.baseName
"""
tofasta $fastq > $prefix.fasta
"""
}
process bar {
input:
path fasta
script:
"""
fastx_check $fasta
"""
}
~~~
{: .language-groovy}

### Use comments

Comments are an important tool to improve readability and maintenance.
Use them to:
- Annotate data structures expected in a channel.
- Describe higher level functionality.
- Describe presence/absence of (un)expected code.

~~~
#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { READ_QC      } from 'workflows/read_qc'
include { ALIGN_HISAT2 } from 'workflows/align_hisat2'
include { ALIGN_SALMON } from 'workflows/align_salmon'

workflow {

    RNA_SEQ(
        Channel.fromFilePairs( params.reads , checkIfExists: true ),
        file( params.reference, checkIfExists: true )
    )
}

workflow RNA_SEQ {

    take:
    reads        // Queue type; Data: [ sample_id, [ file(read1), file(read2) ] ]
    reference    // Value type; file( "path/to/reference" )

    main:
    // Quality Check input reads
    READ_QC( reads )

    // Align reads to reference
    Channel.empty()
        .set { aligned_reads_ch }
    if( params.aligner == 'hisat2' ){
        ALIGN_HISAT2( READ_QC.out.reads, reference )
        aligned_reads_ch.mix( ALIGN_HISAT2.out.bam )
    } else if ( params.aligner == 'salmon' ) {
        ALIGN_SALMON( READ_QC.out.reads, reference )
        aligned_reads_ch.mix( ALIGN_SALMON.out.bam )
    }
    aligned_reads_ch.view()

}
~~~
{: .language-groovy}

### Report tool versions

Software packaging is a hard problem, and it can be difficult for
a package to report the versions of all the tools it has. It
may also be excessive to report the version of everything included
in a package, when only a handful of tools are used. This means
that it's up to us to effectively report the versions of the tools
we use to aid reproducibility.

~~~
process HISAT2_ALIGN {

    input:
    tuple val(sample), path(reads)
    path index

    output:
    tuple val(sample), path("*.bam"), emit: bam
    tuple val(sample), path("*.log"), emit: summary
    path "versions.yml"             , emit: versions

    script:
    def HISAT2_VERSION = '2.2.0' // Version not available using command-line
    """
    hisat2 ... | samtools ...

    cat <<-END_VERSIONS > versions.yml
    HISAT2_ALIGN:
        hisat2: $HISAT2_VERSION
        samtools: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
~~~
{: .language-groovy}

### Name output channels

Output channels from processes and workflows can be named,
which helps readability.

~~~
workflow ALIGN_HISAT2 {

    take:
    reads
    reference

    main:
    HISAT2_INDEX( reference )
    HISAT2_ALIGN( reads, HISAT2_INDEX.out.index )

    emit:
    alignment = HISAT2_ALIGN.out.bam

}

process HISAT2_ALIGN {

    input:
    tuple val(sample), path(reads)
    path index

    output:
    tuple val(sample), path("*.bam"), emit: bam
    tuple val(sample), path("*.log"), emit: summary
    path "versions.yml"             , emit: versions

    script:
    def VERSION = '2.2.0' // Version not available using command-line
    """
    hisat2 ... | samtools ...

    cat <<-END_VERSIONS > versions.yml
    HISAT2_ALIGN:
        hisat2: $HISAT2_VERSION
        samtools: \$(samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
~~~
{: .language-groovy}

### Use params.parameters in workflow blocks, not in process blocks

### All input files/directories should be a process input

### Avoid lots of short running process


{% include links.md %}
