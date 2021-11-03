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
include { ALIGN_STAR   } from 'workflows/align_star'

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
    } else if ( params.aligner == 'star' ) {
        ALIGN_STAR( READ_QC.out.reads, reference )
        aligned_reads_ch.mix( ALIGN_STAR.out.bam )
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
    def HISAT2_VERSION = '2.2.0' // Version not available using command-line
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

The `params` variables are accessible from anywhere in a workflow.
They can be useful to provide a wide variety of properties and
decision making options. For example, one could use a `params.aligner`
variable in a workflow to select a particular alignment tool. This
in turn could be coded like:

~~~
process ALIGN {

    input:
    tuple val(sample), path(reads)
    path index

    output:
    tuple val(sample), path("*.bam"), emit: bam
    ...

    script:
    if ( params.aligner == 'hisat2' ) {
        """
        hisat2 ... | samtools ...

        ...
        """
    } else if ( params.aligner == 'star' ){
        """
        star ...

        ...
        """
    }
}
~~~
{: .language-groovy}

A better practice is to use it as an input value.
~~~
process ALIGN {

    input:
    tuple val(sample), path(reads)
    path index
    val aligner       // params.aligner is provided as a third parameter

    output:
    tuple val(sample), path("*.bam"), emit: bam
    ...

    script:
    if ( aligner == 'hisat2' ) {
        """
        hisat2 ... | samtools ...

        ...
        """
    } else if ( aligner == 'star' ){
        """
        star ...

        ...
        """
    }
}
~~~
{: .language-groovy}

This allows one to see from the `workflow` block where all
parameters are being used, making the workflow easier to maintain.
There is also a danger that one could modify `params` variables during
pipeline execution, potentially leading to unreproducible results
in more complex workflows.

### All input files/directories should be a process input

Depending on the platform you execute your workflow, files
may be easily accessible over the network, or downloadable
from the internet. However not all execution platforms
support this. A strength of Nextflow is file staging, i.e.,
preparing files for use in process tasks. It could be
fine to code as below on your execution platform.

~~~
process DO_SOMETHING {

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.report"), emit: report
    ...

    script:
    """
    wget ftp://path/to/database

    check_reads $reads /local/copy/database > $sample.report
    ...
    """
}
~~~
{: .language-groovy}

Staging files by providing them as process input has several benefits.
- Files are updated only when necessary.
- A single file/folder can be shared, without downloading multiple
times as in the example above.
- Nextflow supports retrieving files from any valid URL, meaning
potentially few lines of code.

~~~
process DO_SOMETHING {

    input:
    tuple val(sample), path(reads)
    path database

    output:
    tuple val(sample), path("*.report"), emit: report
    ...

    script:
    """
    check_reads $reads $database > $sample.report
    ...
    """
}
~~~
{: .language-groovy}

### Avoid lots of short running processes

Many execution platforms are inefficient if a workflow
tries to execute many short running processes. It
can take more time to schedule and request resources
for each small instance than bundling the short processes
into a larger process task. Nextflow provides some
convenient channel operators, such as `buffer` and `collate`
that can help group together inputs into batches that
can run for longer with a given requested resource.

~~~
workflow REFINE_DATA {
    take:
    datapoints

    main:
    BATCH_TASK( datapoints.collate(100) )
}

process BATCH_TASK {

    input:
    val data

    script:
    """
    parallel --jobs $task.cpus "short_task {1}" :::: $data
    """
}
~~~
{: .language-groovy}

{% include links.md %}
