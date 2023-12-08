---
title: Nextflow coding practices
teaching: 30
exercises: 15
---

::::::::::::::::::::::::::::::::::::::: objectives

- Learn how to use whitespace and comments to improve code readability.
- Understand coding pitfalls that reduce portability.
- Understand coding pitfalls that reduce maintainability.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I make my code readable?
- How do I make my code portable?
- How do I make my code maintainable?

::::::::::::::::::::::::::::::::::::::::::::::::::

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

```groovy
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
    FOO ( Channel.fromPath( params.reads, checkIfExists: true ) )
    BAR ( FOO.out )
    // Tip: Use vertical spacing and indentation for many parameters.
    BAZ (
        Channel.fromPath( params.gene_list, checkIfExists: true ),
        FOO.out,
        BAR.out,
        file( params.gene_db, checkIfExists: true )
    )

}

// Tip: Uppercase process names help readability
process FOO {

    // Tip: Separate process parts into distinct blocks
    input:
    path fastq

    output:
    path "*.fasta"

    script:
    prefix = fastq.baseName
    """
    tofasta $fastq > $prefix.fasta
    """
}
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Improve the workflow readability

Use whitespace to improve the readability of the following code.

```groovy
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
path "*.fasta"
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
```

:::::::::::::::  solution

## Solution

```groovy
#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads = ''

workflow {
    FOO ( Channel.fromPath( params.reads ) )
    BAR ( FOO.out )
}

process FOO {

    input:
    path fastq

    output:
    path "*.fasta"

    script:
    prefix = fastq.baseName
    """
    tofasta $fastq > $prefix.fasta
    """
}

process BAR {

    input:
    path fasta

    script:
    """
    fastx_check $fasta
    """
}
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Use comments

Comments are an important tool to improve readability and maintenance.
Use them to:

- Annotate data structures expected in a channel.
- Describe higher level functionality.
- Describe presence/absence of (un)expected code.
- Mandatory and optional process inputs.

```groovy
workflow ALIGN_SEQ {

    take:
    reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
    reference    // file( "path/to/reference" )

    main:
    // Quality Check input reads
    READ_QC ( reads )

    // Align reads to reference
    Channel.empty()
        .set { aligned_reads_ch }
    if( params.aligner == 'hisat2' ){
        ALIGN_HISAT2 ( READ_QC.out.reads, reference )
        aligned_reads_ch.mix( ALIGN_HISAT2.out.bam )
            .set { aligned_reads_ch }
    } else if ( params.aligner == 'star' ) {
        ALIGN_STAR ( READ_QC.out.reads, reference )
        aligned_reads_ch.mix( ALIGN_STAR.out.bam )
            .set { aligned_reads_ch }
    }
    aligned_reads_ch.view()

    emit:
    bam = aligned_reads_ch   // queue channel: [ sample_id, file(bam_file) ]

}

process COUNT_KMERS {

    input:
    // Mandatory
    tuple val(sample), path(reads)  // [ 'sample_id', [ read1, read2 ] ]: Reads in which to count kmers
    // Optional
    path kmer_table                 // 'path/to/kmer_table': Table of k-mers to count

    ...
}
```

### Report tool versions

Software packaging is a hard problem, and it can be difficult for
a package to report the versions of all the tools it has. It
may also be excessive to report the version of everything included
in a package, when only a handful of tools are used. This means
that it's up to us to effectively report the versions of the tools
we use to aid reproducibility.

```groovy
process HISAT2_ALIGN {

    ...

    script:
    def HISAT2_VERSION = '2.2.0' // Version not available using command-line
    """
    hisat2 ... | samtools ...

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: $HISAT2_VERSION
        samtools: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
```

### Name output channels

Output channels from processes and workflows can be named
using the `emit:` keyword, which helps readability.

```groovy
workflow ALIGN_HISAT2 {

    ...

    emit:
    alignment = HISAT2_ALIGN.out.bam

}

process HISAT2_ALIGN {

    ...

    output:
    tuple val(sample), path("*.bam"), emit: bam
    tuple val(sample), path("*.log"), emit: summary
    path "versions.yml"             , emit: versions

    ...
}
```

### Use params.parameters in workflow blocks, not in process blocks

The `params` variables are accessible from anywhere in a workflow.
They can be useful to provide a wide variety of properties and
decision making options. For example, one could use a `params.aligner`
variable in a workflow to select a particular alignment tool. This
in turn could be coded like:

```groovy
process ALIGN {

    input:
    tuple val(sample), path(reads)
    path index

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
```

A better practice is to use it as an input value.

```groovy
process ALIGN {

    input:
    tuple val(sample), path(reads)
    path index
    val aligner       // params.aligner is provided as a third parameter

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
```

This allows one to see from the `workflow` block where all
parameters are being used, making the workflow easier to maintain.
There is also a danger that one could modify `params` variables during
pipeline execution, potentially leading to unreproducible results
in more complex workflows.

### All input files/directories should be a process input

Depending on the platform you execute your workflow, files
may be easily accessible over the network, or downloadable
from the internet. However not all execution platforms
support this. The example below could work well on your
system, but fail on another (for example compute nodes without
internet connection).

```groovy
process READ_CHECK {

    input:
    tuple val(sample), path(reads)

    ...

    script:
    """
    wget ftp://path/to/database

    check_reads $reads /local/copy/database > $sample.report
    ...
    """
}
```

A strength of Nextflow is file staging, i.e.,
preparing files for use in process tasks.
Staging files by providing them as process input has several benefits.

- Files are updated only when necessary.
- A single file/folder can be shared, without downloading multiple
  times as in the example above.
- Nextflow supports retrieving files from any valid URL, meaning
  potentially fewer lines of code.

```groovy
process READ_CHECK {

    input:
    tuple val(sample), path(reads)
    path database

    ...

    script:
    """
    check_reads $reads $database > $sample.report
    ...
    """
}
```

If may be that a file is an optional input depending on other parameters.
In cases when no file should be provided, one can pass an empty list `[]`
instead.

```
workflow {

    COUNT_KMERS ( reads, [] )
}

process COUNT_KMERS {

    input:
    // Mandatory
    tuple val(sample), path(reads)  // [ 'sample_id', [ read1, read2 ] ]: Reads in which to count kmers
    // Optional
    path kmer_table                 // 'path/to/kmer_table': Table of k-mers to count

    ...
}
```

### Avoid lots of short running processes

Many execution platforms are inefficient if a workflow
tries to execute many short running processes. It
can take more time to schedule and request resources
for each small instance than bundling the short processes
into a larger process task. Nextflow provides
convenient channel operators, such as `buffer`, `collate`,
`collect`, and `collectFile`,
that help group together inputs into batches which
can run for longer with a given requested resource.
The short tasks themselves can also be parallelised
inside a process script using the command-line tools
`xargs` or `parallel`.

```groovy
workflow REFINE_DATA {

    take:
    datapoints

    main:
    BATCH_TASK ( datapoints.collate(100) )
}

process BATCH_TASK {

    input:
    val data

    script:
    """
    printf "%s\\n" $data | \
        xargs -P $task.cpus -I {} \
        short_task {}

    # Alternative:
    # parallel --jobs $task.cpus "short_task {1}" :::: $data
    """
}
```

### Include a test profile

A `test` profile is a configuration profile that specifies
a short running test data set to check the functionality of
the whole pipeline. It can also demonstrate to users of your
workflow the kinds of inputs and outputs to expect. Another
benefit is the possibility of automated testing of your
workflow, ensuring the workflow keeps working as you add
new functionality.

```groovy
profiles {
    test {
        params {
           reads = 'https://github.com/my_repo/test/test_reads.fastq.gz'
           reference = 'https://github.com/my_repo/test/test_reference.fasta.gz'
        }
    }
}
```

### Write modules that use existing containers

Using containers for software packaging is strongly recommended,
as they are intended to operate the same, irrespective of the
operating system it runs on. Writing modules which use
existing containers reduces maintenance needed for a workflow,
and minimises the need to resolve package conflicts. A
helpful resource for this is the bioconda channel for the
package manager conda, which provides packaging for many
bioinformatics tools. In addition to this,
[Biocontainers](https://biocontainers.pro/) builds both
Docker and Singularity images for each tool packaged on
the bioconda channel. Multi-package containers (known
as mulled containers) can also be created following the
instructions on the [Multi Package Containers repository](https://github.com/BioContainers/multi-package-containers).

```
process FASTQC {

    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

    ...
}
```

Building your own container images should be used as a last resort.
A preferred option is to write a conda recipe for the tool
to be included in the bioconda channel. This makes the tool available
via a package manager, and containers are automatically built
for the tool.

### Use file compression and temporary disk space when possible

Disk space is often a valuable resource on most compute systems.
When possible, work with compressed files.
There are several useful shell commands and operations
that work well with compressed data, such as `gunzip -c` and `zgrep`.
Operations such as command substitution (`$( command )`), or
process substitution (`>( command_list )`, or `<( command_list)`)
can also help working with compressed data. Lastly,
named pipes can also be used if the above approaches fail.

```groovy
process COUNT_FASTA {

    input:
    path fasta         // reference.fasta.gz

    script:
    """
    zgrep -c '^>' $fasta
    """
}
```

Another good practice is to use local temporary disk space (also
known as scratch space). Often, the workDir is located on a
shared disk space over a network, which slows down processes
that read and write a lot to disk. Using scratch space not only
speeds up disk I/O, but also saves space in the workDir since only
files which match the process output directive are copied back
across for caching. The process directive `process.scratch`
can be provided with either a boolean or the path to use
for scratch space.

```groovy
process {
    scratch = '/tmp'
}
```

### Use consistent naming conventions

Using consistent naming conventions greatly helps readability.
For example using uppercase for process names helps
distinguish them from other workflow components like
channels or operators. Here are other suggestions one
can follow from [Nf-core](https://nf-co.re/developers/adding_modules#naming-conventions).



:::::::::::::::::::::::::::::::::::::::: keypoints

- Nextflow is not sensitive to whitespace. Use it to layout code for readability.
- Use comments and whitespace to group chunks of code to describe big picture functionality.
- Report tool versions in the scripts.
- Name channel outputs using the `emit:` keyword.
- Avoid `params.parameter` in a process. Pass all parameters using input channels.
- Input files should be passed using input channels.
- Group short running commands into a larger process.
- Include a test profile which runs the workflow on a small test data set.
- Write your processes to reuse existing containers/software bundles.
- Use compressed files and temporary disk space when possible.
- Use consistent naming conventions.

::::::::::::::::::::::::::::::::::::::::::::::::::


