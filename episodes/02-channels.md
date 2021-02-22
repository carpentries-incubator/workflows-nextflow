---
title: "Channels"
teaching: 15
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Channels

Channels are a key data structure of Nextflow that allows the implementation of reactive-functional oriented computational workflows based on the [Dataflow programming paradigm](https://en.wikipedia.org/wiki/Dataflow_programming).

They are used to logically connect tasks each other or to implement functional style data transformations.


## Channel types

Nextflow distinguish two different kinds of channels: queue channels and value channels.

### Queue channel

A queue channel is a asynchronous unidirectional FIFO queue which connects two processes or operators.

* What asynchronous means? That operations are non-blocking.

* What unidirectional means? That data flow from a producer to a consumer.

* What FIFO means? That the data is guaranteed to be delivered in the same order as it is produced.

A queue channel is implicitly created by process output definitions or using channel factories methods such as [Channel.from](https://www.nextflow.io/docs/latest/channel.html#from) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).

Try the following snippets:


> ## View Channel contents
> ~~~
> ch = Channel.from(1,2,3)
> ~~~
> 	Use the built-in println function to print the ch variable.
>	  Apply the view method to the ch channel, therefore prints each item emitted by the channels.
>
> > ## Solution
> > ~~~
> > println(ch)  
> > ch.view() 
> > ~~~
> {: .solution}
{: .challenge}

> Queue channels usage
> A queue channel can have one and exactly one producer and one and exactly one consumer.
{: .callout}


### Value channels

A value channel a.k.a. singleton channel by definition is bound to a single value and it can be read unlimited times without consuming its content.

~~~
ch = Channel.value('Hello')
ch.view()
ch.view()
ch.view()
~~~
{: .source}

It prints:
~~~
Hello
Hello
Hello
~~~
{: .output}


## Channel factories

### value

The value factory method is used to create a value channel. An optional not null argument can be specified to bind the channel to a specific value. For example:

~~~
ch1 = Channel.value() 
ch2 = Channel.value( 'Hello there' )
ch2 = Channel.value( [1,2,3,4,5] )
~~~
{: .source}

* Creates an empty value channel.
* Creates a value channel and binds a string to it.
* Creates a value channel and binds a list object to it that will be emitted as a sole emission.

### from

The factory Channel.from allows the creation of a queue channel with the values specified as argumThe factory Channel.from allows the creation of a queue channel with the values specified as argument.


~~~
ch = Channel.from( 1, 3, 5, 7 )
ch.view{ "value: $it" }
~~~
{: .source}

The first line in this example creates a variable ch which holds a channel object. This channel emits the values specified as a parameter in the from method. Thus the second line will print the following:

~~~
value: 1
value: 3
value: 5
value: 7
~~~
{: .output}

> Method Channel.from will be deprecated and replaced by Channel.of (see below).
{: .callout}


### of
The method Channel.of works in a similar manner to Channel.from, though it fixes some inconsistent behavior of the latter and provides a better handling for range of values. For example:

~~~
Channel
    .of(1..21, 'X', 'Y')
    .view()
~~~
{: .source}

### fromList

The method Channel.fromList creates a channel emitting the elements provided by a list objects specified as argument:

~~~
list = ['hello', 'world']

Channel
    .fromList(list)
    .view()
~~~
{: .source}


### fromPath
The fromPath factory method create a queue channel emitting one or more files matching the specified glob pattern.

~~~
Channel.fromPath( '/data/big/*.txt' )
~~~
{: .source}

This example creates a channel and emits as many items as there are files with txt extension in the /data/big folder. Each element is a file object implementing the [Path](https://docs.oracle.com/javase/8/docs/api/java/nio/file/Paths.html) interface.

> Two asterisks, i.e. **, works like * but crosses directory boundaries. This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.
{: .callout}






|Name|	Description|
|-----|----------|
|glob |When true interprets characters *, ?, [] and {} as glob wildcards, otherwise handles them as normal characters (default: true)|
|type | Type of paths returned, either file, dir or any (default: file) |
| hidden | When true includes hidden files in the resulting paths (default: false)|
| maxDepth | Maximum number of directory levels to visit (default: no limit) |
| followLinks | When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true) |
| relative | When true returned paths are relative to the top-most common directory (default: false) |
| checkIfExists | When true throws an exception of the specified path do not exist in the file system (default: false) |

When true throws an exception of the specified path do not exist in the file system (default: false)

> ## Challenge Title
>
>  Use the Channel.fromPath method to create a channel emitting all files with the suffix .fq in the data/ggal/ and any subdirectory, then print the file name.
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}


### fromFilePairs

The fromFilePairs method creates a channel emitting the file pairs matching a glob pattern provided by the user. The matching files are emitted as tuples in which the first element is the grouping key of the matching pair and the second element is the list of files (sorted in lexicographical order).

~~~
Channel
    .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
    .view()
~~~
{: .source}

It will produce an output similar to the following:

~~~
[SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
[SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
[SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
[SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
[SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
[SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]
~~~
{: .output}


> The glob pattern must contain at least a star wildcard character.
{: .callout}




fixme add table

> ## Challenge Title
>
>  Use the fromFilePairs method to create a channel emitting all pairs of fastq read in the data/ggal/ directory and print them.
>  Then use the flat:true option and compare the output with the previous execution.
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

### fromSRA

The Channel.fromSRA method that makes it possible to query of [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) archive and returns a channel emitting the FASTQ files matching the specified selection criteria.

The query can be project ID or accession number(s) supported by the [NCBI ESearch API](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch). For example the following snippet:

~~~
Channel
    .fromSRA('SRP043510')
    .view()
~~~
{: .source}

prints:

~~~
[SRR1448794, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448794/SRR1448794.fastq.gz]
[SRR1448795, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/005/SRR1448795/SRR1448795.fastq.gz]
[SRR1448792, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/002/SRR1448792/SRR1448792.fastq.gz]
[SRR1448793, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/003/SRR1448793/SRR1448793.fastq.gz]
[SRR1910483, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/003/SRR1910483/SRR1910483.fastq.gz]
[SRR1910482, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/002/SRR1910482/SRR1910482.fastq.gz]
(remaining omitted)
~~~
{: .output}

> Read pairs are implicitly managed are returned as a list of files.

It’s straightforward to use this channel as an input using the usual Nextflow syntax. For example:

~~~
params.accession = 'SRP043510'
reads = Channel.fromSRA(params.accession)

process fastqc {
    input:
    tuple sample_id, file(reads_file) from reads

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
}
~~~

The code snippet above creates a channel containing 24 samples from a chromatin dynamics study and runs FASTQC on the resulting files.


{% include links.md %}

