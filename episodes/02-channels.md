---
title: "Channels"
teaching: 30
exercises: 10
questions:
- "What are Nextflow channels?"
- "Define the different types of Nextflow channels?"
- "What are the major differences between queue and values channels?"
- "How do you create a channel?"
- "How do you create a queue channel from a specified glob pattern?"
- "How do you modify the behaviour of a channel factory?"
-
objectives:
- "Understand what a Nextflow Channel is."
- "Understand the differences between value and queue channels."
- "Create a value and queue channel using Channel factory methods."
- "Create a queue channel from a specified glob pattern."
- "Modify the behaviour of a Channel factory"
keypoints:
- "Channels are a key data structure of Nextflow that allows the implementation of reactive-functional oriented computational workflows"
- "Nextflow distinguish two different kinds of channels, queue channels and value channels. Values channels contain a  single value and can be used multiple times in workflow. Queue channels can be used once are consumed when they are used by a process."
- "Channel factory methods are used to create value and queue channels."
- "You can see the contents of a channel using the `.view()` operator"
---

# Channels

In the last episode we learnt that channels are the way in which Nextflow connect workflow tasks `process`.

Channels are a key data structure of Nextflow that allows the implementation of [reactive-functional](https://en.wikipedia.org/wiki/Functional_reactive_programming) oriented computational workflows based on the [Dataflow programming paradigm](https://en.wikipedia.org/wiki/Dataflow_programming).

They are used to logically connect tasks/processes to each other.

> Reactive programming is often explained with an analogy to a spreadsheet: Imagine a cell that calculates the input of two other cells. Once you change one of the inputs, the sum updates as well. The cell reacts to the changes and updates itself.
This is very similar to dataflow programming. Conceptually, the focus here lies on the flow of the data instead of on the flow of control.[more here](https://itnext.io/demystifying-functional-reactive-programming-67767dbe520b)
{: .callout}

![Channel files](../fig/channel-files.png)

Maybe use this "Nextflow is based on the Dataflow programming model in which processes communicate through channels. from [here](https://www.nextflow.io/docs/latest/channel.html#)

## Channel types

Nextflow distinguish two different kinds of channels: **queue** channels and **value** channels.

### Queue channel

A *queue* channel is a *asynchronous* unidirectional FIFO queue which connects two `processes` or `operators`.

* What *asynchronous* means? That operations are non-blocking.

* What unidirectional means? That data flow from a producer to a consumer.

* What *FIFO* means? That the data is guaranteed to be delivered in the same order as it is produced.

A queue channel is  created by process output definitions, which we will cover in the next episode, or using channel factories methods such as [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).


> ## Create and view Channel contents
> Create a file called `channel.nf` and type the following code into it.
> ~~~
> ch = Channel.of(1,2,3)
> ch.view()
> ~~~
> Run the code using
> ~~~~
> nextflow run channel.nf
> ~~~~
> > ## Solution
> > In this example you have created a queue channel with the value 1,2,3 in it.
> > The `.view` operator can be used to view the contents of the  ch channel.
> > ~~~
> > 1
> > 2
> > 3
> > ~~~
> {: .solution}
{: .challenge}

### Queue channels usage

Once a channel is used by an operator or process it can not be used again.

~~~
ch = Channel.of(1,2,3)
ch.view()
ch.view()
~~~
{: .source}


Produces an error as the channel `ch` has been used twice.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `channel.nf` [chaotic_payne] - revision: 63c4f7a5f3
1
2
3
Channel `ch` has been used as an input by more than a process or an operator

 -- Check script 'main.nf' at line: 3 or see '.nextflow.log' file for more details
~~~
{: .output}

We will see in operator episodes how to create multiple channels from the data.

### Value channels

The second type of Nextflow channel is a `value` channel. A **value** channel by is bound to a *single* value (singleton). A value channel can read unlimited times without consuming its content. This can be useful when setting a single value which is used by multiple processes.

~~~
ch = Channel.value('GRCh38')
ch.view()
ch.view()
ch.view()
~~~
{: .source}

It prints:

~~~
GRCh38
GRCh38
GRCh38
~~~
{: .output}


> ## Queue vs Value Channel
>
> Which of the following channels could be used multiple times and why?
>
> 1. `ch = Channel.value('GRCh38')`
> 2. `ch = Channel.of('GRCh38')`
>
> > ## Solution
> > 1. Yes: `ch = Channel.value('GRCh38')` is a value channel which can be used multiple times.
> > 2. No: `ch = Channel.of('GRCh38')` is a queue channel which can only be used once.
> {: .solution}
{: .challenge}

## Creating Channels, Channel factories

Channel factories are the way Nextflow creates the different channel types (queue and value).
A factory is the term used when creating an object of a different type e.g. value or queue channels.

### value

The `value` factory method is used to create a value channel.
An optional  argument can be specified to bind the channel to a specific value. For example:

~~~
ch1 = Channel.value()
ch2 = Channel.value('GRCh38')
ch2 = Channel.value( ['chr1','chr2','chr3','chr4','chr5'] )
~~~
{: .source}

1. Creates an empty value channel.
1. Creates a value channel and binds a string to it.
1. Creates a value channel and binds a list object to it that will be emitted as a sole emission.

### of

The factory `Channel.of` allows the creation of a `queue` channel with the values specified as argument.


~~~
ch = Channel.of( 'chr1','chr3','chr5','chr7' )
ch.view()
~~~
{: .source}


The first line in this example creates a variable `ch` which holds a queue channel object.
This channel emits the values specified as a parameter in the `of` method. Thus the second line will print the following:

~~~
chr1
chr3
chr5
chr7
~~~
{: .output}

~~~
ch= Channel.of(1..21, 'X', 'Y')
ch.view()
~~~
{: .source}

> ## Channel.from
> You may see the method `Channel.from` in older nextflow scripts, this performs a similar function but will is deprecated so you should use `Channel.of` instead.
{: .callout}


### fromList

The method `Channel.fromList` creates a channel emitting the elements in a list objects specified as argument. A List object can be defined by placing the items in square brackets `[]` separated by a comma.

~~~
aligner_list = ['salmon', 'kallisto']

aligner_ch = Channel.fromList(aligner_list)

aligner_ch.view()
~~~
{: .source}

This would produce.

~~~
salmon
kallisto
~~~
{: .source}

> ## Channel.fromList vs Channel.of
> In the above example the channel has two elements. If you has used the Channel.of(list) it would have  contained only 1 element `[salmon, kallisto]`.
{: .callout}

### fromPath

The previous channel factories deal with sending values. A special  channel factory `fromPath` is used when dealing with files.

The `fromPath` factory method create a **queue channel** emitting one or more files matching the specified  pattern (glob) specifying the location of files.

~~~
ch = Channel.fromPath( 'data/ggal/*.fq' )
ch.view()
~~~
{: .source}


This example creates a channel and emits as many items as there are files with `fq` extension in the `/data/ggal` folder.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [disturbed_pike] - revision: 34d3884f06
/home/ec2-user/environment/data/ggal/lung_2.fq
/home/ec2-user/environment/data/ggal/lung_1.fq
/home/ec2-user/environment/data/ggal/gut_1.fq
/home/ec2-user/environment/data/ggal/liver_1.fq
/home/ec2-user/environment/data/ggal/liver_2.fq
/home/ec2-user/environment/data/ggal/gut_2.fq
~~~
{: .output}

> ## asterisks
> Two asterisks, i.e. `**`, works like * but crosses directory boundaries. This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns.
{: .callout}

We can change the behaviour of Channel.fromPath by changing the options.

Table 1. Available options

|Name|Description|
|-----|----------|
|glob |When true interprets characters *, ?, [] and {} as glob wildcards, otherwise handles them as normal characters (default: true)|
|type | Type of paths returned, either file, dir or any (default: file) |
| hidden | When true includes hidden files in the resulting paths (default: false)|
| maxDepth | Maximum number of directory levels to visit (default: no limit) |
| followLinks | When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true) |
| relative | When true returned paths are relative to the top-most common directory (default: false) |
| checkIfExists | When true throws an exception of the specified path do not exist in the file system (default: false) | When true throws an exception of the specified path do not exist in the file system (default: false)|

Learn more about the glob patterns syntax at this [link](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob).


We can change the default options to give an error if the file doesn't exist using the `checkIfExists` option. In Nextflow option parameters are seprated by a `,` and defined with a colon `:`.

~~~
ch = Channel.fromPath( 'data/chicken/*.fq', checkIfExists: true )
ch.view()
~~~
{: .source}

Gives an error as there is no data/chicken directory.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [intergalactic_mcclintock] - revision: d2c138894b
No files match pattern `*.fq` at path: /chicken/ggal/
~~~
{: .output}

> ## Using Channel.fromPath
>
>  Use the `Channel.fromPath` method to create a channel emitting all files with the suffix `.fq` in the `data/ggal/` and any subdirectory, then print the file name using the view operator.
>
> > ## Solution
> >
> > ~~~
> > ch = Channel.fromPath('data/data/ggal/*.fq')
> > ch.view()
> > ~~~
> >
> {: .solution}
{: .challenge}


### fromFilePairs

We have seen how to process files individually using `fromPath`. In Bioinformatics we often want to process files in pairs or larger groups such as read pairs in sequencing.

Nextflow provides helper method for these common bioinformatics use cases. The `fromFilePairs` method create a channel containing a named list (tuple) for every file matching a specific pattern.

The `fromFilePairs` method creates a channel emitting the file pairs matching a glob pattern provided by the user. The matching files are emitted as tuples in which the first element is the grouping key of the matching pair and the second element is the list of files (sorted in lexicographical order).

~~~
filepair_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')
filepair_ch.view()
~~~
{: .source}

It will produce an output channel of three elements, a value the glob pattern match and the set of file pairs (tuple):

~~~
[SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
[SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
[SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
[SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
[SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
[SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]
~~~
{: .output}

You can alter the behaviour of the `fromFilePairs` like `fromPath` by changing the options.

e.g.
~~~
filepair_ch = Channel.fromFilePairs('/my/data/SRR*_{1,2}.fastq', flat: true)
filepair_ch.view()
~~~
{: .source}

The option `flat:true`  will produce a list not a tuple.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [nice_kimura] - revision: 7ba6ecfed9
[lung, /home/ec2-user/environment/data/ggal/lung_1.fq, /home/ec2-user/environment/data/ggal/lung_2.fq]
[liver, /home/ec2-user/environment/data/ggal/liver_1.fq, /home/ec2-user/environment/data/ggal/liver_2.fq]
[gut, /home/ec2-user/environment/data/ggal/gut_1.fq, /home/ec2-user/environment/data/ggal/gut_2.fq]
~~~
{: .output}
see more [here](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)

> ## More complex patterns
> If you need to match more complex patterns you should create a sample sheet specifying the files and create a channel from that.
{: .callout}


> ## glob pattern
> The glob pattern must contain at least a star wildcard character.
{: .callout}



|Name|Description|
|-----|----------|
|type|Type of paths returned, either file, dir or any (default: file)|
|hidden|When true includes hidden files in the resulting paths (default: false)|
|maxDepth|Maximum number of directory levels to visit (default: no limit)|
|followLinks|When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true)|
|size|Defines the number of files each emitted item is expected to hold (default: 2). Set to -1 for any.|
|flat|When true the matching files are produced as sole elements in the emitted tuples (default: false).|
|checkIfExists|When true throws an exception of the specified path do not exist in the file system (default: false)|

> ## fromFilePairs
>
>  Use the `fromFilePairs` method to create a channel emitting all pairs of fastq read in the `data/ggal/` directory and print them.
>  Then use the `flat:true` option and compare the output with the previous execution.
>
> > ## Solution
> >
> > ~~~
> > Channel.fromFilePairs('data/ggal/*_{1,2}.fq').view()
> > Channel.fromFilePairs('data/ggal/*_{1,2}.fq', flat:true).view()
> > ~~~
> >
> {: .solution}
{: .challenge}

### fromSRA

Another useful helper method is `Channel.fromSRA`. The `Channel.fromSRA` method that makes it possible to query of [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) archive and returns a channel emitting the FASTQ files matching the specified selection criteria.

The query can be project ID or accession number(s) supported by the [NCBI ESearch API](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch). For example the following snippet:

~~~
sra_ch =Channel.fromSRA('SRP043510')
sra_ch.view()
~~~
{: .source}

This will prints a line containing a named list ,tuple, for every fastq file associated with that accession.

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

Multiple accession IDs can be specified using a list object:

~~~
ids = ['ERR908507', 'ERR908506', 'ERR908505']
Channel.fromSRA(ids).view()
~~~
{: .source}

~~~
[ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
~~~
{: .output}  

> ## Read pairs
> Read pairs are implicitly managed are returned as a list of files.
{: .callout}

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
