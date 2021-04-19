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
- "Channels are how nextflow communicates data between processes. "
- "Nextflow distinguish two different kinds of channels, queue channels and value channels. Values channels contain a  single value and can be used multiple times in workflow. Queue channels can be used once are consumed when they are used by a process or an operator."
- "Channel factory methods are used to create value and queue channels."
- "Channel factory methods have optional arguments, e.g. `checkIfExists` , that can be used to alter the creation and behaviour of a channel. "
---

# Channels

In the last episode we learnt that channels are the way in which Nextflow sends data around a workflow. Channels connect processes to each other, via their inputs and outputs. Channels can store multiple items, such as files (e.g. fastq files) or values. The number of items a channel stores determines how many times a process runs. It is channels that enable Nextflow to implement reactive workflows based on the Dataflow programming paradigm.

> Reactive programming is often explained with an analogy to a spreadsheet: Imagine a cell that calculates the input of two other cells. Once you change one of the inputs, the sum updates as well. The cell reacts to the changes and updates itself. More info [here](https://itnext.io/demystifying-functional-reactive-programming-67767dbe520b)
{: .callout}

![Channel files](../fig/channel-files.png)


## Channel types

Nextflow distinguish two different kinds of channels: **queue** channels and **value** channels.

### Queue channel

Queue channels are the way in which Nextflow sends data that is consumed or produced by a process.
A *queue* channel has three properties.

* That operations are non-blocking, which means processes run as soon as they receive input from a channel.

* That data flow is in one direction from a producer to a consumer.

* The first element of a queue is the first out of the queue (First in- First out).

A queue channel can be created by a process in their output definitions, which we will cover in the next episode, or created using channel factory methods such as [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath) which will we cover in this episode.


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

Once a channel is used/consumed by an operator or process it can not be used again.

If we add another `ch.view()` operation on the same channel object.

~~~
ch = Channel.of(1,2,3)
ch.view()
ch.view()
~~~
{: .source}

Then run using :

~~~
nextflow run channel.nf
~~~
{: .language-bash}

It will produce an error as the queue channel `ch` is consumed in by the first `.view` operator.

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

We will see in the operator episodes how to create multiple channels from the data.

### Value channels

The second type of Nextflow channel is a `value` channel. A **value** channel by is bound to a *single* value. A value channel can read unlimited times without consuming its content. This can be useful when setting a parameter value which is used by multiple processes.

The Nextflow script below:

~~~
ch = Channel.value('GRCh38')
ch.view()
ch.view()
ch.view()
~~~
{: .source}

Will print:

~~~
GRCh38
GRCh38
GRCh38
~~~
{: .output}


As it is a value channel can be used unlimited times.

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

Channel factories are a way Nextflow can create  different channel types (queue and value).
A factory is the term used when creating an object of a different type e.g. value or queue channels.

### value

The `value` factory method is used to create a value channel.
An optional  argument , values between the `()`, can be specified to bind the channel to a specific value. For example:

~~~
ch1 = Channel.value()
ch2 = Channel.value('GRCh38')
ch2 = Channel.value( ['chr1','chr2','chr3','chr4','chr5'] )
~~~
{: .source}

1. Creates an empty value channel.
1. Creates a value channel and binds a string to it.
1. Creates a value channel and binds a list object to it that will be emitted as a single item.

### of

The factory `Channel.of` allows the creation of a `queue` channel with the values specified as argument.


~~~
ch = Channel.of( 'chr1','chr3','chr5','chr7' )
ch.view()
~~~
{: .source}

The first line in this example creates a variable `ch` which holds a queue channel object.
This channel contains the four values specified as a parameters in the `of` method. Therefore the second line will print the following:

~~~
chr1
chr3
chr5
chr7
~~~
{: .output}

You can specify a range of number, as a single argument,  using the range operator `...`. More information [here](https://www.logicbig.com/tutorials/misc/groovy/range-operator.html).

~~~
ch= Channel.of(1..21, 'X', 'Y')
ch.view()
~~~
{: .source}

> ## Channel.from
> You may see the method `Channel.from` in older nextflow scripts, this performs a similar function but will is deprecated so you should use `Channel.of` instead.
{: .callout}


### fromList

The method `Channel.fromList` creates a channel emitting the elements in a list objects specified as argument. A List object can be defined by placing the values in square brackets `[]` separated by a comma.

~~~
aligner_list = ['salmon', 'kallisto']

aligner_ch = Channel.fromList(aligner_list)

aligner_ch.view()
~~~
{: .source}

This would produce two lines.

~~~
salmon
kallisto
~~~
{: .source}

> ## Channel.fromList vs Channel.of
> In the above example the channel has two elements. If you has used the Channel.of(list) it would have  contained only 1 element `[salmon, kallisto]` and any operator or process using the channel would run once .
{: .callout}

### fromPath

The previous channel factories deal with sending values. A special  channel factory `fromPath` is used when dealing with files.

The `fromPath` factory method create a **queue channel** emitting one or more files matching the specified  pattern specifying the location of files. This pattern can be the location of a single file or a glob pattern that matches multiple files or directories. The location of the file should be specified as a relative path to the location of nextflow script you are running.

The script below creates a channel with a single file as its' content.

~~~
ch = Channel.fromPath( 'data/yeast/reads/ref1_2.fq.gz' )
ch.view()
~~~
{: .source}

~~~
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/ref1_2.fq.gz
~~~
{: .output}

Whilst the script below creates a channel that contains as many items as there are files with `.fq.gz` extension in the `/data/yeast/reads` folder.

~~~
ch = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )
ch.view()
~~~
{: .source}

~~~
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/ref1_2.fq.gz
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/etoh60_3_2.fq.gz
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/temp33_1_2.fq.gz
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/temp33_2_1.fq.gz
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/ref2_1.fq.gz
/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/temp33_3_1.fq.gz
[..truncated..]
~~~
{: .output}

> ## asterisks
> Two asterisks, i.e. `**`, works like * but crosses directory boundaries. This syntax is generally used for matching complete paths. Curly brackets specify a collection of sub-patterns. Learn more about the glob patterns syntax at this [link](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob).
{: .callout}

You can change the behaviour of `Channel.fromPath` method by changing its options. A list of `.fromPath` options is shown below.

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


We can change the default options to give an error if the file doesn't exist using the `checkIfExists` option. In Nextflow option  are separated by a `,` and defined with a colon `:`.

If we execute a nextflow script with the contents below . It will run and not produce an output.
This likely not what we want.
~~~
ch = Channel.fromPath( 'data/chicken/reads/*.fq.gz' )
ch.view()
~~~
{: .source}

It we add the the option `checkIfExists: true `

~~~
ch = Channel.fromPath( 'data/chicken/reads/*.fq.gz', checkIfExists: true )
ch.view()
~~~
{: .source}

This will give an error as there is no data/chicken directory.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [intergalactic_mcclintock] - revision: d2c138894b
No files match pattern `*.fq` at path: /chicken/ggal/
~~~
{: .output}

> ## Using Channel.fromPath
>
>  Use the `Channel.fromPath` method to create a channel emitting all files in the `data/yeast/` then print the file name using the view operator.
>
> > ## Solution
> > You need Two asterisks, i.e. `**`, to  crosses directory boundaries
> > ~~~
> > ch = Channel.fromPath('data/yeast/**')
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
filepair_ch = Channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
filepair_ch.view()
~~~
{: .source}

This will produce a queue channel containing three elements, a value (the glob pattern match), and the set of file pairs. A named list is called a tuple.

The asterisk, `*`, matches any number of characters (including none), and the `{}` braces specify a collection of subpatterns. Therefore the `*_{1,2}` pattern matches any file name ending in `_1.fq.gz` or `_2.fq.gz` .

~~~
[etoh60_3, [/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/etoh60_3_1.fq.gz, /Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/etoh60_3_2.fq.gz]]
[temp33_1, [/Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/temp33_1_1.fq.gz, /Users/ggrimes2/Documents/nextflow-training/data/yeast/reads/temp33_1_2.fq.gz]]
~~~
{: .output}


See more information about the channel factory  `fromFilePairs` [here](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)

> ## More complex patterns
> If you need to match more complex patterns you should create a sample sheet specifying the files and create a channel from that.
{: .callout}

> ## glob pattern
> The glob pattern must contain at least a star wildcard character.
{: .callout}


> ## fromFilePairs
>
>  Use the `fromFilePairs` method to create a channel emitting all pairs of fastq read in the `data/yeast/reads` directory and print them.
>  Then use the `flat:true` option and compare the output with the previous execution.
>
> > ## Solution
> >
> > ~~~
> > Channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz').view()
> > Channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz', flat:true).view()
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
