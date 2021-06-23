---
title: "Channels"
teaching: 30
exercises: 10
questions:
- "How do I get data into Nextflow?"
- "How do I handle different types of input, e.g. files and parameters?"
- "How do I create a Nextflow Channel?"
- "How can I use pattern matching to select input files?"
- "How do I change the way inputs are handled?"
objectives:
- "Understand how Nextflow manages data using Channels."
- "Understand why channels are useful."
- "Understand the different types of Nextflow Channels."
- "Create a value and queue channel using Channel factory methods."
- "Select files as input based on a glob pattern."
- "Edit Channel factory arguments to alter how data is read in."
keypoints:
- "Channels must be used to import data into Nextflow."
- "Nextflow has two different kinds of channels, queue channels and value channels."
- "Data in value channels can be used multiple times in workflow."
- "Data in queue channels are consumed when they are used by a process or an operator."
- "Channel factory methods, such as `Channel.of`, are used to create Channels."
- "Channel factory methods have optional parameters e.g., `checkIfExists`, that can be used to alter the creation and behaviour of a channel."
---

TODO: Describe data used in this episode.

# Channels

Earlier we learnt that channels are the way in which Nextflow sends data around a workflow. Channels connect processes via their inputs and outputs. Channels can store multiple items, such as files (e.g., fastq files) or values. The number of items a channel stores determines how many times a process runs.

## Why use Channels?

Channels let Nextflow handle file management, allowing complex tasks to be split up, run in parallel, and reduces 'admin' required to get the right inputs to the right parts of the pipeline.

![Channel files](../fig/channel-files.png)

Channels are asynchronous, which means that output data from a set of processes will not necessarily be output in the same order as they went in.
However, the first element into a queue is the first out of the queue (First in- First out). This allows processes to run as soon as they receive input from a channel. Channels only send data in onedirection, from a producer (a process/operator), to a consumer (another process/operator).

## Channel types

Nextflow distinguishes between two different kinds of channels: **queue** channels and **value** channels.

### Queue channel

Queue channels are a type of channel in which data is consumed (used up) to make input for a process/operator. Queue channels can be created in two ways,

1. As the outputs of a process.
1. Or a queue channels can be explicitly created using channel factory methods such as [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).


> ## DSL1
> In Nextflow DSL1 queue channels can only be used once in a workflow, either connecting workflow input to process input, or process output to input for another process.
> In DSL2 we can uses a queue channel mutltiple times.
{: .challenge}


### Value channels

The second type of Nextflow channel is a `value` channel. A **value** channel is bound to a *single* value. A value channel can be used an unlimited number times since its content is not consumed. This is also useful for processes that need to reuse input from a channel.


> ## Queue vs Value Channel
>
> What is the difference between a queue and value channel?
> {: .language-groovy }
>
> > ## Solution
> > 1. A queue channels is consumable and can store multiple values. fixme add better description.
> > 1. A value channel a.k.a. singleton channel by definition is bound to a single value
> {: .solution}
{: .challenge}

## Creating Channels using Channel factories

Channel factories are used to explicitly create channels. In programming,
factory methods (functions) are a programming design pattern used
to create different types of objects (in this case, different types
of channels). They are implemented for things that represent more
generalised concepts, such as a `Channel`. Channel factories are
called using the `Channel.<method>` syntax, and return a specific instance
of a `Channel`.

### The value Channel factory

The `value` factory method is used to create a value channel.
Values are put inside  parentheses `()`  to assign them to a channel.

For example:

~~~
ch1 = Channel.value()
ch2 = Channel.value('GRCh38')
ch2 = Channel.value( ['chr1','chr2','chr3','chr4','chr5'] )
~~~
{: .language-groovy }

1. Creates an empty value channel.
1. Creates a value channel and binds a string to it.
1. Creates a value channel and binds a list object to it that will be emitted as a single item.

The value method can only take 1 argument, however, this can be a single list containing several elements.

A [List object](https://www.tutorialspoint.com/groovy/groovy_lists.htm) can be defined by placing the values in square brackets `[]` separated by a comma.

~~~
myList = [1776, -1, 33, 99, 0, 928734928763]
~~~
{: .language-groovy }

## Queue channel factory

Queue (consumable)  channels can be created using the following channel factory methods.

* of
* fromList
* fromPath
* fromFilePairs
* fromSRA

### The of Channel factory

When you want to create a channel containing multiple values you can use the channel factory `Channel.of`. `Channel.of` allows the creation of a `queue` channel with the values specified as arguments, separated by a `,`.

~~~
chromosome_ch = Channel.of( 'chr1','chr3','chr5','chr7' )
chromosome_ch.view()
~~~
{: .language-groovy }

~~~
chr1
chr3
chr5
chr7
~~~
{: .output}

The first line in this example creates a variable `chromosome_ch`. `chromosome_ch` is a queue channel containing the four values specified as arguments in the `of` method. The `view` operator will print one line per item in a list. Therefore the `view` operator on the second line will print four lines, one for each element in the channel:


You can specify a range of numbers as a single argument using the Groovy range operator `..`. This creates each value in the range (including the start and end values) as a value in the channel. The Groovy range operator can also produce ranges of dates, letters, or time.
More information on the range operator can be found [here](https://www.logicbig.com/tutorials/misc/groovy/range-operator.html).

~~~
chromosome_ch= Channel.of(1..22, 'X', 'Y')
chromosome_ch.view()
~~~
{: .language-groovy }

Arguments passed to the `of` method can be of varying types e.g., combinations of numbers, strings, or objects. In the above examples we have examples of both string and number data types.

> ## Channel.from
> You may see the method `Channel.from` in older nextflow scripts. This performs a similar function but is now deprecated (no longer used), and so `Channel.of` should be used instead.
{: .callout}

> ## Create and view Channel contents
> Create a file called `channel.nf` and type the following code into it.
> ~~~
> ch = Channel.of(1,2,3)
> ch.view()
> ~~~
> {: .language-groovy }
> Run the code using
> ~~~~
> nextflow run channel.nf
> ~~~~
> {: .language-groovy }
> How many lines of output do you get?
> > ## Solution
> > In this example you have created a queue channel with three values 1,2,3 in it.
> > This will produce three lines of output, one for each value.
> > The `.view` operator can be used to view the contents of the channel object `ch`.
> > ~~~
> > 1
> > 2
> > 3
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}


### The fromList Channel factory

You can use the `Channel.fromList` method to create a queue channel from a list object.

~~~
aligner_list = ['salmon', 'kallisto']

aligner_ch = Channel.fromList(aligner_list)

aligner_ch.view()
~~~
{: .language-groovy }

This would produce two lines.

~~~
salmon
kallisto
~~~
{: .output}

> ## Channel.fromList vs Channel.of
> In the above example, the channel has two elements. If you has used the Channel.of(list) it would have  contained only 1 element `[salmon, kallisto]` and any operator or process using the channel would run once.
{: .callout}


> ## Creating channels from a list
>
>  Write a Nextflow script that creates both a `queue` and `value` channel
>  for the list
> ~~~
> ids = ['ERR908507', 'ERR908506', 'ERR908505']`
> ~~~
> {: .language-groovy }
>  Then print the contents of the channels using the `view` operator.
>   How many lines does the queue and value channel print ?
>
>  **Hint:** Use the `fromList()` and `value()` Channel factory methods.
> > ## Solution
> >
> > ~~~
> > ids = ['ERR908507', 'ERR908506', 'ERR908505']
> >
> > queue_ch = Channel.fromList(ids)
> > value_ch = Channel.value(ids)
> > queue_ch.view()
> > value_ch.view()
> > ~~~
> > {: .language-groovy }
> > The queue channel `queue_ch` will print three lines.
> >
> > The value channel `value_ch` will print one line.
> {: .solution}
{: .challenge}

### The fromPath Channel factory

The previous channel factory methods dealt with sending general values
in a channel. A special channel factory method `fromPath` is used when wanting to pass files.

The `fromPath` factory method creates a **queue channel** emitting one or more files matching a file path.
The file path (written as a quoted string) can be the location of a single file or a "glob pattern" that matches multiple files or directories. The file path can be a relative path ( path to the file from the current directory), or an absolute path ( path to the file from the system root directory - starts with `/`).

The script below creates a queue channel with a single file as its content.

~~~
read_ch = Channel.fromPath( 'data/yeast/reads/ref1_2.fq.gz' )
read_ch.view()
~~~
{: .language-groovy }

~~~
data/yeast/reads/ref1_2.fq.gz
~~~
{: .output}

The script below creates a queue channel that contains as many items as there are files with `_1.fq.gz` or `_2.fq.gz` extension in the `data/yeast/reads` folder.

~~~
read_ch = Channel.fromPath( 'data/yeast/reads/*_{1,2}.fq.gz' )
read_ch.view()
~~~
{: .language-groovy }

~~~
data/yeast/reads/ref1_2.fq.gz
data/yeast/reads/etoh60_3_2.fq.gz
data/yeast/reads/temp33_1_2.fq.gz
data/yeast/reads/temp33_2_1.fq.gz
data/yeast/reads/ref2_1.fq.gz
data/yeast/reads/temp33_3_1.fq.gz
[..truncated..]
~~~
{: .output}

Two asterisks in the file path, i.e. `**`, works like `*` but will also search sub directories. This syntax is generally used for matching complete paths. Curly brackets `{}` specify a collection of sub-patterns. Learn more about the glob patterns syntax at this [link](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob).


You can change the behaviour of `Channel.fromPath` method by changing its options. A list of `.fromPath` options is shown below.

Available fromPath options:

|Name|Description|
|-----|----------|
| glob | When true, the characters `*`, `?`, `[]` and `{}` are interpreted as glob wildcards, otherwise they are treated as literal characters (default: true)|
| type | The type of file paths matched by the string, either `file`, `dir` or `any` (default: file) |
| hidden | When true, hidden files are included in the resulting paths (default: false)|
| maxDepth | Maximum number of directory levels to visit (default: no limit) |
| followLinks | When true, symbolic links are followed during directory tree traversal, otherwise they are managed as files (default: true) |
| relative | When true returned paths are relative to the top-most common directory (default: false) |
| checkIfExists | When true throws an exception if the specified path does not exist in the file system (default: false)|


We can change the default options for the `fromPath` method to give an error if the file doesn't exist using the `checkIfExists` parameter. In Nextflow, method parameters are separated by a `,` and parameter values specified with a colon `:`.

If we execute a Nextflow script with the contents below . It will run and not produce an output. This is likely not what we want.

~~~
read_ch = Channel.fromPath( 'data/chicken/reads/*.fq.gz' )
read_ch.view()
~~~
{: .language-groovy }

It we add the argument `checkIfExists` with the value `true`.

~~~
read_ch = Channel.fromPath( 'data/chicken/reads/*.fq.gz', checkIfExists: true )
read_ch.view()
~~~
{: .output}

This will give an error as there is no data/chicken directory.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [intergalactic_mcclintock] - revision: d2c138894b
No files match pattern `*.fq` at path: data/chicken/reads/
~~~
{: .output}

> ## Using Channel.fromPath
>
> Use the `Channel.fromPath` method to create a channel containing all files in the `data/yeast/` directory, including the subdirectories.
> Add the parameter `hidden` to include hidden files also.
> Then print all file names using the `view` operator.
>
> **Hint:** You need two asterisks, i.e. `**`, to search subdirectories.
> > ## Solution
> > ~~~
> > all_files_ch = Channel.fromPath('data/yeast/**', hidden: true)
> > all_files_ch.view()
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}


### The fromFilePairs Channel factory

We have seen how to process files individually using `fromPath`. In Bioinformatics we often want to process files in pairs or larger groups, such as read pairs in sequencing.

Nextflow provides a convenient factory method for this common bioinformatics use case. The `fromFilePairs` method creates a queue channel containing a `tuple` for every set of files matching a specific pattern (e.g., `/path/to/*_{1,2}.fq.gz`). A `tuple` is a grouping of data, represented as a Groovy List.
The first element of the tuple emitted from `fromFilePairs` is a string based on the shared part of the filenames (i.e., the `*` part of the glob pattern). The second element is the list of files matching the remaining part of the glob pattern (i.e., the `<string>_{1,2}.fq.gz` pattern).

~~~
read_pair_ch = Channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
read_pair_ch.view()
~~~
{: .language-groovy }

~~~
[etoh60_3, [data/yeast/reads/etoh60_3_1.fq.gz, data/yeast/reads/etoh60_3_2.fq.gz]]
[temp33_1, [data/yeast/reads/temp33_1_1.fq.gz, data/yeast/reads/temp33_1_2.fq.gz]]
[ref1, [data/yeast/reads/ref1_1.fq.gz, data/yeast/reads/ref1_2.fq.gz]]
[ref2, [data/yeast/reads/ref2_1.fq.gz, data/yeast/reads/ref2_2.fq.gz]]
[temp33_2, [data/yeast/reads/temp33_2_1.fq.gz, data/yeast/reads/temp33_2_2.fq.gz]]
[ref3, [data/yeast/reads/ref3_1.fq.gz, data/yeast/reads/ref3_2.fq.gz]]
[temp33_3, [data/yeast/reads/temp33_3_1.fq.gz, data/yeast/reads/temp33_3_2.fq.gz]]
[etoh60_1, [data/yeast/reads/etoh60_1_1.fq.gz, data/yeast/reads/etoh60_1_2.fq.gz]]
[etoh60_2, [data/yeast/reads/etoh60_2_1.fq.gz, data/yeast/reads/etoh60_2_2.fq.gz]]
~~~
{: .output}

This will produce a queue channel, `read_pair_ch` , containing nine elements. Each element is a tuple that has a string value (the file prefix matched, e.g `temp33_1`) and a list with the two files e,g. `[data/yeast/reads/etoh60_2_1.fq.gz, data/yeast/reads/etoh60_2_2.fq.gz]` .

The asterisk character `*`, matches any number of characters (including none), and the `{}` braces specify a collection of subpatterns. Therefore the `*_{1,2}.fq.gz` pattern matches any file name ending in `_1.fq.gz` or `_2.fq.gz` .


#### What if you want to capture more than a pair?

If you want to capture more than two files for a pattern you will need to change the default `size` argument to the number of expected matching files.

For Example:

~~~
read_group_ch = Channel.fromFilePairs('data/yeast/reads/ref{1,2,3}*',size:6)
read_group_ch.view()
~~~
{: .language-groovy }

The code above will create a queue channel containing one element. The element is a tuple of which contains a string value, that is the pattern **ref**, and a list of six files matching the pattern.

~~~
[ref, [data/yeast/reads/ref1_1.fq.gz, data/yeast/reads/ref1_2.fq.gz, data/yeast/reads/ref2_1.fq.gz, data/yeast/reads/ref2_2.fq.gz, data/yeast/reads/ref3_1.fq.gz, data/yeast/reads/ref3_2.fq.gz]]
~~~
{: .output}

See more information about the channel factory `fromFilePairs` [here](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)

> ## More complex patterns
> If you need to match more complex patterns you should create a sample sheet specifying the files and create a channel from that. This will be covered in the operator episode.
{: .callout}

> ## The glob pattern
> The pattern must contain at least a star wildcard character.
{: .callout}


> ## fromFilePairs
>
>  1. Use the `fromFilePairs` method to create a channel containing all pairs of fastq read in the `data/yeast/reads` directory and print them.
>  1. Then use the `fromFilePairs` `size` argument with parameter value 6 and the pattern `data/yeast/reads/*_{1,2,3}_{1,2}.fq.gz`.
> 2. Compare the output with the previous execution.
>
> > ## Solution
> >
> > ~~~
> > pairs_ch = Channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
> > pairs_ch.view()
> > multi_ch = Channel.fromFilePairs('data/yeast/reads/*_{1,2,3}_{1,2}.fq.gz', size:6)
> > multi_ch.view()
> > ~~~
> > {: .language-groovy }
> > The first method will create a channel with 9 elements.
> > The second method with create a channel with 2 elements. The pattern will match files starting temp33 and etoh60.
> {: .solution}
{: .challenge}

### The fromSRA Channel factory

Another useful factory method is `fromSRA`. The `fromSRA` method makes it possible to query the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) archive and returns a queue channel emitting the FASTQ files matching the specified selection criteria.

The queries can be project IDs or accession numbers supported by the [NCBI ESearch API](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch).

~~~
sra_ch =Channel.fromSRA('SRP043510')
sra_ch.view()
~~~
{: .language-groovy }

This will print a tuple for every fastq file associated with that SRA project accession.

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
sra_ch = Channel.fromSRA(ids)
sra_ch.view()
~~~
{: .language-groovy }

~~~
[ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
~~~
{: .output}  

> ## Read pairs
> Read pairs are implicitly managed, and are returned as a list of files.
{: .callout}


{% include links.md %}
