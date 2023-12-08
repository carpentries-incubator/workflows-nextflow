---
title: Channels
teaching: 30
exercises: 10
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand how Nextflow manages data using channels.
- Understand the different types of Nextflow channels.
- Create a value and queue channel using channel factory methods.
- Select files as input based on a glob pattern.
- Edit channel factory arguments to alter how data is read in.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I get data into Nextflow?
- How do I handle different types of input, e.g. files and parameters?
- How do I create a Nextflow channel?
- How can I use pattern matching to select input files?
- How do I change the way inputs are handled?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Channels

Earlier we learnt that channels are the way in which Nextflow sends data around a workflow. Channels connect processes via their inputs and outputs. Channels can store multiple items, such as files (e.g., fastq files) or values. The number of items a channel stores determines how many times a process will run using that channel as input.  
*Note:* When the process runs using one item from the input channel, we will call that run a `task`.

## Why use Channels?

Channels are how Nextflow handles file management, allowing complex tasks to be split up, run in parallel, and reduces 'admin' required to get the right inputs to the right parts of the pipeline.

![](fig/channel-files.png){alt='Channel files'}

Channels are asynchronous, which means that outputs from a set of processes will not necessarily be produced in the same order as the corresponding inputs went in.
However, the first element into a channel queue is the first out of the queue (First in - First out). This allows processes to run as soon as they receive input from a channel. Channels only send data in one direction, from a producer (a process/operator), to a consumer (another process/operator).

## Channel types

Nextflow distinguishes between two different kinds of channels: **queue** channels and **value** channels.

### Queue channel

Queue channels are a type of channel in which data is consumed (used up) to make input for a process/operator. Queue channels can be created in two ways:

1. As the outputs of a process.
2. Explicitly using channel factory methods such as [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).

:::::::::::::::::::::::::::::::::::::::::  callout

## DSL1

In Nextflow DSL1 queue channels can only be used once in a workflow, either connecting workflow input to process input, or process output to input for another process.
In DSL2 we can use a queue channel multiple times.


::::::::::::::::::::::::::::::::::::::::::::::::::

### Value channels

The second type of Nextflow channel is a `value` channel. A **value** channel is bound to a **single** value. A value channel can be used an unlimited number times since its content is not consumed. This is also useful for processes that need to reuse input from a channel, for example, a reference genome sequence file that is required by multiple steps within a process, or by more than one process.

:::::::::::::::::::::::::::::::::::::::  challenge

## Queue vs Value Channel.

What type of channel would you use to store the following?

1. Multiple values.
2. A list with one or more values.
3. A single value.

:::::::::::::::  solution

## Solution

1. A queue channels is used to store multiple values.
2. A value channel is used to store a single value, this can be a list with multiple values.
3. A value channel is used to store a single value.
  
  

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Creating Channels using Channel factories

Channel factories are used to explicitly create channels. In programming,
factory methods (functions) are a programming design pattern used
to create different types of objects (in this case, different types
of channels). They are implemented for things that represent more
generalised concepts, such as a `Channel`.

Channel factories are
called using the `Channel.<method>` syntax, and return a specific instance
of a `Channel`.

### The value Channel factory

The `value` factory method is used to create a value channel.
Values are put inside  parentheses `()`  to assign them to a channel.

For example:

```groovy 
ch1 = Channel.value( 'GRCh38' )
ch2 = Channel.value( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] )
ch3 = Channel.value( ['chr1' : 248956422, 'chr2' : 242193529, 'chr3' : 198295559] )
```

1. Creates a value channel and binds a string to it.
2. Creates a value channel and binds a list object to it that will be emitted as a single item.
3. Creates a value channel and binds a map object to it that will be emitted as a single item.

The value method can only take 1 argument, however, this can be a single list or map containing several elements.

*Reminder:*

- A [List object](https://www.tutorialspoint.com/groovy/groovy_lists.htm) can be defined by placing the values in square brackets `[]` separated by a comma.
- A [Map object](https://www.tutorialspoint.com/groovy/groovy_maps.htm) is similar, but with `key:value pairs` separated by commas.

```groovy 
myList = [1776, -1, 33, 99, 0, 928734928763]
myMap = [ p1 : "start", q2 : "end" ]
```

## Queue channel factory

Queue (consumable) channels can be created using the following channel factory methods.

- Channel.of
- Channel.fromList
- Channel.fromPath
- Channel.fromFilePairs
- Channel.fromSRA

### The **of** Channel factory

When you want to create a channel containing multiple values you can use the channel factory `Channel.of`. `Channel.of` allows the creation of a `queue` channel with the values specified as arguments, separated by a `,`.

```groovy 
chromosome_ch = Channel.of( 'chr1', 'chr3', 'chr5', 'chr7' )
chromosome_ch.view()
```

```output
chr1
chr3
chr5
chr7
```

The first line in this example creates a variable `chromosome_ch`. `chromosome_ch` is a queue channel containing the four values specified as arguments in the `of` method. The `view` operator will print one line per item in a list. Therefore the `view` operator on the second line will print four lines, one for each element in the channel:

You can specify a range of numbers as a single argument using the Groovy range operator `..`. This creates each value in the range (including the start and end values) as a value in the channel. The Groovy range operator can also produce ranges of dates, letters, or time.
More information on the range operator can be found [here](https://www.logicbig.com/tutorials/misc/groovy/range-operator.html).

```groovy 
chromosome_ch = Channel.of(1..22, 'X', 'Y')
chromosome_ch.view()
```

Arguments passed to the `of` method can be of varying types e.g., combinations of numbers, strings, or objects. In the above examples we have examples of both string and number data types.

:::::::::::::::::::::::::::::::::::::::::  callout

## Channel.from

You may see the method `Channel.from` in older nextflow scripts. This performs a similar function but is now deprecated (no longer used), and so `Channel.of` should be used instead.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Create a value and Queue and view Channel contents

1. Create a Nextflow script file called `channel.nf` .
2. Create a Value channel `ch_vl` containing the String `'GRCh38'`.
3. Create a Queue channel `ch_qu` containing the values  1 to 4.
4. Use `.view()` operator on the channel objects to view the contents of the channels.
5. Run the code using

```bash 
$ nextflow run channel.nf
```

:::::::::::::::  solution

## Solution

```groovy 
ch_vl = Channel.value('GRCh38')
ch_qu = Channel.of(1,2,3,4)
ch_vl.view()
ch_qu.view()
```

```output
 N E X T F L O W  ~  version 21.04.0
 Launching `channel.nf` [condescending_dalembert] - revision: c80908867b
 GRCh38
 1
 2
 3
 4
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### The **fromList** Channel factory

You can use the `Channel.fromList` method to create a queue channel from a list object.

```groovy 
aligner_list = ['salmon', 'kallisto']

aligner_ch = Channel.fromList(aligner_list)

aligner_ch.view()
```

This would produce two lines.

```output
salmon
kallisto
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Channel.fromList vs Channel.of

In the above example, the channel has two elements. If you has used the Channel.of(aligner\_list) it would have  contained only 1 element `[salmon, kallisto]` and any operator or process using the channel would run once.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Creating channels from a list

Write a Nextflow script that creates both a `queue` and `value` channel
for the list

```groovy 
ids = ['ERR908507', 'ERR908506', 'ERR908505']
```

Then print the contents of the channels using the `view` operator.
How many lines does the queue and value channel print?

**Hint:** Use the `fromList()` and `value()` Channel factory methods.

:::::::::::::::  solution

## Solution

```groovy 
ids = ['ERR908507', 'ERR908506', 'ERR908505']

queue_ch = Channel.fromList(ids)
value_ch = Channel.value(ids)
queue_ch.view()
value_ch.view()
```

```output 
N E X T F L O W  ~  version 21.04.0
Launching `channel_fromList.nf` [wise_hodgkin] - revision: 22d76be151
ERR908507
ERR908506
ERR908505
[ERR908507, ERR908506, ERR908505]
```

The queue channel `queue_ch` will print three lines.

The value channel `value_ch` will print one line.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### The **fromPath** Channel factory

The previous channel factory methods dealt with sending general values in a channel. A special channel factory method `fromPath` is used when wanting to pass files.

The `fromPath` factory method creates a **queue channel** containing one or more files matching a file path.

The file path (written as a quoted string) can be the location of a single file or a "glob pattern" that matches multiple files or directories.

The file path can be a relative path (path to the file from the current directory), or an absolute path (path to the file from the system root directory - starts with `/`).

The script below creates a queue channel with a single file as its content.

```groovy 
read_ch = Channel.fromPath( 'data/yeast/reads/ref1_2.fq.gz' )
read_ch.view()
```

```output
data/yeast/reads/ref1_2.fq.gz
```

You can also use  glob syntax to specify pattern-matching behaviour for files.
A glob pattern is specified as a string and is matched against directory or file names.

- An asterisk, `*`, matches any number of characters (including none).
- Two asterisks, `**`, works like \*  but will also search sub directories. This syntax is generally used for matching complete paths.
- Braces  `{}` specify a collection of subpatterns. For example:
  `{bam,bai}` matches "bam" or "bai"

For example the script below uses the `*.fq.gz` pattern to create a queue channel that contains as many items as there are files with `.fq.gz` extension in the `data/yeast/reads` folder.

```groovy 
read_ch = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )
read_ch.view()
```

```output
data/yeast/reads/ref1_2.fq.gz
data/yeast/reads/etoh60_3_2.fq.gz
data/yeast/reads/temp33_1_2.fq.gz
data/yeast/reads/temp33_2_1.fq.gz
data/yeast/reads/ref2_1.fq.gz
data/yeast/reads/temp33_3_1.fq.gz
[..truncated..]
```

**Note** The pattern must contain at least a star wildcard character.

You can change the behaviour of `Channel.fromPath` method by changing its options. A list of `.fromPath` options is shown below.

Available fromPath options:

| Name          | Description                                                                                                                 | 
| ------------- | --------------------------------------------------------------------------------------------------------------------------- |
| glob          | When true, the characters `*`, `?`, `[]` and `{}` are interpreted as glob wildcards, otherwise they are treated as literal characters (default: true)                                                                                                  | 
| type          | The type of file paths matched by the string, either `file`, `dir` or `any` (default: file)                                                                       | 
| hidden        | When true, hidden files are included in the resulting paths (default: false)                                                | 
| maxDepth      | Maximum number of directory levels to visit (default: no limit)                                                             | 
| followLinks   | When true, symbolic links are followed during directory tree traversal, otherwise they are managed as files (default: true) | 
| relative      | When true returned paths are relative to the top-most common directory (default: false)                                     | 
| checkIfExists | When true throws an exception if the specified path does not exist in the file system (default: false)                      | 

We can change the default options for the `fromPath` method to give an error if the file doesn't exist using the `checkIfExists` parameter. In Nextflow, method parameters are separated by a `,` and parameter values specified with a colon `:`.

If we execute a Nextflow script with the contents below, it will run and not produce an output. This is likely not what we want.

```groovy 
read_ch = Channel.fromPath( 'data/chicken/reads/*.fq.gz' )
read_ch.view()
```

Add the argument `checkIfExists` with the value `true`.

```groovy
read_ch = Channel.fromPath( 'data/chicken/reads/*.fq.gz', checkIfExists: true )
read_ch.view()
```

This will give an error as there is no data/chicken directory.

```output
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [intergalactic_mcclintock] - revision: d2c138894b
No files match pattern `*.fq.gz` at path: data/chicken/reads/
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Using Channel.fromPath

1. Create a Nextflow script `channel_fromPath.nf`
2. Use the `Channel.fromPath` method to create a channel containing all files in the `data/yeast/` directory, including the subdirectories.
3. Add the parameter to include any hidden files.
4. Then print all file names using the `view` operator.

**Hint:** You need two asterisks, i.e. `**`, to search subdirectories.

:::::::::::::::  solution

## Solution

```groovy 
all_files_ch = Channel.fromPath('data/yeast/**', hidden: true)
all_files_ch.view()
```

```output 
N E X T F L O W  ~  version 21.04.0
Launching `nf-training/scripts/channels/channel_fromPath.nf` [reverent_mclean] - revision: cf02269bcb
data/yeast/samples.csv
data/yeast/bams/ref1.bam.bai
data/yeast/bams/ref3.bam
data/yeast/bams/etoh60_3.
[..truncated..]
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### The **fromFilePairs** Channel factory

We have seen how to process files individually using `fromPath`. In Bioinformatics we often want to process files in pairs or larger groups, such as read pairs in sequencing.

For example is the `data/yeast/reads` directory we have nine groups of read pairs.

| Sample group  | read1                                                                                                                       | read2                             | 
| ------------- | --------------------------------------------------------------------------------------------------------------------------- | --------------------------------- |
| ref1          | data/yeast/reads/ref1\_1.fq.gz                                                                                               | data/yeast/reads/ref1\_2.fq.gz     | 
| ref2          | data/yeast/reads/ref2\_1.fq.gz                                                                                               | data/yeast/reads/ref2\_2.fq.gz     | 
| ref3          | data/yeast/reads/ref3\_1.fq.gz                                                                                               | data/yeast/reads/ref3\_2.fq.gz     | 
| temp33\_1      | data/yeast/reads/temp33\_1\_1.fq.gz                                                                                           | data/yeast/reads/temp33\_1\_2.fq.gz | 
| temp33\_2      | data/yeast/reads/temp33\_2\_1.fq.gz                                                                                           | data/yeast/reads/temp33\_2\_2.fq.gz | 
| temp33\_3      | data/yeast/reads/temp33\_3\_1.fq.gz                                                                                           | data/yeast/reads/temp33\_3\_2.fq.gz | 
| etoh60\_1      | data/yeast/reads/etoh60\_1\_1.fq.gz                                                                                           | data/yeast/reads/etoh60\_1\_2.fq.gz | 
| etoh60\_2      | data/yeast/reads/etoh60\_2\_1.fq.gz                                                                                           | data/yeast/reads/etoh60\_2\_2.fq.gz | 
| etoh60\_3      | data/yeast/reads/etoh60\_3\_1.fq.gz                                                                                           | data/yeast/reads/etoh60\_3\_2.fq.gz | 

Nextflow provides a convenient Channel factory method for this common bioinformatics use case. The `fromFilePairs` method creates a queue channel containing a `tuple` for every set of files matching a specific glob pattern (e.g., `/path/to/*_{1,2}.fq.gz`).

A `tuple` is a grouping of data, represented as a Groovy List.

1. The first element of the tuple emitted from `fromFilePairs` is a string based on the shared part of the filenames (i.e., the `*` part of the glob pattern).
2. The second element is the list of files matching the remaining part of the glob pattern (i.e., the `<string>_{1,2}.fq.gz` pattern). This will include any files ending `_1.fq.gz` or `_2.fq.gz`.

```groovy 
read_pair_ch = Channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
read_pair_ch.view()
```

```output
[etoh60_3, [data/yeast/reads/etoh60_3_1.fq.gz, data/yeast/reads/etoh60_3_2.fq.gz]]
[temp33_1, [data/yeast/reads/temp33_1_1.fq.gz, data/yeast/reads/temp33_1_2.fq.gz]]
[ref1, [data/yeast/reads/ref1_1.fq.gz, data/yeast/reads/ref1_2.fq.gz]]
[ref2, [data/yeast/reads/ref2_1.fq.gz, data/yeast/reads/ref2_2.fq.gz]]
[temp33_2, [data/yeast/reads/temp33_2_1.fq.gz, data/yeast/reads/temp33_2_2.fq.gz]]
[ref3, [data/yeast/reads/ref3_1.fq.gz, data/yeast/reads/ref3_2.fq.gz]]
[temp33_3, [data/yeast/reads/temp33_3_1.fq.gz, data/yeast/reads/temp33_3_2.fq.gz]]
[etoh60_1, [data/yeast/reads/etoh60_1_1.fq.gz, data/yeast/reads/etoh60_1_2.fq.gz]]
[etoh60_2, [data/yeast/reads/etoh60_2_1.fq.gz, data/yeast/reads/etoh60_2_2.fq.gz]]
```

This will produce a queue channel, `read_pair_ch` , containing nine elements.

Each element is a tuple that has;

1. string value (the file prefix matched, e.g `temp33_1`)
2. and a list with the two files e,g. `[data/yeast/reads/temp33_1_1.fq.gz, data/yeast/reads/temp33_1_2.fq.gz]` .

The asterisk character `*`, matches any number of characters (including none), and the `{}` braces specify a collection of subpatterns. Therefore the `*_{1,2}.fq.gz` pattern matches any file name ending in `_1.fq.gz` or `_2.fq.gz` .

### What if you want to capture more than a pair?

If you want to capture more than two files for a pattern you will need to change the default `size` argument (the default value is 2) to the number of expected matching files.

For example in the directory `data/yeast/reads` there are six files with the prefix `ref`.
If we want to group (create a tuple) for all of these files we could write;

```groovy 
read_group_ch = Channel.fromFilePairs('data/yeast/reads/ref{1,2,3}*',size:6)
read_group_ch.view()
```

The code above will create a queue channel containing one element. The element is a tuple of which contains a string value, that is the pattern **ref**, and a list of six files matching the pattern.

```output
[ref, [data/yeast/reads/ref1_1.fq.gz, data/yeast/reads/ref1_2.fq.gz, data/yeast/reads/ref2_1.fq.gz, data/yeast/reads/ref2_2.fq.gz, data/yeast/reads/ref3_1.fq.gz, data/yeast/reads/ref3_2.fq.gz]]
```

See more information about the channel factory `fromFilePairs` [here](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)

:::::::::::::::::::::::::::::::::::::::::  callout

## More complex patterns

If you need to match more complex patterns you should create a sample sheet specifying the files and create a channel from that. This will be covered in the operator episode.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Create a channel containing groups of files

1. Create a Nextflow script file `channel_fromFilePairs.nf` .
2. Use the `fromFilePairs` method to create a channel containing three tuples. Each tuple will contain the pairs of fastq reads for the three temp33 samples in the `data/yeast/reads` directory

:::::::::::::::  solution

## Solution

```groovy 
pairs_ch = Channel.fromFilePairs('data/yeast/reads/temp33*_{1,2}.fq.gz')
pairs_ch.view()
```

```output 
N E X T F L O W  ~  version 21.04.0
Launching `channels.nf` [stupefied_lumiere] - revision: a3741edde2
[temp33_1, [data/yeast/reads/temp33_1_1.fq.gz, data/yeast/reads/temp33_1_2.fq.gz]]
[temp33_3, [data/yeast/reads/temp33_3_1.fq.gz, data/yeast/reads/temp33_3_2.fq.gz]]
[temp33_2, [data/yeast/reads/temp33_2_1.fq.gz, data/yeast/reads/temp33_2_2.fq.gz]]
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### The **fromSRA** Channel factory

Another useful factory method is `fromSRA`. The `fromSRA` method makes it possible to query the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) archive and returns a queue channel emitting the FASTQ files matching the specified selection criteria.

The queries can be project IDs or accession numbers supported by the [NCBI ESearch API](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch).

If you want to use this functionality, you will need an [NCBI API KEY](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/), and to set the environment variable `NCBI_API_KEY` to its value.

```groovy 
sra_ch =Channel.fromSRA('SRP043510')
sra_ch.view()
```

This will print a tuple for every fastq file associated with that SRA project accession.

```output
[SRR1448794, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448794/SRR1448794.fastq.gz]
[SRR1448795, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/005/SRR1448795/SRR1448795.fastq.gz]
[SRR1448792, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/002/SRR1448792/SRR1448792.fastq.gz]
[SRR1448793, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/003/SRR1448793/SRR1448793.fastq.gz]
[SRR1910483, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/003/SRR1910483/SRR1910483.fastq.gz]
[SRR1910482, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/002/SRR1910482/SRR1910482.fastq.gz]
(remaining omitted)
```

Multiple accession IDs can be specified using a list object:

```groovy 
ids = ['ERR908507', 'ERR908506', 'ERR908505']
sra_ch = Channel.fromSRA(ids)
sra_ch.view()
```

```output
[ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Read pairs from SRA

Read pairs are implicitly managed, and are returned as a list of files.


::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::: keypoints

- Channels must be used to import data into Nextflow.
- Nextflow has two different kinds of channels: queue channels and value channels.
- Data in value channels can be used multiple times in workflow.
- Data in queue channels are consumed when they are used by a process or an operator.
- Channel factory methods, such as `Channel.of`, are used to create channels.
- Channel factory methods have optional parameters e.g., `checkIfExists`, that can be used to alter the creation and behaviour of a channel.

::::::::::::::::::::::::::::::::::::::::::::::::::


