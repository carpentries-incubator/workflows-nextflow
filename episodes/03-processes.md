---
title: "Processes"
teaching: 50 min
exercises: 20 min
questions:
- "How do I run tasks/processes in Nextflow?"
- "How do I get data, files and values, into and out of processes?"
- "How do can I control when a process is executed?"
- "How do I control resources, such as number of CPUs and memory, available to processes?"
- "How do I save output/results from a process?"
objectives:
- "Understand how Nextflow uses processes to execute tasks."
- "Create a Nextflow process."
- "Define inputs and outputs to a process."
- "Understand how to use conditionals control process execution."
- "Use process directives to control execution of a process."
- "Use the publishDir directive to save result files to a directory."
keypoints:
- "A Nextflow process is an independent task/step in a workflow"
- "Processes contain up to five definition blocks including, directives, inputs, outputs, when clause and finally a script block."
- "The script block contains the commands you would like to run."
- "Inputs and Outputs to a process are defined using the input and output blocks."
- "The execution of a process can be controlled using conditional statements."
- "Files produced within a process can be saved to a directory using the `publishDir` directive."
---


# Processes

We now know how to create and use Channels to send data around a workflow. We will now see how to run tasks within a workflow.

A `process` is the way Nextflow execute commands you would run on the command line or custom scripts.


Processes can be thought of as a particular task/steps in a workflow, e.g. an alignment step in RNA-Seq analysis. Processes  are independent of each other (don't require another processes to execute) and can not communicate/write to each other . It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output channels.

The process definition starts with keyword the `process`, followed by process name and finally the process `body` delimited by curly brackets `{}`. The process body must contain a string  which represents the command or, more generally, a script that is executed by it.

A basic process with no input or output channels looks like the following example:

~~~
process index {

  script:
  """
  salmon index -t ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmerLen 31
  """
}
~~~
{: .language-groovy }

This example build a salmon index for the yeast transcriptome. We use the nextflow variable `${projectDir}` to specify the directory where the main script is located .


### definition blocks

A process may contain five definition blocks, that control how a command is executed within a process, these  are:

1. **directives**: allow the definition of optional settings that affect the execution of the current process e.g. the number of cpus a task uses and the amount of memory allocated.
1. **inputs**: Define the input dependencies, usually channels, which determines the number of times a process is executed.
1. **outputs**: Defines the output channels used by the process to send results/data produced by the process.
1. **when clause**: Allow you to define a condition that must be verified in order to execute the process.
1. **The script block**: The script block is a statement within quotes that defines the command that is executed by the process to carry out its task.


The syntax is defined as follows:

~~~
process < name > {
  [ directives ]        
  input:                
  < process inputs >
  output:               
  < process outputs >
  when:                 
  < condition >
  [script|shell|exec]:  
  < user script to be executed >
}
~~~
{: .language-groovy }

* Zero, one or more process directives, e.g cpus
* Zero, one or more process inputs
* Zero, one or more process outputs
* An optional boolean conditional to trigger the process execution
* The command to be executed


## Script

At minimum a process block must contain a `script` block.

The `script` block is a string "statement" that defines the command that is executed by the process to carry out its task. These are normally the commands you would run on a terminal.

A process contains only one script block, and it must be the last statement when the process contains input and output declarations.

The script block can be a simple string in single quotes e.g. `"samtools index ref1.sorted.bam"` .
For commands that span multiple lines you can encase the command in the triple quote.

For example:

~~~
process multi_line {
    script:
    """
    samtools sort -o ref1.sorted.bam ${projectDir}/data/yeast/bams/ref1.bam
    samtools index ref1.sorted.bam
    samtools flagstat ref1.sorted.bam
    """
}
~~~
{: .language-groovy }

By default the process command is interpreted as a **Bash** script. However any other scripting language can be used just simply starting the script with the corresponding [Shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) declaration. For example:

~~~
process pyStuff {
  script:
  """
  #!/usr/bin/env python
  import gzip

  reads = 0
  bases = 0

  with gzip.open('data/yeast/reads/ref1_1.fq.gz', 'rb') as read:
      for id in read:
          seq = next(read)
          reads += 1
          bases += len(seq.strip())
          next(read)
          next(read)

  print("reads", reads)
  print("bases", bases)
  """
}
~~~
{: .language-groovy }

~~~
process rStuff {
  script:
  """
  #!/usr/bin/env Rscript
  library("ShortRead")
  countFastq(dirPath="data/yeast/reads/ref1_1.fq.gz")
  """
}
~~~
{: .language-groovy }

This allows the the use of a different programming languages which may better fit a particular job. However for large chunks of code is suggested to save them into separate files and invoke them from the process script.

~~~
process pyStuff {
  script:
  """
  python myscript.py
  """
}
~~~
{: .language-groovy }

### Script parameters

The command in the script block can be defined dynamically using Nextflow variables e.g. `${projectDir}`.
To reference a variable in the script block you can use the `$` in front of the nextflow variable name, and additionally you can add `{}` around the variable name e.g. `${projectDir}`.


In the example below the Nextflow variable `kmer` is set to the value 31.

~~~
kmer = 31

process index {

  script:
  """
  salmon index \
  -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  -i index \
  --kmer $kmer
  echo "kmer size is $kmer"
  """
}
~~~
{: .language-groovy }

In most cases we do not want to hard code parameter values.
A special Nextflow map variable `params` can be used to assign values from the command line.

In the example below we define the variable `params.kmer` with a default value of 31 in the Nextflow script.
~~~
params.kmer = 31

process index {

  script:
  """
  salmon index \
  -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  -i index \
  --kmer $params.kmer
  echo "kmer size is $params.kmer"
  """
}
~~~
{: .language-groovy }

We can change the default value of `kmer` to 11 by running the Nextflow script using the command below. *Remember* parameter have two hyphens `--` .

~~~
nextflow run script.nf --kmer 11
~~~
{: .bash }


> ## Script parameters
>
> Change the k-mer value used to create the salmon index from the command line.
> ~~~
> nextflow run process_script_params.nf --kmer <some value>
> ~~~
> Note the kmer values must not be greater than 31 and an odd number.
>
> > ## Solution
> > ~~~
> > nextflow run process_script_params.nf --kmer 27
> > ~~~
> {: .solution}
{: .challenge}


### Bash variables

 Nextflow uses the same Bash syntax for variable substitutions, `$variable`, in strings, Bash variables need to be escaped using `\` character.

 In the example below we will use the bash `PWD` variable.


~~~
params.kmer = 31

process index {

  script:
  """
  salmon index -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer $params.kmer
  echo "kmer size is $params.kmer"
  echo "index is located in" \$PWD
  """
}
~~~
{: .language-groovy }


This will print the location of the working directory using the bash environment variable `PWD`.


Another alternative is to use a `shell` statement instead of `script` which uses a different syntax for Nextflow variable: `!{..}`. This allow enables you to use both Nextflow and Bash variables in the same script.

```
params.aligner = 'salmon'

process aligner_log {
  shell:
  '''
  X='Align using'
  echo $X !{params.aligner}
  '''
}
```
{: .language-groovy }


### Conditional execution

Sometimes you want to change how a process is run depending on some parameter. In Nextflow scripts we can use conditional statements such as the `if` statement or any other expression evaluating to string value. For example, the Nextflow script below will change what index is created depending on the Nextflow variable `params.aligner`.

~~~
params.aligner = 'kallisto'
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.kmer = 31

process index {
  script:
  if( params.aligner == 'kallisto' )
    """
    echo index using kallisto
    kallisto index -i index  -k $params.kmer $params.transcriptome
    """
  else if( params.aligner == 'salmon' )
    """
    echo index using salmon
    salmon index -t $params.transcriptome -i index --kmer $params.kmer
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.aligner")
}
~~~
{: .language-groovy }


## Inputs

Processes are isolated from each other but can communicate by sending values and files via Nextflow channels into input and output blocks.

The input block determine the dependency and the number of time a process executes.

![Process Flow](../fig/channel-process.png)


The `input` block defines which channels the process is expecting to receive input from. You can only define one input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below:

~~~
input:
  <input qualifier> <input name> from <source channel>
~~~
{: .language-groovy }

The input qualifier declares the type of data to be received.

> ## input qualifiers
> * `val`: Lets you access the received input value by its name in the process script.
> * `env`: Lets you use the received value to set an environment variable named as > the specified input name.
> * `path`: Lets you handle the received value as a path, staging the file properly in the execution context.
> * `stdin`: Lets you forward the received value to the process stdin special file.
> * `tuple`: Lets you handle a group of input values having one of the above qualifiers.
> * `each`: Lets you execute the process for each entry in the input collection.
{: .callout}

### Input values

The `val` qualifier allows you to receive value data as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

~~~
chr_ch = Channel.of( 1..22,'X','Y' )

process printChr {

  input:
  val chr from chr_ch

  script:
  """
  echo processing chromosome $chr
  """
}
~~~
{: .language-groovy }

In the above example the process is executed 24 times; each time a value is received from the queue channel `chr_ch` it is used to run process. Thus, it results in an output similar to the one shown below:

~~~
processing chromosome 3
processing chromosome 1
processing chromosome 2
..truncated...
~~~
{: .output}

> ## Channel order
> The channel guarantees that items are delivered in the same order as they have been sent,  but,  since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received.
{: .callout}

### Input files

The `path` qualifier allows the handling of files. This means that Nextflow will stage it in the process execution directory, and it can be access in the script by using the name specified in the input declaration.

The input name can be either defined as user specified filename inside quotes as in the example below where specific the name of the file as `sample.fastq`.

~~~
reads = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )

process listFiles {

    input:
    path 'sample.fastq' from reads

    script:
    """
    ls -l sample.fastq
    """
}
~~~
{: .language-groovy }


Or more commonly the input file name can be defined dynamically by defining the input name as a Nextflow variable, `sample`, as shown below:

~~~
reads = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )

process listFiles {
    input:
    path sample from reads
    script:
    """
    ls -l $sample
    """
}
~~~
{: .language-groovy }


> # File Objects as inputs
> When a process declares an input file the corresponding channel elements must be file objects i.e. created with the path helper function from the file specific channel factories e.g. `Channel.fromPath` or `Channel.fromFilePairs`.
{: .callout}


>  Exercise
>  Write a Nextflow script `fastqc.nf` that creates a queue channel, `Channel.fromPath`, containing all read files matching the pattern `data/yeast/reads/*_1.fq.gz` followed by a process with input and script directives that has the commands.
> `mkdir fastqc_out`
> `fastqc -o fastqc_out ${reads}`
> > Solution
> > ~~~
> > reads_ch = Channel.fromPath( 'data/yeast/reads/*_1.fq.gz' )
> >
> > process fastqc {
> >    input:
> >    path fastq from reads_ch
> >    script:
> >    """
> >    mkdir fastqc_out
> >    fastqc -o fastqc_out ${reads}
> >    """
> >}
> >~~~
> {: .solution}
{: .challenge}

### Combining input channels


A key feature of processes is the ability to handle inputs from multiple channels.
However it’s important to understands how the content of channel and affect the execution of a process.

Consider the following example:

~~~
ch_num = Channel.of(1,2,3)
ch_letters = Channel.of('a','b','c')

process combine_channels {
  input:
  val x from ch_num
  val y from ch_letters
  script:
   """
   echo $x and $y
   """
}
~~~
{: .language-groovy }

Both channels contain three value, therefore the process is executed three times, each time with a different pair:

~~~
2 and b

1 and a

3 and c
~~~
{: .output}

What is happening is that the process waits until it receives an input value from all the queue channels declared as input.

When this condition is verified, it uses up the input values coming from the respective queue channels, runs the task. This logic repeats until one or more queue channels have no more content. The process then stops.

What happens when not all channels have the same number of elements?

For example:

~~~

ch_num = Channel.of(1,2)
ch_letters = Channel.of('a','b','c','d')

process combine_channels {
  echo true
  input:
  val x from ch_num
  val y from ch_letters
  script:
   """
   echo $x and $y
   """
}
~~~
{: .language-groovy }

In the above example the process is executed only two time, because when a queue channel has no more data to be processed it stops the process execution.

~~~
2 and b

1 and a
~~~
{: .output}

### Value channels and process termination

Note however that value channels, `Channel.value` , do not affect the process termination.

To better understand this behaviour compare the previous example with the following one:

~~~
ch_num = Channel.value(1)
ch_letters = Channel.of('a','b','c')

process combine_channels_val_queue {
  echo true
  input:
  val x from ch_num
  val y from ch_letters
  script:
   """
   echo $x and $y
   """
}
~~~
{: .language-groovy }

In this example the process is run three times.

~~~
1 and b
1 and a
1 and c
~~~
{: .output}


> ## Exercise Combining input channels
> Write a nextflow script `salmon_index.nf` that combines two input channels
> 1. transcriptome_ch = channel.value('data/yeast/transcriptome/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
> 2. kmer_ch = channel.value(21)
>  And run the command `salmon index -t $transcriptome -i index -k $kmer`
> > Solution
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.value(21)
> > process combine {
> >  input:
> >  path transcriptome from transcriptome_ch
> >  val kmer from kmer_ch
> >  script:
> >   """
> >   salmon index -t $transcriptome -i index -k $kmer
> >   """
> > }
> > ~~~
> >
> {: .solution}
{: .challenge}

# Input repeaters

We saw previously that by default the number of time a process run is defined by the queue channel with the fewest items. However, the `each` qualifier allows you to repeat the execution of a process for each item in a list or a queue channel, every time new data is received. For example if we can the previous example by using the input qualifer `each` for the letters queue channel:

~~~
ch_num = Channel.of(1,2)
ch_letters = Channel.of('a','b','c','d')

process combine_channels {
  input:
  val x from ch_num
  each y from ch_letters
  script:
   """
   echo $x and $y
   """
}
~~~
{: .language-groovy }

The process will run eight times.

~~~
2 and d
1 and a
1 and c
2 and b
2 and c
1 and d
1 and b
2 and a
~~~
{: .output}

> ## Exercise: Input repeaters
> Extend the previous Combining example by adding more values in the `kmer` queue channel  
> `kmer_ch = channel.of(21,26,36)`
> and changing the `kmer` input qualifer to `each`.
>> ## Solution
>> ~~~
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.of(21,26,36)
> > process combine {
> >  input:
> >  path transcriptome from transcriptome_ch
> >  each kmer from kmer_ch
> >  script:
> >   """
> >   echo salmon index -t $transcriptome -i index -k $kmer
> >   """
> > }
>>
> {: .solution}
{: .challenge}

## Outputs

We have seen how to input data into a process now we will see how to output file and values from a process.

The `output` declaration block allows us to define the channels used by the process to send out the file and values produced.

An output block is not required, but if it is present it can contain one or more outputs declarations. The output block follows the syntax shown below:

~~~
output:
  <output qualifier> <output name> into <target channel>[,channel,..]
  <output qualifier> <output name> into <target channel>[,channel,..]
  ...
~~~  
{: .language-groovy }

As you can see from above the output block use `into` whilst the input block uses the `from` keyword.

### Output values

The type of output data is defined using output qualifiers.

The `val` qualifier allows us to output a value defined in the script.
Because Nextflow processes can only communicate through channels if we want to share a value input into one process as input to another process we would  need to define that value in the output declaration block as shown in the following example:


~~~
methods_ch = channel.of('salmon','kallisto')

process method_type {
  input:
  val x from methods_ch

  output:
  val x into out_ch

  script:
  """
  echo $x > method.txt
  """
}

// use the view operator to display contents of the channel
out_ch.view({ "Received: $it" })
~~~
{: .language-groovy }

### Output files

If we want to capture a file instead of a value we can use the
`path` qualifier that can capture one or more files produced by the process, over the specified channel.

~~~
methods_ch = channel.of('salmon','kallisto')

process method_type {
  input:
  val x from methods_ch

  output:
  path 'method.txt' into out_ch

  """
  echo $x > method.txt
  """
}

// use the view operator to display contents of the channel
out_ch.view({ "Received: $it" })
~~~
{: .language-groovy }

In the above example the process `method_type` creates a file named `method.txt` containing the method name.

Since a file parameter using the same name, `method.txt`, is declared in the output block , when the task is completed that file is sent over the `out_ch` channel. A downstream `operator` or `process` declaring the same channel as input will be able to receive it.

### Multiple output files

When an output file name contains a `*` or `?` character it is interpreted as a pattern match.
This allows to capture multiple files into a list and output them as a one item channel.

For example here we will capture `sample.bam.bai` in the output channel. Note that `sample.bam` is not captured in the output channel as input files are not included in the list of possible matches

>>> do we need to include this example
~~~

bam_ch = channel.fromPath("data/yeast/bams/*.bam")

process index {
    input:
    path 'sample.bam' from bam_ch

    output:
    path "sample.bam*" into index_out_ch

    script:
    """
    samtools index sample.bam
    """
}
/*
*The flatMap operator applies a function to every item emitted by a channel, and returns the items so obtained as a new channel
*/
index_out_ch
    .flatMap()
    .view({ "File: ${it.name}" })
~~~
{: .language-groovy }

it prints:

~~~
File: sample.bam.bai

File: sample.bam.bai

File: sample.bam.bai

File: sample.bam.bai

File: sample.bam.bai

File: sample.bam.bai

File: sample.bam.bai
[truncated]
~~~
{: .output}

Some caveats on glob pattern behaviour:

* Input files are not included in the list of possible matches.
* Glob pattern matches against both files and directories path.
* When a two stars pattern `**` is used to recourse across directories, only file paths are matched i.e. directories are not included in the result list.



### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic evaluated string which references values defined in the input declaration block or in the script global context. For example::

~~~
bam_ch = channel.fromPath("data/yeast/bams/*.bam")

process index {
    input:
    path bam from bam_ch

    output:
    path "${bam}*" into index_out_ch

    script:
    """
    samtools index $bam
    """
}

index_out_ch.view({ "File: ${it.name}" })
~~~
{: .language-groovy }

In the above example, each time the process is executed an index file is produced whose name depends on the actual value of the `bam` input.

fix this as removed reference to flatmap operator
~~~
File: ref3.bam.bai

File: ref1.bam.bai

File: etoh60_3.bam.bai

File: ref2.bam.bai

File: etoh60_1.bam.bai

File: etoh60_2.bam.bai
[truncated]
~~~
{: .output}

> # Exercise: output channels
> Modify the nextflow script process_exercise_output.nf to include an output blocks that captures the different index folders. use the view operator on the output channel.
> > solution
> > ~~~
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.of(21,26,36)
> > process combine {
>  input:
>  path transcriptome from transcriptome_ch
>  each kmer from kmer_ch
>  output:
>  path "index_${kmer}" into out_ch
>  script:
>  """
>   echo salmon index -t $transcriptome -i index_$kmer -k $kmer
>  """
> }
> out_ch.view()
> ~~~
> {: .solution}
{: .challenge}

### Composite inputs and outputs

So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the `tuple` qualifiers.

In tuples the first element is the grouping key and the second element is the list of files.

~~~
[group_key,[file1,file2,...]]
~~~

When using channel containing a tuple, the corresponding input declaration must be declared with a `tuple` qualifier, followed by definition of each element in the tuple.

~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
process tuple_input{
input:
  tuple val(sample_id), path(reads) from reads_ch
script:
"""
echo $sample_id
echo $reads
"""
}
~~~
{: .language-groovy }

outputs

~~~
ref1
ref1_1.fq.gz ref1_2.fq.gz
~~~
{: .output }

In the same manner output channel emitting tuple of values can be declared using the `tuple` qualifier following by the definition of each tuple element in the tuple.

In the code snippet below the output channel `bam_ch` would contain a tuple with
the grouping key value as the Nextflow variable `sample_id` and a single `sample.bam` file stored as a value list.

~~~
output:
  tuple val(sample_id), path('sample.bam') into bam_ch
~~~
{: .language-groovy }

An example can be seen in this script below.

~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

process pseudo_align {
  input:
    tuple val(sample_id), path(reads) from reads_ch
    each index from index_ch
  output:
    tuple val(sample_id), path("${sample_id}.bam") into bam_ch
  script:
  """
   salmon quant  -i $index \ -1 ${reads[0]} -2 ${reads[1]} -o $sample_id -l A \
  --writeMappings |samtools sort |samtools view -bS -o ${sample_id}.bam
  """
}

bam_ch.view()
~~~
{: .language-groovy }

> # Composite inputs and outputs
> Fill in the blanks for process_exercise_tuple.nf.
> > Solution
> >
> > reads_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')
> > process fastqc {
> > input:
> >   tuple val(sample_id), path(reads) from reads_ch
> > output:
> >    tuple val(sample_id), path("fastqc_out") into bam_ch
> >  script:
> >  """
> >  mkdir fastqc_out
> >  fastqc $reads -o fastqc_out -t 1
> >  """
> >}
> >
> > bam_ch.view()> >
> >~~~
> {: .solution}
{: .challenge}


## When

The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value; `true` or `false`.

It is useful to enable/disable the process execution depending the state of various inputs and parameters. For example:

~~~

chr_ch = channel.of(1..22,'X','Y')

process conditional {

  input:
  val chr from chr_ch

  when:
  chr <= 22

  script:
  """
  echo $chr
  """
}
~~~
{: .language-groovy }


## Directives

Directive declarations allow the definition of optional settings, like resources or file admin, that affect the execution of the current process without affecting the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e. `input`, `output`, etc).

Directives are commonly used to define the amount of computing resources to be used or extra information for configuration or logging purpose. For example:

~~~
chr_ch = channel.of(1..22,'X','Y')

process printchr {
  tag "$chr"
  cpus 1

  input:
  val chr from chr_ch

  script:
  """
  echo $chr
  echo number of cpus $task.cpus
  """
}
~~~
{: .language-groovy }

Uses the `tag` directive to allow you to associate each process execution with a custom label, so that it will be easier to identify them in the log file or in the trace execution report and the `cpus` directive allows you to define the number of CPU required for each the process’ task.

Another commonly used directive is memory specification `memory`; a complete list of directives is available at this [link](https://www.nextflow.io/docs/latest/process.html#directives).

> # Directives
> Modify the Nextflow script `process_exercise_directives.nf`
> Add a `tag` directive logging the sample_id in the execution output.
> Add a cpus directive to specify the number of cpus as 2.
> Change the -t option value to $task.cpus in the script directive.
> > ## solution
> > ~~~
> > reads_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')
> > process fastqc {
> > tag "$sample_id"
> > cpus 2
> >  input:
> >    tuple val(sample_id), path(reads) from reads_ch
> >  output:
> >     tuple val(sample_id), path("fastqc_out") into bam_ch
> >   script:
> >   """
> >   mkdir fastqc_out
> >   fastqc $reads -o fastqc_out -t $task.cpus
> >   """
> > }
> >
> > ~~~
> {: .solution}
{: .challenge}

## Organising outputs

### PublishDir directive

Nextflow manages intermediate results from the pipelines expected outputs independently.

Files created by a process are stored in a task specific working directory which is considered as a temporary. Normally this is under the `work` directory , that can be deleted upon completion.

The files you want the workflow to return as results need to be defined in the process output block and then the output directory specified using the directive [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

~~~
publishDir <directory>, parameter: value, parameter2: value ...
~~~

For example:

~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

process quant {
  publishDir "results/bams"
  input:
    tuple val(sample_id), path(reads) from reads_ch
    path index from index_ch
  output:
    tuple val(sample_id), path("${sample_id}.bam") into bam_ch
  script:
  """
  salmon quant  -i $index \ -1 ${reads[0]} -2 ${reads[1]} -o $sample_id -l A \
  --writeMappings |samtools sort |samtools view -bS -o ${sample_id}.bam
  """
}
bam_ch.view()
~~~
{: .language-groovy }

In the above example, the `publishDir "results/bams"`, will create a symbolic link to all the bam files output by the process `quant` to the directory path `results/bams`.

**The publishDir output is relative to the directory of the main Nextflow script.**

### publishDir parameters

publishDir directive can take optional parameters, for example the `mode` parameter can take the value `copy` to specify that you wish to copy the file to output directory rather than just a symbolic link to the files in the working directory. Full list [here](https://www.nextflow.io/docs/latest/process.html#publishdir).


###  Manage semantic sub-directories

You can use more then one `publishDir` to keep different outputs in separate directories. To specify which files to put in which output directory use the parameter `pattern` with the a glob  pattern that selects which files to publish from the overall set of output files.

For example:


~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

process pseudo_align {
  publishDir "results/bams", pattern: "*.bam"
  publishDir "results/quant", pattern: "**.sf"
  input:
    tuple val(sample_id), path(reads) from reads_ch
    path index from index_ch
  output:
    tuple val(sample_id), path("${sample_id}.bam") into bam_ch
    path "${sample_id}/quant.sf" into quant_ch
  script:
  """
  salmon quant  -i $index \ -1 ${reads[0]} -2 ${reads[1]} -o $sample_id -l A \
  --writeMappings |samtools sort |samtools view -bS -o ${sample_id}.bam
  """
}
~~~
{: .language-groovy }

The above example will create an output folder structure in the directory results, which contains a separate sub-directory for bam files, `pattern:"*.bam"` ,  and quant files, `pattern:"**.sf"`.


> # directives
>  Add a `publishDir` directive to the nextflow script `process_publishDir_exercise.nf` that saves
>  the index directory to the results folder .
> > ~~~
> >
> > params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
> > transcriptome_ch = channel.fromPath(params.transcriptome,checkIfExists: true)
> > process index {
> >  publishDir "results"
> >  input:
> >  path transcriptome from transcriptome_ch
> >
> >  output:
> >  path "index" into index_ch
> >
> >  script:
> >  """
> >    salmon index -t $transcriptome -i index
> >  """
> >}
> > ~~~
> {: .solution}
{: .challenge}

> ## Nextflow Patterns
> If you want to find out common structures of Nextflow process the [Nextflow Patterns page](http://nextflow-io.github.io/patterns/) collects some recurrent implementation patterns used in Nextflow applications.
{: .callout}

{% include links.md %}
