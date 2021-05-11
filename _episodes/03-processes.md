---
title: "Processes"
teaching: 50 min
exercises: 20 min
questions:
- "How do I run commands/tasks in Nextflow?"
- "How do I get data, files and values, into and out of processes?"
- "How do I specify conditions for a process in order for it to execute?"
- "How do I control resources, such as number of CPUs and memory, for a process?"
- "How do I save output/results from a process?"
objectives:
- "Understand how Nextflow uses processes to execute tasks."
- "Create a Nextflow process."
- "Define inputs and outputs to a process."
- "Use the when declaration to define a condition for process execution."
- "Use process directives to control execution of a process."
- "Use the publishDir directive to save results in a directory."
keypoints:
- "A Nextflow process is an independent task/step in a workflow"
- "Processes contain up to five definition blocks including, directives, inputs, outputs, when clause and finaly a script block."
- "Inputs and Outputs to a process are defined using the input and output blocks."
- "The execution of a process can be controlled using conditional statements like `if`."
- "Task output files are output from a process using the `PublishDir` directive."
---


# Processes

We now know how to create and use Channels to send data around a workflow. We will now see how to run tasks within a workflow.

A `process` is the way Nextflow execute commands you would run on the command line or custom scripts.


Processes can be thought of as a particular task/steps in a workflow, e.g. an alignment step in RNA-Seq analysis. Processes  are independent of each other (don't require another processes to execute) and can not communicate/write to each other . It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output which are Channels

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

This example build a salmon index for the yeast transcriptome. We use the nextflow variable `${projectDir}` to specify the path to the directory where the script is run.


### definition blocks

A process may contain five definition blocks, that control how a command is executed within a process, these  are:

1. **directives**: allow the definition of optional settings that affect the execution of the current process e.g. the number of cpus a task uses and the amount of memory allocated.
1. **inputs**: Define the input dependencies, usually channels, which determines the number of times a process is executed.
1. **outputs**: Defines the output channels used by the process to send results/data produced by the process.
1. **when clause**: Allow you to define a condition that must be verified in order to execute the process.
1. **The script block**: The script block is a string statement that defines the command that is executed by the process to carry out its task.


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

The `script` block is a string statement that defines the command that is executed by the process to carry out its task. These are normally the commands you would run on a terminal.

A process contains only one script block, and it must be the last statement when the process contains input and output declarations.

The script block can be a simple string  or multi-line string. Multi-line strings simplifies the writing of non trivial scripts composed by multiple commands spanning over multiple lines. For example:


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

> This allows the the use of a different programming languages which may better fit a particular job. However for large chunks of code is suggested to save them into separate files and invoke them from the process script.
{: .callout}


### Script parameters

The command in the script block can be defined dynamically using Nextflow variables e.g. `projectDir`.
To reference a variable in the script block you can use the `$` in front of the nextflow variable name, and additionally you can add `{}` around the variable name e.g. `${projectDir}`.

A special nextflow variable `params` can be used to define variable from the command line.

~~~
params.kmer = 31

process index {

  script:
  """
  salmon index -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer $params.kmer
  echo "kmer size is $params.kmer"
  """
}
~~~
{: .language-groovy }

In this example we set the Nextflow variable `params.kmer` to value `31` and the process `index` uses this to set the  size of k-mers that should be used for the quasi index.

> ## Script parameters
>
> Change the kmer value used to create the salmon index from the command line.
> ~~~
> nextflow run process_script_params.nf --kmer
> Note the kmer values must not be greater than 31 and an odd number.
> ~~~
>
> > ## Solution
> > ~~~
> > nextflow run process_script_params.nf --kmer 27
> > ~~~
> {: .solution}
{: .challenge}


A process script can contain any string format supported by the Groovy programming language. In practice this allows us to use string interpolation, such as variable substitutions and use multiline string as in the script above. Refer to [String interpolation](https://seqera.io/training/#_string_interpolation) for more information.


Since Nextflow uses the same Bash syntax for variable substitutions, `$variable`, in strings, Bash variables need to be escaped using `\` character.

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


This will print the location of the working directory using the bash environment variable variable `PWD`.


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

Sometimes you want to change how a process is run depending on some parameter. In Nextflow scripts we can use conditional statements such as the `if` statement or any other expression evaluating to string value. For example, the nextflow script below will change what index is created depending on the nextflow variable `params.aligner`.

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

Nextflow processes are isolated from each other but can communicate by sending values through channels.

Inputs implicitly determine the dependency and the number of time a process executes. A process is run each time new data is ready to be consumed from the input channel:

![Process Flow](../fig/channel-process.png)


The `input` block defines which channels the process is expecting to receive input data from. You can only define one input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below:

~~~
input:
  <input qualifier> <input name> from <source channel>
~~~
{: .language-groovy }


The input qualifier declares the type of data to be received.

* `val`: Lets you access the received input value by its name in the process script.
* `env`: Lets you use the received value to set an environment variable named as the specified input name.
* `path`: Lets you handle the received value as a path, staging the file properly in the execution context.
* `stdin`: Lets you forward the received value to the process stdin special file.
* `tuple`: Lets you handle a group of input values having one of the above qualifiers.
* `each`: Lets you execute the process for each entry in the input collection.


### Input values

The `val` qualifier allows you to receive value data as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

~~~
chr_ch = Channel.of( 1..22,'X','Y' )

process printChr {

  input:
  val chr from chr_ch

  """
  echo process chromosome $chr
  """
}
~~~
{: .language-groovy }

In the above example the process is executed 24 times, each time a value is received from the channel `chr_ch` and used to process the script. Thus, it results in an output similar to the one shown below:

~~~
process chromosome 3
process chromosome 1
process chromosome 2
..truncated...
~~~
{: .output}

> ## Channel order
> The channel guarantees that items are delivered in the same order as they have been sent,  but,  since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received.
{: .callout}

### Input files

The `path` qualifier allows the handling of files. This means that Nextflow will stage it in the process execution directory, and it can be access in the script by using the name specified in the input declaration. In the example below we specific the name of the file as `sample.fastq`.

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


The input file name can also be defined dynamically using a Nextflow variable reference, `path sample`, as shown below:

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
>  Write a nextflow script `fastqc.nf` that creates a queue channel containing all read files matching the pattern `data/yeast/reads/*_1.fq.gz` followed by a process the commands.
> `mkdir fastqc_out`
> `fastqc -o fastqc_out ${reads}`
> > Solution
> > ~~~
> > reads = Channel.fromPath( 'data/yeast/reads/*_1.fq.gz' )
> >
> > process foo {
> >    input:
> >    file sample from reads
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
However it’s important to understands how the content of channel and their semantic affect the execution of a process.

Consider the following example:

~~~
ch_num = Channel.of(1,2,3)
ch_letters = Channel.of('a','b','c')

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

Both channels emit three value, therefore the process is executed three times, each time with a different pair:

~~~
2 and b

1 and a

3 and c
~~~
{: .output}

What is happening is that the process waits until it receives an input value from all the channels declared as input.

When this condition is verified, it consumes the input values coming from the respective channels, runs the task. This logic repeats until one or more channels have no more content.

This means items within a channel are consumed  one after another and the first empty channel cause the process execution to stop even if there are other values in other channels.

**What happens when not all channels have the same number of elements?**

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

In the above example the process is executed only two time, because when a channel has no more data to be processed it stops the process execution.

~~~
2 and b

1 and a
~~~
{: .output}

### Value channels and process termination

Note however that value channels do not affect the process termination.

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

> ## Exercise Combing input channels
> Write a nextflow script `salmon_index.nf` that combines two input channels
> 1. transcriptome_ch = channel.value('data/yeast/transcriptome/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
> 2. kmer_ch = channel.value(21)
>  And run the command salmon index -t $transcriptome -i index -k $kmer
> > Solution
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.value(21)
> > process combine {
> >  input:
> >  path transcriptome from transcriptome_ch
> >  val kmer from kmer_ch
> >  script:
> >   """
> >   echo salmon index -t $transcriptome -i index -k $kmer
> >   """
> > }
> > ~~~
> >
> {: .solution}
{: .challenge}

# Input repeaters

We saw previously that by default the number of time a process run is defined by the queue channel with the fewest items.
The `each` qualifier allows you to repeat the execution of a process for each item in a list or a queue channel, every time new data is received. For example if we can the previous example by using the input qualifer `each` for the letters queue channel:

~~~
ch_num = Channel.of(1,2)
ch_letters = Channel.of('a','b','c','d')

process combine_channels {
  echo true
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

> ## Exercise: Input repeaters
> Extend the previous combing example by adding more values in the `kmer` queue channel  
> `kmer_ch = channel.value(21,26,36)`
> and changing the `kmer` input qualifer to `each`.
>> ## Solution
>> ~~~
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.value(21)
> > process combine {
> >  input:
> >  path transcriptome from transcriptome_ch
> >  val kmer from kmer_ch
> >  script:
> >   """
> >   echo salmon index -t $transcriptome -i index -k $kmer
> >   """
> > }
>>
> {: .solution}
{: .challenge}

## Outputs

We have seen how to input data into a process now we will see how to output results from a process.

The `output` declaration block allows us to define the channels used by the process to send out the results produced.

There can be defined at most one output block and it can contain one or more outputs declarations. The output block follows the syntax shown below:

~~~
output:
  <output qualifier> <output name> into <target channel>[,channel,..]
~~~  
{: .language-groovy }

### Output values

The type of output data is defined using output qualifiers.

The `val` qualifier allows to output a value defined in the script context. In a common usage scenario, this is a value which has been defined in the *input* declaration block, as shown in the following example::


~~~
methods_ch = channel.of('salmon','kallisto')

process method_type {
  input:
  val x from methods_ch

  output:
  val x into out_ch

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

When an output file name contains a `*` or `?` wildcard character it is interpreted as a pattern match.
This allows to capture multiple files into a list object and output them as a sole emission.

For example here we will capture sample.bam.bai in the output channel. Note that sample.bam is captured in the output channel
as input files are not included in the list of possible matches

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
File: sample.bam
File: sample.bam.bai
~~~
{: .output}

Some caveats on glob pattern behaviour:

* Input files are not included in the list of possible matches.
* Glob pattern matches against both files and directories path.
* When a two stars pattern `**` is used to recourse across directories, only file paths are matched i.e. directories are not included in the result list.



### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic evaluated string which references values defined in the input declaration block or in the script global context. For example::
fixme better example
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
/*
*The flatMap operator applies a function to every item emitted by a channel, and returns the items so obtained as a new channel
*/
index_out_ch
    .flatMap()
    .view({ "File: ${it.name}" })
~~~
{: .language-groovy }

In the above example, each time the process is executed an index file is produced whose name depends on the actual value of the `bam` input.


> # Exercise: output channels
> Modify the nextflow script process_exercise_3.nf to include an output defintion blocks that captures the different index folders.
> > solution
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.of(21,26,36)
> > process combine {
> .  input:
  path transcriptome from transcriptome_ch
  each kmer from kmer_ch
  script:
  """
   echo salmon index -t $transcriptome -i index_$kmer -k $kmer
  """
 }
> {: .solution}
{: .challenge}

### Composite inputs and outputs

So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the `tuple` qualifiers.

When using channel emitting a tuple, the corresponding input declaration must be declared with a `tuple` qualifier, followed by definition of each single element in the tuple.


~~~
input:
  tuple val(sample_id), path(reads) from reads_ch
~~~
{: .language-groovy }


In the same manner output channel emitting tuple of values can be declared using the `tuple` qualifier following by the definition of each tuple element in the tuple.

~~~
output:
  tuple val(sample_id), path('sample.bam') into bam_ch
~~~
{: .language-groovy }

~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

process pseudo_align {
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

> #Exercise
> Modify the script of the previous exercise so that the bam file is named as the given sample_id.
> > Solution
> >
> > reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')
> > process foo {
> >  input:
> >    tuple val(sample_id), file(sample_files) from reads_ch
> >  output:
> >    tuple val(sample_id), file('sample.bam') into bam_ch
> >  script:
> >  """
> >    echo align --reads $sample_id > ${sample_id}.bam
> >  """
> >}
> >
> >bam_ch.view()
> >
> >~~~
> {: .solution}
{: .challenge}

## When


The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value, true or false.

It is useful to enable/disable the process execution depending the state of various inputs and parameters. For example:

~~~

chr_ch = channel.of(1..22,'X','Y')

process conditional {
  echo true
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

Directive declarations allow the definition of optional settings that affect the execution of the current process without affecting the semantic of the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e. `input`, `output`, etc).

Directives are commonly used to define the amount of computing resources to be used or other meta directives like that allows the definition of extra information for configuration or logging purpose. For example:

~~~
chr_ch = channel.of(1..22,'X','Y')

process printchr {
  tag "$chr"
  echo true
  cpus 1

  input:
  val chr from chr_ch

  script:
  """
  sleep 2
  echo $chr
  """
}
~~~
{: .language-groovy }

The complete list of directives is available at this [link](https://www.nextflow.io/docs/latest/process.html#directives).

> # Exercise
> Modify the script of the previous exercise adding a `tag` directive logging the sample_id in the execution output.
> > ## solution
> > ~~~
> > chr_ch = channel.of(1..22,'X','Y')
> >
> > process printchr {
> >  tag "process chromsome : $chr"
> >  label 'big_mem'
> >  cpus 2
> >  memory 2.GB
> >  //container 'image/name'
> >
> >  input:
> >  val chr from chr_ch
> >
> >  script:
> >  """
> >  echo $chr
> >  """
> >}
> > ~~~
> {: .solution}
{: .challenge}

## Organising outputs

### PublishDir directive

Nextflow manages independently workflow execution intermediate results from the pipeline expected outputs. Task output files are created in the task specific execution directory which is considered as a temporary directory that can be deleted upon completion.

The pipeline result files need to be marked explicitly using the directive [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) in the process that’s creating such file. For example:

~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

process quant {
  publishDir "results/bams", mode: "copy"
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

The above example will copy all bam files created by the  process quant to the directory path `results/bams`.

###  Manage semantic sub-directories

You can use more then one publishDir to keep different outputs in separate directory. For example:


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

The above example will create an output structure in the directory results, which contains a separate sub-directory for bam files  and quant files.


> # directives
>  Add a `publishDir` directive to the nextflow script `process_publishDir_exercise.nf` that saves
>  index directory to the results folder .
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
> The [Nextflow Patterns page](http://nextflow-io.github.io/patterns/) collects some recurrent implementation patterns used in Nextflow applications.
{: .callout}

{% include links.md %}
