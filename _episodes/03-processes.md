---
title: "Processes"
teaching: 0
exercises: 0
questions:
- "What is a Nextflow process?"
- "How do I create a Nextflow process?"
- "How do I input data into processes|"
- "How do I output data from a process?"
- "How do I sepcify conditions for a process in order for it to execute?"
- "What are process directives?"
- "How do i save output from a process?"
objectives:
- "Understand a Nextflow process."
- "Create a nextflow process."
- "Use values and files as inputs to a process."
- "Use the when declaration to define a condition for process execution."
- "Understand what process directives are." 
keypoints:
- "A Nextflow process is an independent task/step in a workflow"
- "Processes"
- "Task output files are output from a process using the `PublishDir` directive"
---


# Processes

We now know how to create and use Channels to control data flows in Nextflow. We will now see how to process independent tasks within a workflow.

A `process` is the basic Nextflow computing primitive to execute foreign function i.e. custom scripts or tools.

> ## primitives
> In computing, language primitives are the simplest elements available in a programming language. A primitive is the smallest 'unit of processing' available to a programmer of a given machine, or can be an atomic element of an expression in a language.
> {: .callout}

*Processes can be thought of as a particular task/steps in a workflow, e.g. an alignment step in RNA-Seq analysis. Processes and are independent of each other (don't require another processes to execute) and can not communicate/write to each other . It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output which are Channels*

The process definition starts with keyword the `process`, followed by process name and finally the process `body` delimited by curly brackets `{}`. The process body must contain a string which represents the command or, more generally, a script that is executed by it.

A basic process looks like the following example:

~~~
process sayHello {
  """
  echo 'Hello world!'
  """
}
~~~
{: .source}

A process may contain five definition blocks, respectively: directives, inputs, outputs, when clause and finally the process script. The syntax is defined as follows:

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
{: .source}

* Zero, one or more process directives
* Zero, one or more process inputs
* Zero, one or more process outputs
* An optional boolean conditional to trigger the process execution
* The command to be executed


## Script

The `script` block is a string statement that defines the command that is executed by the process to carry out its task.

A process contains one and only one script block, and it must be the last statement when the process contains input and output declarations.

The script block can be a simple string or multi-line string. The latter simplifies the writing of non trivial scripts composed by multiple commands spanning over multiple lines. For example::


~~~
process example {
    script:
    """
    blastp -db /data/blast -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    blastdbcmd -db /data/blast -entry_batch top_hits > sequences
    """
}
~~~
{: .source}

By default the process command is interpreted as a **Bash** script. However any other scripting language can be used just simply starting the script with the corresponding [Shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) declaration. For example:

~~~
process pyStuff {
  script:
  """
  #!/usr/bin/env python

  x = 'Hello'
  y = 'world!'
  print "%s - %s" % (x,y)
  """
}
~~~

> This allows the compositing in the same workflow script of tasks using different programming languages which may better fit a particular job. However for large chunks of code is suggested to save them into separate files and invoke them from the process script.
{: .callout}


### Script parameters

Process script can be defined dynamically using variable values like in other string.

~~~
params.data = 'World'

process foo {
  script:
  """
  echo Hello $params.data
  """
}
~~~
{: .source}

> # String interpolation
> A process script can contain any string format supported by the Groovy programming language. This allows us to use string interpolation or multiline string as in the script above. Refer to [String interpolation](https://seqera.io/training/#_string_interpolation) for more information.
{: .callout}

> Since Nextflow uses the same Bash syntax for variable substitutions in strings, Bash environment variables need to be escaped using `\` character.
{: .callout}

~~~
process foo {
  script:
  """
  echo "The current directory is \$PWD"
  """
}
~~~
{: .source}

> ## Escape Bash
>
> Try to modify the above script using $PWD instead of \$PWD and check the difference.
.
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

However this won’t allow any more the usage of Nextflow variables in the command script.

Another alternative is to use a `shell` statement instead of `script` which uses a different syntax for Nextflow variable: `!{..}`. This allow to use both Nextflow and Bash variables in the same script.

```
params.data = 'le monde'

process baz {
  shell:
  '''
  X='Bonjour'
  echo $X !{params.data}
  '''
}
```
{: .source}

### Conditional script

The process script can also be defined in a complete dynamic manner using a if statement or any other expression evaluating to string value. For example:

~~~
params.aligner = 'kallisto'
params.transcriptome = "$baseDir/data/yeast/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa"

process foo {
  script:
  if( params.aligner == 'kallisto' )
    """
    kallisto index -i index $transcriptome
    """
  else if( params.aligner == 'salmon' )
    """
    salmon index -t $transcriptome -i index
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.aligner")
}
~~~
{: .source}

> ## Conditional Exercise
>
> Write a custom function that given the aligner name as parameter returns the command string to be executed. Then use this function as the process script body.
.
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

## Inputs

Nextflow processes are isolated from each other but can communicate between themselves sending values through channels.

Inputs implicitly determine the dependency and the parallel execution of the process. The process execution is fired each time a new data is ready to be consumed from the input channel:

![Process Flow](../figs/channel-process.png)


The input block defines which channels the process is expecting to receive inputs data from. You can only define one input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below:

~~~
input:
  <input qualifier> <input name> from <source channel>
~~~
{: .source}

### Input values

The val qualifier allows you to receive data of any type as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

~~~
num = Channel.from( 1, 2, 3 )

process basicExample {
  input:
  val x from num

  """
  echo process job $x
  """
}
~~~
{: .source}

In the above example the process is executed three times, each time a value is received from the channel num and used to process the script. Thus, it results in an output similar to the one shown below:

~~~
process job 3
process job 1
process job 2
~~~
{: .output}

> ## Channel order
> The channel guarantees that items are delivered in the same order as they have been sent - but - since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received.
{: .callout}

### Input files

The `file` qualifier allows the handling of file values in the process execution context. This means that Nextflow will stage it in the process execution directory, and it can be access in the script by using the name specified in the input declaration.

~~~
reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file 'sample.fastq' from reads
    script:
    """
    your_command --reads sample.fastq
    """
}
~~~
{: .source}

The input file name can also be defined using a variable reference as shown below:

~~~
reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file sample from reads.collect()
    script:
    """
    your_command --reads $sample
    """
}
~~~
{: .source}

> # File Objects as inputs
> When a process declares an input file the corresponding channel elements must be file objects i.e. created with the file helper function from the file specific channel factories e.g. `Channel.fromPath` or `Channel.fromFilePairs`.
{: .callout}

Consider the following snippet:


~~~
params.genome = 'data/ggal/transcriptome.fa'

process foo {
    input:
    file genome from params.genome
    script:
    """
    your_command --reads $genome
    """
}
~~~
{: .source}

The above code creates a temporary file named input.1 with the string data/ggal/transcriptome.fa as content. That likely is not what you wanted to do.

### Input path

As of version 19.10.0, Nextflow introduced a new `path` input qualifier that simplifies the handling of cases such as the one shown above. In a nutshell the input path automatically handles string values as file objects. The following example works as expected:


~~~
params.genome = "$baseDir/data/ggal/transcriptome.fa"

process foo {
    input:
    path genome from params.genome
    script:
    """
    your_command --reads $genome
    """
}
~~~
{: .source}

> # Path qualifier
> The path qualifier should be preferred over file to handle process input files when using Nextflow 19.10.0 or later.
{: .callout}



> > Exercise
> > Write a script that creates a channel containing all read files matching the pattern `data/ggal/*_1.fq` followed by a process that concatenates them into a single file and prints the first 20 lines.

~~~
reads = Channel.fromPath( 'data/ggal/*_1.fq' )

process foo {
    input:
    file sample from reads.collect()
    script:
    """
    head -n 20 $sample > combined_n20.txt
    """
}
~~

### Combine input channels

A key feature of processes is the ability to handle inputs from multiple channels. However it’s important to understands how the content of channel and their semantic affect the execution of a process.

Consider the following example:
~~~
process foo {
  echo true
  input:
  val x from Channel.from(1,2,3)
  val y from Channel.from('a','b','c')
  script:
   """
   echo $x and $y
   """
}
~~~
{: .callout}

Both channels emit three value, therefore the process is executed three times, each time with a different pair:

~~~
2 and b

1 and a

3 and c
~~~
{: .output}

What is happening is that the process waits until there’s a complete input configuration i.e. it receives an input value from all the channels declared as input.

When this condition is verified, it consumes the input values coming from the respective channels, and spawns a task execution, then repeat the same logic until one or more channels have no more content.

This means channel values are consumed serially one after another and the first empty channel cause the process execution to stop even if there are other values in other channels.

**What does it happen when not all channels have the same cardinality (i.e. they emit a different number of elements)?**

For example:

~~~
process foo {
  echo true
  input:
  val x from Channel.from(1,2)
  val y from Channel.from('a','b','c','d')
  script:
   """
   echo $x and $y
   """
}
~~~
{: .source}

In the above example the process is executed only two time, because when a channel has no more data to be processed it stops the process execution.

~~~
2 and b

1 and a
~~~
{: .output}

> ## Value channels and process termination
> Note however that value channel do not affect the process termination.
> {: .output}

To better understand this behavior compare the previous example with the following one:

~~~
process bar {
  echo true
  input:
  val x from Channel.value(1)
  val y from Channel.from('a','b','c')
  script:
   """
   echo $x and $y
   """
}
~~~

## Exercise Combine input channels
Write a process that is executed for each read file matching the pattern data/ggal/*_1.fq and use the same data/ggal/transcriptome.fa in each execution.
remember value channels can be read multiple times.
~~~
reads_ch = Channel.fromPath('data/ggal/*_1.fq')
transcriptome_ch = channel.value('ggal/transcriptome.fa')
process combine {
  echo true
  input:
  path(y) from reads_ch
  path(x) from transcriptome_ch
  script:
   """
   echo $x and $y
   """
}
~~~

{: .source}

# Input repeaters

The `each` qualifier allows you to repeat the execution of a process for each item in a collection, every time a new data is received. For example:

~~~
sequences = Channel.fromPath('data/prots/*.tfa')
methods = ['regular', 'expresso', 'psicoffee']

process alignSequences {
  input:
  path seq from sequences
  each mode from methods

  """
  t_coffee -in $seq -mode $mode
  """
}
~~~
{: .source}

In the above example every time a file of sequences is received as input by the process, it executes three tasks running an alignment with a different value for the `mode` option. This is useful when you need to repeat the same task for a given set of parameters.

>> ##Exercise
>> Extend the previous example so a task is executed for each read file matching the pattern data/ggal/*_1.fq and repeat the same task both with salmon and kallisto.
>> ~~~
>> sequences = Channel.fromPath('data/raw_reads/SRR4204500/*.fastq.gz')
>> kmers = [21, 19, 31]
>>
>> process alignSequences {
>>  input:
>>  path seq from sequences
>>  each kmer from kmers
>>
>>  """
>>  echo $kmer $seq
>>  """
>> }
>> ~~~
methods = ['regular', 'expresso', 'psicoffee']

## Outputs

The `output` declaration block allows to define the channels used by the process to send out the results produced.

There can be defined at most one output block and it can contain one or more outputs declarations. The output block follows the syntax shown below:

~~~
output:
  <output qualifier> <output name> into <target channel>[,channel,..]
~~~  
{: .source}

### Output values

The `val` qualifier allows to output a value defined in the script context. In a common usage scenario, this is a value which has been defined in the *input* declaration block, as shown in the following example::


~~~
methods = ['prot','dna', 'rna']

process foo {
  input:
  val x from methods

  output:
  val x into receiver

  """
  echo $x > file
  """
}

receiver.view { "Received: $it" }
~~~
{: .source}

### Output files

The `file` qualifier allows to output one or more files, produced by the process, over the specified channel.

~~~
process randomNum {

    output:
    file 'result.txt' into numbers

    '''
    echo $RANDOM > result.txt
    '''
}

numbers.view { "Received: " + it.text }
~~~

In the above example the process r`andomNum` creates a file named `result.txt` containing a random number.

Since a file parameter using the same name is declared in the output block, when the task is completed that file is sent over the `numbers` channel. A downstream `process` declaring the same channel as input will be able to receive it.

### Multiple output files

When an output file name contains a `*` or `?` wildcard character it is interpreted as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matcher. This allows to capture multiple files into a list object and output them as a sole emission. For example:

~~~
process splitLetters {

    output:
    file 'chunk_*' into letters
    
    script:
    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}
/*
*The flatMap operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel
*/
letters
    .flatMap()
    .view { "File: ${it.name} => ${it.text}" }
~~~
{: .source}

it prints:

~~~
File: chunk_aa => H
File: chunk_ab => o
File: chunk_ac => l
File: chunk_ad => a
~~~
{: .output}

Some caveats on glob pattern behavior:

* Input files are not included in the list of possible matches.
* Glob pattern matches against both files and directories path.
* When a two stars pattern `**` is used to recourse across directories, only file paths are matched i.e. directories are not included in the result list.

Exercise
Remove the flatMap operator and see out the output change. The documentation for the flatMap operator is available at this [link](https://www.nextflow.io/docs/latest/operator.html#flatmap).

### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic evaluated string which references values defined in the input declaration block or in the script global context. For example::

~~~
process align {
  input:
  val species_name from species
  file seq from sequences

  output:
  file "${species_name}.aln" into genomes
  
  script:
  """
  t_coffee -in $seq > ${species_name}.aln
  """
}
~~~
{: .source}

In the above example, each time the process is executed an alignment file is produced whose name depends on the actual value of the `species_name` input.

### Composite inputs and outputs

So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the `tuple` qualifiers.

When using channel emitting tuple of values the corresponding input declaration must be declared with a `tuple` qualifier followed by definition of each single element in the tuple.

In the same manner output channel emitting tuple of values can be declared using the `tuple` qualifier following by the definition of each tuple element in the tuple.




Exercise
Modify the script of the previous exercise so that the bam file is named as the given sample_id.

## When

The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending the state of various inputs and parameters. For example:


~~~
params.dbtype = 'nr'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process find {
  input:
  file fasta from proteins
  val type from params.dbtype

  when:
  fasta.name =~ /^BB11.*/ && type == 'nr'

  script:
  """
  blastp -query $fasta -db nr
  """
}
~~~
{: .source}


## Directives

Directive declarations allow the definition of optional settings that affect the execution of the current process without affecting the semantic of the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e. `input`, `output`, etc).

Directives are commonly used to define the amount of computing resources to be used or other meta directives like that allows the definition of extra information for configuration or logging purpose. For example:

~~~
process foo {
  cpus 2
  memory 8.GB
  container 'image/name'

  script:
  """
  your_command --this --that
  """
}
~~~
{: .source}

The complete list of directives is available at this [link](https://www.nextflow.io/docs/latest/process.html#directives).

Exercise

Modify the script of the previous exercise adding a `tag` directive logging the sample_id in the execution output.


## Organise outputs

## PublishDir directive

Nextflow manages independently workflow execution intermediate results from the pipeline expected outputs. Task output files are created in the task specific execution directory which is considered as a temporary directory that can be deleted upon completion.

The pipeline result files need to be marked explicitly using the directive [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) in the process that’s creating such file. For example:

~~~
process makeBams {
    publishDir "/some/directory/bam_files", mode: 'copy'

    input:
    file index from index_ch
    tuple val(name), file(reads) from reads_ch

    output:
    tuple val(name), file ('*.bam') into star_aligned

    """
    STAR --genomeDir $index --readFilesIn $reads
    """
}
~~~
{: .source}

> ## Nextflow Patterns
> The [Nextflow Patterns page](http://nextflow-io.github.io/patterns/) collects some recurrent implementation patterns used in Nextflow applications. 
{: .callout}

{% include links.md %}

