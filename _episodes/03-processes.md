---
title: "Processes"
teaching: 50 min
exercises: 20 min
questions:
- "How do I run commands/tasks in Nextflow?"
- "How do I define a process in Nextflow?"
- "How do I get data, files and values, into and out of processes?"
- "How do I specify conditions for a process in order for it to execute?"
- "How do i control resources for a process?"
- "How do i save output from a process?"
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
- "Task output files are output from a process using the `PublishDir` directive."
---


# Processes

We now know how to create and use Channels in Nextflow. We will now see how to run tasks within a workflow.

A `process` is the way Nextflow execute commands you would run on the terminal or custom scripts.


*Processes can be thought of as a particular task/steps in a workflow, e.g. an alignment step in RNA-Seq analysis. Processes and are independent of each other (don't require another processes to execute) and can not communicate/write to each other . It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output which are Channels*

The process definition starts with keyword the `process`, followed by process name and finally the process `body` delimited by curly brackets `{}`. The process body must contain a string which represents the command or, more generally, a script that is executed by it.

A basic process with no input or output channels looks like the following example:

~~~
process printWorkingDirectory {
  echo true
  script:
  """
  pwd
  ls -a
  """
}
~~~
{: .language-groovy }

This example runs two unixs command `pwd` to show the directory the command was executed in and `ls -a` to list all the files in the directory . The directive `echo true` is added so that command is printed to the screen. This Produces

~~~
N E X T F L O W  ~  version 20.10.0
Launching `delme.nf` [distraught_miescher] - revision: cdc9fc9d97
executor >  local (1)
[1f/11a6b7] process > printWorkingDirectory [100%] 1 of 1 ✔
work/1f/11a6b7071a4f7852c4ca399166f58a
.
..
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
~~~
{: .output}

A process may contain five definition blocks, these are:

1. **directives**: allow the definition of optional settings that affect the execution of the current process e.g. the number of cpus a task uses.
1. **inputs**: Define the input dependencies which determines the parallel execution of the process.
1. **outputs**: Define the channels used by the process to send out the results produced.
1. **when clause**: Allow you to define a condition that must be verified in order to execute the process
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

The `script` block is a string statement that defines the command that is executed by the process to carry out its task.

A process contains only one script block, and it must be the last statement when the process contains input and output declarations.

The script block can be a simple string or multi-line string. The latter simplifies the writing of non trivial scripts composed by multiple commands spanning over multiple lines. For example:


~~~
process example {
    script:
    """
    samtools sort  ${prefix}.sorted.bam -T $name $bam
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
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

  x = 'Hello'
  y = 'world!'
  print ("%s - %s" % (x,y))
  """
}
~~~
{: .language-groovy }

> This allows the the use of a different programming languages which may better fit a particular job. However for large chunks of code is suggested to save them into separate files and invoke them from the process script.
{: .callout}


### Script parameters

The Process script block can be defined dynamically using Nextflow variables.
To reference a variable in the script block you can use the `$` in front of the variable name and additionally add `{}` around the variable name e.g. `${params.genome}`.

~~~
//nextflow params variable
params.genome = 'GRCh38'

process print_nf_variable {
  echo true
  script:
  """
  echo $params.genome
  """
}
~~~
{: .language-groovy }

In this example we set the Nextflow variable `params.genome` to `GRCh38` and the process
print_nf_variable echos the contents of the variable to the screen.
~~~
N E X T F L O W  ~  version 20.10.0
Launching `delme.nf` [nice_almeida] - revision: d50659d3ab
executor >  local (1)
[f1/378634] process > print_nf_variable [100%] 1 of 1 ✔
GRCh38
~~~
{: .output }

A process script can contain any string format supported by the Groovy programming language. This allows us to use string interpolation or multiline string as in the script above. Refer to [String interpolation](https://seqera.io/training/#_string_interpolation) for more information.


Since Nextflow uses the same Bash syntax for variable substitutions in strings, Bash environment variables need to be escaped using `\` character.


~~~
process printWorkingDirectory {
  echo true
  script:
  """
  echo \$PWD
  """
}
~~~
{: .language-groovy }


This will print the location of the current working directory using the bash environment variable variable `PWD`.


~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [fervent_varahamihira] - revision: 2268ae3939
executor >  local (1)
[1b/202631] process > listFileInCWD [100%] 1 of 1 ✔
work/1b/202631ace5f3647972e8ddbdb0331c
~~~
{: .output}

> ## Escape Bash
>
> Try to modify the above script using the nextflow variable `$PWD` instead of `\$PWD` and check the difference.
>
> > ## Solution
> > If you do not escape the BASH variable PWD it will use the
> > This is the body of the solution.
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `delme.nf` [small_cantor] - revision: 0252a6773a
> > executor >  local (1)
> > [47/18a9f8] process > printWorkingDirectory [100%] 1 of 1 ✔
> > /Users/ggrimes2/Documents/nextflow-training>
> > ~~~
> > This will print the contents of the Nextflow variable `PWD` which is the location od the directory of the nextflow script.
> {: .solution}
{: .challenge}

However this won’t allow any more the usage of Nextflow variables in the command script.

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


> ## process parameters
>
> 1. Write a nextflow script `process_ex1.nf` that executes the following bash command in the script block.
> ~~~
> echo "The nextflow parameters are $params"
~ ~~~
> Run the command below
> ~~~~
> nextflow run process_ex1.nf --genome "GRCh38"
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

### Conditional execution

Sometimes you want to change how a process is run depending on some parameter. In Nextflow scripts we can use conditional statements such as the `if` statement or any other expression evaluating to string value. For example:

~~~
params.aligner = 'kallisto'

process align {
  script:
  if( params.aligner == 'kallisto' )
    """
    echo Align using kallisto
    """
  else if( params.aligner == 'salmon' )
    """
    echo Align using kallisto
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.aligner")
}
~~~
{: .language-groovy }



## Inputs

Nextflow processes are isolated from each other but can communicate between themselves sending values through channels.

Inputs implicitly determine the dependency and the parallel execution of the process. The process execution is run each time new data is ready to be consumed from the input channel:

![Process Flow](../fig/channel-process.png)


The `input` block defines which channels the process is expecting to receive inputs data from. You can only define one input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below:

~~~
input:
  <input qualifier> <input name> from <source channel>
~~~
{: .language-groovy }


The input qualifier declares the type of data to be received.

* val: Lets you access the received input value by its name in the process script.
* env: Lets you use the received value to set an environment variable named as the specified input name.
* path: Lets you handle the received value as a path, staging the file properly in the execution context.
* stdin: Lets you forward the received value to the process stdin special file.
* tuple: Lets you handle a group of input values having one of the above qualifiers.
* each: Lets you execute the process for each entry in the input collection.


### Input values

The `val` qualifier allows you to receive data of any type as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

~~~
chr_ch = Channel.of( 1..21,'X','Y' )

process basicExample {
  echo true
  input:
  val chr from chr_ch

  """
  echo process chromosome $chr
  """
}
~~~
{: .language-groovy }

In the above example the process is executed 23 times, each time a value is received from the channel `chr_ch` and used to process the script. Thus, it results in an output similar to the one shown below:

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

The `path` qualifier allows the handling of file in the process execution context.
This means that Nextflow will stage it in the process execution directory, and it can be access in the script by using the name specified in the input declaration. In the example below we specific the name of the file as `sample.fastq`.

~~~
reads = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )

process foo {
    echo true

    input:
    path 'sample.fastq' from reads
    script:
    """
    ls -l sample.fastq
    """
}
~~~
{: .language-groovy }


The input file name can also be defined dynamically using a variable reference, `path sample`, as shown below:

~~~
reads = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )

process foo {
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
>  Write a nextflow script that creates a channel containing all read files matching the pattern `data/yeast/reads/*_1.fq.gz` followed by a process the commands.
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

fixme add a diagram

A key feature of processes is the ability to handle inputs from multiple channels.
However it’s important to understands how the content of channel and their semantic affect the execution of a process.

Consider the following example:

~~~
ch_num = Channel.of(1,2,3)
ch_letters = Channel.of('a','b','c')

process foo {
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

What is happening is that the process waits until there’s a complete input configuration i.e. it receives an input value from all the channels declared as input.

When this condition is verified, it consumes the input values coming from the respective channels, and spawns a task execution, then repeat the same logic until one or more channels have no more content.

This means channel values are consumed  one after another and the first empty channel cause the process execution to stop even if there are other values in other channels.

**What happens when not all channels have the same number of elements?**

For example:

~~~

ch_num = Channel.of(1,2)
ch_letters = Channel.of('a','b','c','d')

process foo {
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
ch_letters = Channel.from('a','b','c')

process bar {
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

> ## Exercise Combine input channels
> Write a process that is executed for each read file matching the pattern data/ggal/*_1.fq and
> use the same data/ggal/transcriptome.fa in each execution.
> Remember value channels can be read multiple times.
> > Solution
> > reads_ch = Channel.fromPath('data/ggal/*_1.fq')
> > transcriptome_ch = channel.value('ggal/transcriptome.fa')
> > process combine {
> >  echo true
> >  input:
> >  path y from reads_ch
> >  path x from transcriptome_ch
> >  script:
> >   """
> >   echo $x and $y
> >   """
> > }
> > ~~~
> >
> {: .solution}
{: .challenge}

# Input repeaters


The `each` qualifier allows you to repeat the execution of a process for each item in a list or channel, every time new data is received. For example:

~~~
transcriptome_ch = Channel.fromPath('data/yeast/transcriptome/*')
kmers = [12,31,45]

process index {
 echo true

 input:
  path transcriptome from transcriptome_ch
  each kmer from kmers  

 script:
 """
 echo salmon index -t $transcriptome -i index -k $kmer
 """

}
~~~
{: .language-groovy }

In the above example every time a file of sequences is received as input by the process, it executes three tasks running an alignment with a different value for the `kmer` option. This is useful when you need to repeat the same task for a given set of parameters.

> ## Exercise
> Extend the previous example so a task is executed for each read file matching the pattern data/yeast/reads/*.fq.gz and repeat the same task both with salmon and kallisto.
>> ## Solution
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
>>
> {: .solution}
{: .challenge}

## Outputs

The `output` declaration block allows us to define the channels used by the process to send out the results produced.

There can be defined at most one output block and it can contain one or more outputs declarations. The output block follows the syntax shown below:

~~~
output:
  <output qualifier> <output name> into <target channel>[,channel,..]
~~~  
{: .language-groovy }

### Output values

The `val` qualifier allows to output a value defined in the script context. In a common usage scenario, this is a value which has been defined in the *input* declaration block, as shown in the following example::


~~~
methods = ['salmon','kallisto']

process method_type {
  input:
  val x from methods

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
methods = ['salmon','kallisto']

process method_type {
  input:
  val x from methods

  output:
  path 'method.txt' into out_ch

  """
  echo $x > method.txt
  """
}

// use the view operator to display contents of the channel
out_ch.view({ "Received: $it" })
~~~

In the above example the process `method_type` creates a file named `method.txt` containing the method name.

Since a file parameter using the same name, `method.txt`, is declared in the output block , when the task is completed that file is sent over the `out_ch` channel. A downstream `operator` or `process` declaring the same channel as input will be able to receive it.

### Multiple output files

When an output file name contains a `*` or `?` wildcard character it is interpreted as a pattern match. This allows to capture multiple files into a list object and output them as a sole emission.

For example here we will capture sample.bam and sample.bam.bai in the output channel.

~~~

bam_ch = channel.fromPath("data/yeast/bams/*.bam").take(1)

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
> # Exercise
> Remove the flatMap operator and see out the output change. The documentation for the flatMap operator is available at this [link](https://www.nextflow.io/docs/latest/operator.html#flatmap).
> > solution
> >
> >
> {: .solution}
{: .challenge}

### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic evaluated string which references values defined in the input declaration block or in the script global context. For example::
fixme better example
~~~
species = ['human','mouse']

process align {
  input:
  val species_name from species

  output:
  file "${species_name}.aln" into genomes

  script:
  """
  echo ${species_name} > ${species_name}.aln
  """
}

genomes.view()
~~~
{: .language-groovy }

In the above example, each time the process is executed an alignment file is produced whose name depends on the actual value of the `species_name` input.

### Composite inputs and outputs

So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the `tuple` qualifiers.

When using channel emitting tuple the corresponding input declaration must be declared with a `tuple` qualifier followed by definition of each single element in the tuple.

In the same manner output channel emitting tuple of values can be declared using the `tuple` qualifier following by the definition of each tuple element in the tuple.


~~~
reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

process pseudo_align {
  input:
    tuple val(sample_id), path(reads) from reads_ch
  output:
    tuple val(sample_id), path('sample.bam') into bam_ch
  script:
  """
  salmon quant  -i $index \ -1 ${reads[0]} -2 ${reads[1]} -o $pair_id \
  --writeMappings |samtools sort |samtools view -bS -o sample.bam
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


The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending the state of various inputs and parameters. For example:
fixme
fixme better example> also there is an argument about if you should use when or filter the channel, see next episode on operators
may used autosome exmaple e.g. chr=channel.of(1..21,'X','Y')
~~~
//fix doesn't work
chr_ch = channel.of(1..21,'X','Y')
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process find {
  input:
  file fasta from proteins
  val chr from chr_ch

  when:
  chr % 2 == 0

  script:
  """
  echo blastp -query $fasta -db nr
  """
}
~~~
{: .language-groovy }


## Directives

Directive declarations allow the definition of optional settings that affect the execution of the current process without affecting the semantic of the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e. `input`, `output`, etc).

Directives are commonly used to define the amount of computing resources to be used or other meta directives like that allows the definition of extra information for configuration or logging purpose. For example:

~~~
chr_ch = channel.of(1..21,'X','Y')

process printchr {
  label 'big_mem'
  tag "$chr"
  echo true
  cpus 1
  memory 2.GB

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
> > chr_ch = channel.of(1..21,'X','Y')
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

## Organise outputs

### PublishDir directive

Nextflow manages independently workflow execution intermediate results from the pipeline expected outputs. Task output files are created in the task specific execution directory which is considered as a temporary directory that can be deleted upon completion.

The pipeline result files need to be marked explicitly using the directive [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) in the process that’s creating such file. For example:

~~~
chr_ch = channel.of(1..21,'X','Y')

process printchr {
  publishDir "results/chr"
  echo true


  input:
  val chr from chr_ch

  output:
  path "${chr}.txt"

  script:
  """
  sleep 2
  echo $chr > ${chr}.txt
  """
}

~~~
{: .language-groovy }

The above example will copy all `${chr}.txt` files created by the `printchr` process to the directory path `results/chr`.

###  Manage semantic sub-directories

You can use more then one publishDir to keep different outputs in separate directory. For example:




~~~
chr_ch = channel.of(1..21,'X','Y')

process printchr {
  publishDir "results/chr/autosomes",pattern:"[0-9]*.txt"
  publishDir "results/chr/sex",pattern:"{X,Y}*.txt"
  echo true


  input:
  val chr from chr_ch

  output:
  path "${chr}.txt"

  script:
  """
  sleep 2
  echo $chr > ${chr}.txt
  """
}
~~~
{: .language-groovy }

The above example will create an output structure in the directory my-results, which contains a separate sub-directory for each given sample ID each of which contain the folders counts and outlooks.


fixme Create an exercise



> ## Nextflow Patterns
> The [Nextflow Patterns page](http://nextflow-io.github.io/patterns/) collects some recurrent implementation patterns used in Nextflow applications.
{: .callout}

{% include links.md %}
