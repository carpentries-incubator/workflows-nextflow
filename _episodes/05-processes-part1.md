---
title: "Processes"
teaching: 60
exercises: 20
questions:
- "How do I run tasks/processes in Nextflow?"
- "How do I get data, files and values, into a processes?"
objectives:
- "Understand how Nextflow uses processes to execute tasks."
- "Create a Nextflow process."
- "Define inputs to a process."
keypoints:
- "A Nextflow process is an independent task/step in a workflow"
- "Processes contain up to five definition blocks including, directives, inputs, outputs, when clause and finally a script block."
- "The script block contains the commands you would like to run."
- "Inputs are defined using the input blocks."
---


# Processes

We now know how to create and use Channels to send data around a workflow. We will now see how to run tasks within a workflow using processes.

A `process` is the way Nextflow execute commands you would run on the command line or custom scripts.

Processes can be thought of as a particular task/steps in a workflow, e.g. an alignment step in RNA-Seq analysis. Processes  are independent of each other (don't require another processes to execute) and can not communicate/write to each other . It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output channels.

For example, below is a command you would run to create a index directory for the [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) aligner on the command line:

~~~
$ salmon index -t data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i data/yeast/salmon_index --kmerLen 31
~~~
{: .language-bash }

Below we will show how to convert this into a simple Nextflow process.

## Process definition

The process definition starts with keyword the `process`, followed by process name and finally the process `body` delimited by curly brackets `{}`. The process body must contain a string  which represents the command or, more generally, a script that is executed by it.

~~~
process INDEX {
  script:
  "salmon index -t ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i data/yeast/salmon_index --kmerLen 31"
}

~~~
{: .language-groovy }

This process would run once.

> ## projectDir
We use the special Nextflow variable `${projectDir}` to specify the directory where the main script is located. This is important as Nextflow scripts are executed in a separate working directory.
 {: .callout }

To add the process to a workflow add a `workflow` block below the process.
We will learn more about the `workflow` block in the next episode.

**Note:** As we are using DSL2 we need to include `nextflow.enable.dsl=2` in the script.

~~~
//process_index.nf
nextflow.enable.dsl=2

process INDEX {
  script:
  "salmon index -t ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i data/yeast/salmon_index --kmerLen 31"
}

workflow {
  //process is called like a function in the workflow block
  INDEX()
}
~~~
{: .language-groovy }

We can now run the process

~~~
$ nextflow run index.nf
~~~
{: .language-bash }


> ## A Simple Process
>
> Create a Nextflow script `simple_process.nf` that has one process `SALMON_VERSION` that runs the command.
> ~~~
> salmon --version
> ~~~
> {: .language-bash}
>
> > ## Solution
> > ~~~
> > nextflow.enable.dsl=2
> >
> > process SALMON_VERSION {
> >    
> >   script:
> >   """
> >    salmon --version
> >    """
> > }
> >
> > workflow {
> >   SALMON_VERSION()
> > }
> > ~~~
> > {: .language-groovy }
> > **Note** We need to add the Nextflow run option `-process.echo` to print the output to the terminal.
> > ~~~
> > $ nextflow run simple_process.nf -process.echo
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}


### Definition blocks

The previous example was a simple `process` with no defined inputs and outputs that ran only once. To control inputs, outputs and how a command is executed a process may contain five definition blocks:

1. **directives**: allow the definition of optional settings that affect the execution of the current process e.g. the number of cpus a task uses and the amount of memory allocated.
1. **inputs**: Define the input dependencies, usually channels, which determines the number of times a process is executed.
1. **outputs**: Defines the output channels used by the process to send results/data produced by the process.
1. **when clause**: Allow you to define a condition that must be verified in order to execute the process.
1. **The script block**: The script block is a statement within quotes that defines the command that is executed by the process to carry out its task.


The syntax is defined as follows:

~~~
process < NAME > {
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

The `script` block is a String "statement" that defines the command that is executed by the process to carry out its task. These are normally the commands you would run on a terminal.

A process contains only one `script` block, and it must be the last statement when the process contains `input` and `output` declarations.

The `script` block can be a simple one line string in quotes e.g. `"samtools index ref1.sorted.bam"` .
Or, for commands that span multiple lines you can encase the command in  triple quotes `"""`.

For example:

~~~
//process_multi_line.nf
nextflow.enable.dsl=2

process PROCESSBAM {
    script:
    """
    samtools sort -o ref1.sorted.bam ${projectDir}/data/yeast/bams/ref1.bam
    samtools index ref1.sorted.bam
    samtools flagstat ref1.sorted.bam
    """
}

workflow {
  PROCESSBAM()
}
~~~
{: .language-groovy }

By default the process command is interpreted as a **Bash** script. However any other scripting language can be used just simply starting the script with the corresponding [Shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) declaration. For example:

~~~
//process_python.nf
nextflow.enable.dsl=2

process PYSTUFF {
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

workflow {
  PYSTUFF()
}
~~~
{: .language-groovy }

~~~
//process_rscript.nf
nextflow.enable.dsl=2

process RSTUFF {
  script:
  """
  #!/usr/bin/env Rscript
  library("ShortRead")
  countFastq(dirPath="data/yeast/reads/ref1_1.fq.gz")
  """
}

workflow {
  RSTUFF()
}
~~~
{: .language-groovy }

This allows the the use of a different programming languages which may better fit a particular job. However, for large chunks of code is suggested to save them into separate files and invoke them from the process script.

~~~
nextflow.enable.dsl=2

process PYSTUFF {

  script:
  """
  python myscript.py
  """
}

workflow {
  PYSTUFF()
}
~~~
{: .language-groovy }

### Script parameters

The command in the `script` block can be defined dynamically using Nextflow variables e.g. `${projectDir}`.
To reference a variable in the script block you can use the `$` in front of the Nextflow variable name, and additionally you can add `{}` around the variable name e.g. `${projectDir}`.

> ##  Variable substitutions
> Similar to bash scripting Nextflow uses the "$" character to introduces variable substitutions. The variable name to be expanded may be enclosed in braces `{variable_name}`, which are optional but serve to protect the variable to be expanded from characters immediately following it which could be interpreted as part of the name.
{: .callout }

In the example below the variable `kmer` is set to the value 31 at the top of the Nextflow script.
The variable is referenced using the `$kmer` syntax within the multi-line string statement in the `script` block.
A Nextflow variable can be used multiple times in the script block.

~~~
//process_script.nf
nextflow.enable.dsl=2

kmer = 31

process INDEX {

  script:
  """
  salmon index \
  -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  -i index \
  --kmer $kmer
  echo "kmer size is $kmer"
  """
}

workflow {
  INDEX()
}
~~~
{: .language-groovy }

In most cases we do not want to hard code parameter values. We saw in the parameter episode the use of a special Nextflow  variable `params` that can be used to assign values from the command line. You would do this by adding a key name to the params variable and specifying a value, like `params.keyname = value`

In the example below we define the variable `params.kmer` with a default value of 31 in the Nextflow script.
~~~
//process_script_params.nf
nextflow.enable.dsl=2

params.kmer = 31

process INDEX {

  script:
  """
  salmon index \
  -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  -i index \
  --kmer $params.kmer
  echo "kmer size is $params.kmer"
  """
}

workflow {
  INDEX()
}
~~~
{: .language-groovy }

Remember, we can change the default value of `kmer` to 11 by running the Nextflow script using the command below. **Note:** parameter have two hyphens `--` .

~~~
nextflow run process_script_params.nf --kmer 11
~~~
{: .language-bash }


> ## Script parameters
>
> For the Nextflow script below.
> ~~~
> //process_script_params.nf
> nextflow.enable.dsl=2
> params.kmer = 31
>
> process INDEX {
>
>  script:
>  """
>  salmon index -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer $params.kmer
>  echo "kmer size is" $params.kmer
>  """
> }
>
> workflow {
>   INDEX()
> }
> ~~~
> {: .language-groovy}
>
> Change the k-mer value used to create the salmon index on the command line using the `--kmer` command line option.
>
> ~~~
> $ nextflow run process_script_params.nf --kmer <some value> -process.echo
> ~~~
> {: .language-bash}
> **Note:** The kmer value must not be greater than 31 and an odd number.
> **Note:** The Nextflow option `-process.echo` will print the process' stdout to the terminal.
>
> > ## Solution
> > ~~~
> > nextflow run process_script_params.nf --kmer 27 -process.echo
> > ~~~
> > {: .language-bash }
> ~~~
> > kmer size is 27
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}


### Bash variables

Nextflow uses the same Bash syntax for variable substitutions, `$variable`, in strings.
However, Bash variables need to be escaped using `\` character in front of `\$variable` name.

In the example below we will use set the bash variable KMERSIZE to the value of `$params.kmer`.


~~~
//process_escape_bash.nf
nextflow.enable.dsl=2

process INDEX {

  script:
  """
  #set bash variable KMERSIZE
  KMERSIZE=$params.kmer
  salmon index -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer \$KMERSIZE
  echo "kmer size is $params.kmer"
  """
}

params.kmer = 31

workflow {
  INDEX()
}
~~~
{: .language-groovy }

### Shell

Another alternative is to use a `shell` block definition instead of `script`.
When using the `shell` statement Bash variables are referenced in the normal way `$my_bash_variable`;
However, the `shell` statement uses a different syntax for Nextflow variable substitutions: `!{nextflow_variable}`, which is needed to use both Nextflow and Bash variables in the same script.

For example in the script below that uses the `shell statment `
we reference the Nextflow variables as such, `!{projectDir}` and `!{params.kmer}` and the Bash variable
as `$PWD`.

```
//process_shell.nf
nextflow.enable.dsl=2

params.kmer = 31

process INDEX {

  shell:
  '''
  #set bash variable KMERSIZE
  KMERSIZE=!{params.kmer}
  salmon index -t !{projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer ${KMERSIZE}
  echo "kmer size is  !{params.kmer}"
  '''
}

workflow {
  INDEX()
}
```
{: .language-groovy }


### Conditional script execution

Sometimes you want to change how a process is run depending on some condition. In Nextflow scripts we can use conditional statements such as the `if` statement or any other expression evaluating to boolean value `true` or `false`.

### If statement

The `if` statement uses the same syntax common other programming lang such Java, C, JavaScript, etc.

~~~
if( < boolean expression > ) {
    // true branch
}
else if ( < boolean expression > ) {
    // true branch
}
else {
    // false branch
}
~~~
{: .language-groovy }


For example, the Nextflow script below will use the `if` statement to change what index is created depending on the Nextflow variable `params.aligner`.

~~~
//process_conditional.nf
nextflow.enable.dsl=2

params.aligner = 'kallisto'
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.kmer = 31

process INDEX {
  script:
  if( params.aligner == 'kallisto' ) {
    """
    echo indexed using kallisto
    kallisto index -i index  -k $params.kmer $params.transcriptome
    """
  }  
  else if( params.aligner == 'salmon' ) {
    """
    echo indexed using salmon
    salmon index -t $params.transcriptome -i index --kmer $params.kmer
    """
  }  
  else {
    """
    echo Unknown aligner $params.aligner"
    """
  }  
}

workflow {
  INDEX()
}
~~~
{: .language-groovy }

~~~
nextflow run main.nf -process.echo --aligner kallisto
~~~
{: .language-bash }

~~~
indexed using kallisto
~~~
{: .language-groovy }

## Inputs

Processes are isolated from each other but can communicate by sending values and files via Nextflow channels into `input` and `output` blocks.

The `input` block defines which channels the process is expecting to receive input from.
The number of elements in input channels determine the process dependencies and the number of time a process executes.

![Process Flow](../fig/channel-process.png)


You can only define one input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below:

~~~
input:
  <input qualifier> <input name>
~~~
{: .language-groovy }

The input qualifier declares the type of data to be received.

> ## Input qualifiers
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
//process_input_value.nf
nextflow.enable.dsl=2

process PRINTCHR {

  input:
  val chr

  script:
  """
  echo processing chromosome $chr
  """
}

chr_ch = Channel.of( 1..22,'X','Y' )

workflow {

  PRINTCHR(chr_ch)
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

The input file name can be defined dynamically by defining the input name as a Nextflow variable and referenced in the script using  `$variable_name` syntax.

For example in the script below we assign the variable name `read` to the input files and it
is referenced using the variable substitution syntax `${read}` in the script block:

~~~
//process_input_file.nf
nextflow.enable.dsl=2

process NUMLINES {
    input:
    path read
    script:
    """
    printf '${read} '
    gunzip -c ${read} | wc -l
    """
}

reads_ch = Channel.fromPath( 'data/yeast/reads/ref*.fq.gz')

workflow {
  NUMLINES(reads_ch)
}

~~~
{: .language-groovy }

~~~
[cd/77af6d] process > NUMLINES (1) [100%] 6 of 6 ✔
ref1_1.fq.gz 58708

ref3_2.fq.gz 52592

ref2_2.fq.gz 81720

ref2_1.fq.gz 81720

ref3_1.fq.gz 52592

ref1_2.fq.gz 58708
~~~
{: .output }

The input name can also be defined as user specified filename inside quotes as in the example below where specific the name of the file as `'sample.fq.gz'`.

~~~
//process_input_file_02.nf
nextflow.enable.dsl=2

process NUMLINES {
    input:
    path 'sample.fq.gz'
    script:
    """
    printf 'sample.fq.gz'
    gunzip -c sample.fq.gz | wc -l
    """
}

reads_ch = Channel.fromPath( 'data/yeast/reads/ref*.fq.gz')

workflow {
  NUMLINES(reads_ch)
}

~~~
{: .language-groovy }


~~~
[d2/eb0e9d] process > NUMLINES (1) [100%] 6 of 6 ✔
sample.fq.gz58708

sample.fq.gz52592

sample.fq.gz81720

sample.fq.gz81720

sample.fq.gz52592

sample.fq.gz58708
~~~
{: .output }


> ## File Objects as inputs
> When a process declares an input file the corresponding channel elements must be file objects i.e. created with the path helper function from the file specific channel factories e.g. `Channel.fromPath` or `Channel.fromFilePairs`.
{: .callout}


> ## Add input channel
> Add an input channel to the script below that takes the reads channel as input.
> [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality control tool for high throughput sequence data.
> ~~~
> //process_exercise_input.nf
> nextflow.enable.dsl=2
>
> process FASTQC {
>    //add input channel
>    
>    script:
>    """
>    mkdir fastqc_out
>    fastqc -o fastqc_out ${reads}
>    """
> }
> reads_ch = Channel.fromPath( 'data/yeast/reads/*_1.fq.gz' )
>
> workflow {
>   FASTQC(reads_ch)
> }
> ~~~
> {: .language-groovy }
> Then run your script using
> ~~~
> nextflow run fastqc.nf
> ~~~
> {: .language-bash }
> > ## Solution
> > ~~~
> > //process_exercise_input_answer.nf
> > nextflow.enable.dsl=2
> > process FASTQC {
> >    input:
> >    path reads
> >    script:
> >    """
> >    mkdir fastqc_out
> >    fastqc -o fastqc_out ${reads}
> >    """
> > }
> > reads_ch = Channel.fromPath( 'data/yeast/reads/*_1.fq.gz' )
> >
> > workflow {
> >   FASTQC(reads_ch)
> > }
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

### Combining input channels

A key feature of processes is the ability to handle inputs from multiple channels.
However it’s important to understands how the content of channel and affect the execution of a process.

Consider the following example:

~~~
//process_combine.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  val y
  script:
   """
   echo $x and $y
   """
}

num_ch = Channel.of(1,2,3)
letters_ch = Channel.of('a','b','c')

workflow {
  COMBINE(num_ch,letters_ch)
}
~~~
{: .language-groovy }

Both channels contain three elements, therefore the process is executed three times, each time with a different pair:

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
//process_combine_02.nf

nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  val y
  script:
   """
   echo $x and $y
   """
}

ch_num = Channel.of(1,2)
ch_letters = Channel.of('a','b','c','d')

workflow {
  COMBINE(ch_num,ch_letters)
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
//process_combine_03.nf

nextflow.enable.dsl=2

process COMBINE {
  echo true
  input:
  val x
  val y
  script:
   """
   echo $x and $y
   """
}
ch_num = Channel.value(1)
ch_letters = Channel.of('a','b','c')

workflow {
  COMBINE(ch_num,ch_letters)
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


> ##  Combining input channels
> Write a nextflow script `salmon_index_combine.nf` that combines two input channels
> ~~~
> transcriptome_ch = channel.value('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
> kmer_ch = channel.of(21)
> ~~~
> {: .language-groovy }
And include the command below in the script directive
>
> ~~~~
  script:
   """
   salmon index -t $transcriptome -i index -k $kmer` .
   """
> ~~~~
> {: .language-groovy }
> > ## Solution
> > ~~~
> > // process_exercise_combine_answer.nf
> > nextflow.enable.dsl=2
> > process COMBINE {
> >  input:
> >  path transcriptome
> >  val kmer
> >  script:
> >   """
> >   salmon index -t $transcriptome -i index -k $kmer
> >   """
> > }
> >
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.of(21)
> >
> > workflow {
> >   COMBINE(transcriptome_ch,kmer_ch)
> > }
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

# Input repeaters

We saw previously that by default the number of time a process run is defined by the queue channel with the fewest items. However, the `each` qualifier allows you to repeat the execution of a process for each item in a list or a queue channel, every time new data is received.

For example if we can fix the previous example by using the input qualifer `each` for the letters queue channel:

~~~
//process_repeat.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  each y
  script:
   """
   echo $x and $y
   """
}

ch_num = Channel.of(1,2)
ch_letters = Channel.of('a','b','c','d')

workflow {
  COMBINE(ch_num,ch_letters)
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

> ## Input repeaters
>  ~~~
>  nextflow.enable.dsl=2
>  process COMBINE {
>   input:
>   path transcriptome
>   val kmer
>   script:
>    """
>    salmon index -t $transcriptome -i index -k $kmer
>    """
>  }
>
>  transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
>  kmer_ch = channel.of(21)
>
>  workflow {
>    COMBINE(transcriptome_ch,kmer_ch)
>  }
>  ~~~
>  {: .language-groovy }
> Extend the previous  exercise script `salmon_index.nf` above, by adding more values to the `kmer` queue channel
> ~~~  
> kmer_ch = channel.of(21,27,31)
> ~~~
> {: .language-groovy }
>
> and changing the `transcriptome` input qualifer from `path` to `each`.
> How many times does this process run ?
>
> > ## Solution
> > ~~~
> > //process_exercise_repeat_answer.nf
> > nextflow.enable.dsl=2
> >
> > process COMBINE {
> >  input:
> >  each transcriptome from transcriptome_ch
> >  path kmer from kmer_ch
> >  script:
> >   """
> >   echo salmon index -t $transcriptome -i index -k $kmer
> >   """
> > }
> >
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz',checkIfExists: true)
> > kmer_ch = channel.of(21,27,31)
> >
> >  workflow {
> >   COMBINE(transcriptome_ch,kmer_ch)
> > }
> > ~~~
> > {: .language-groovy }
> > This process runs three times.
> {: .solution}
{: .challenge}


{: .output }
{% include links.md %}