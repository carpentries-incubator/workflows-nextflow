---
title: "Processes Part 2"
teaching: 30
exercises: 10
questions:
- "How do I get data, files, and values,  out of processes?"
- "How do I handle grouped input and output?"
- "How can I control when a process is executed?"
- "How do I control resources, such as number of CPUs and memory, available to processes?"
- "How do I save output/results from a process?"
objectives:
- "Define outputs to a process."
- "Understand how to handle grouped input and output using the tuple qualifier."
- "Understand how to use conditionals to control process execution."
- "Use process directives to control execution of a process."
- "Use the `publishDir` directive to save result files to a directory."
keypoints:
- "Outputs to a process are defined using the output blocks."
- "You can group input and output data from a process using the tuple qualifer."
- "The execution of a process can be controlled using the `when` declaration and conditional statements."
- "Files produced within a process and defined as `output` can be saved to a directory using the `publishDir` directive."
---

## Outputs

We have seen how to input data into a process; now we will see how to output files and values from a process.

The `output` declaration block allows us to define the channels used by the process to send out the files and values produced.

An output block is not required, but if it is present it can contain one or more output declarations.

The output block follows the syntax shown below:

~~~
output:
  <output qualifier> <output name>
  <output qualifier> <output name>
  ...
~~~  
{: .language-groovy }


### Output values

Like the input, the type of output data is defined using type qualifiers.

The `val` qualifier allows us to output a value defined in the script.


Because Nextflow processes can only communicate through channels, if we want to share a value input into one process as input to another process we would need to define that value in the output declaration block as shown in the following example:


~~~
//process_output_value.nf
nextflow.enable.dsl=2

process METHOD {
  input:
  val x

  output:
  val x

  script:
  """
  echo $x > method.txt
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
methods_ch = channel.of('salmon', 'kallisto')

workflow {
  METHOD(methods_ch)
  // use the view operator to display contents of the channel
  METHOD.out.view({ "Received: $it" })
}
~~~
{: .language-groovy }

~~~
[d8/fd3f9d] process > METHOD (2) [100%] 2 of 2 ✔

Received: salmon
Received: kallisto
~~~
{: .output }

### Output files

If we want to capture a file instead of a value as output we can use the
`path` qualifier that can capture one or more files produced by the process, over the specified channel.

~~~
//process_output_file.nf
nextflow.enable.dsl=2

methods_ch = channel.of('salmon', 'kallisto')

process METHOD {
  input:
  val x

  output:
  path 'method.txt'

  """
  echo $x > method.txt
  """
}

workflow {
  METHOD(methods_ch)
  // use the view operator to display contents of the channel
  METHOD.out.view({ "Received: $it" })
}
~~~
{: .language-groovy }

~~~
executor >  local (2)
[13/85ecb6] process > METHOD (1) [100%] 2 of 2 ✔

Received: /Users/ggrimes2/Downloads/nf-training/work/ea/f8998abfdcbfc6038757a60806c00a/method.txt
Received: /Users/ggrimes2/Downloads/nf-training/work/13/85ecb69a9bc85370e749254646984d/method.txt
~~~
{: .output }

In the above example the process `METHOD` creates a file named `method.txt` in the work directory containing the method name.

Since a file parameter using the same name, `method.txt`, is declared in the output block, when the task is completed that file is sent over the output channel.

A downstream `operator`, such as `.view` or a `process` declaring the same channel as input will be able to receive it.

### Multiple output files

When an output file name contains a `*` or `?` metacharacter it is interpreted as a pattern match.
This allows us to capture multiple files into a list and output them as a one item channel.

For example, here we will capture the files `fastqc.html` and directory `fastqc.zip` produced as results from FastQC in the output channel.

~~~
//process_output_multiple.nf
nextflow.enable.dsl=2

process FASTQC {
  input:
  path read

  output:
  path "fqc_res/*"

  script:
  """
  mkdir fqc_res
  fastqc $read -o fqc_res
  """
}

read_ch = channel.fromPath("data/yeast/reads/ref1*.fq.gz")

workflow {
  FASTQC(read_ch)
  FASTQC.out.view()
}
~~~
{: .language-groovy }

~~~
$ nextflow run process_output_multiple.nf
~~~
{: .language-bash }

~~~
[3e/86de98] process > FASTQC (2) [100%] 2 of 2 ✔
[work/64/9de164199568800a4609c3af78cf71/fqc_res/ref1_2_fastqc.html, work/64/9de164199568800a4609c3af78cf71/fqc_res/ref1_2_fastqc.zip]
Analysis complete for ref1_2.fq.gz

[work/3e/86de9868ecf321702f2df0b8ccbfd3/fqc_res/ref1_1_fastqc.html, work/3e/86de9868ecf321702f2df0b8ccbfd3/fqc_res/ref1_1_fastqc.zip]
Analysis complete for ref1_1.fq.gz
~~~
{: .output}

**Note:** There are some caveats on glob pattern behaviour:

* Input files are not included in the list of possible matches.
* Glob pattern matches against both files and directories path.
* When a two stars pattern `**` is used to recurse through subdirectories, only file paths are matched i.e. directories are not included in the result list.


> ## Output channels
> Modify the nextflow script `process_exercise_output.nf` to include an output block that captures the different index folders `index_$kmer`.
> Use the `view` operator on the output channel.
> ~~~
> //process_exercise_output.nf
> nextflow.enable.dsl=2
>
> process INDEX {
>   input:
>   path transcriptome
>   each kmer
> 
>   //add output block here to capture index folders "index_$kmer"
> 
>   script:
>   """
>   salmon index -t $transcriptome -i index_$kmer -k $kmer
>   """
> }
> 
> transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
> kmer_ch = channel.of(21, 27, 31)
>
> workflow {
>   INDEX(transcriptome_ch, kmer_ch)
>   INDEX.out.view()
> }
> ~~~
> {: .language-groovy }
> > ## Solution
> > ~~~
> > //process_exercise_output_answer.nf
> > nextflow.enable.dsl=2
> >
> > process INDEX {
> >   input:
> >   path transcriptome
> >   each kmer
> >  
> >   output:
> >   path "index_${kmer}"
> >  
> >   script:
> >   """
> >   salmon index -t $transcriptome -i index_$kmer -k $kmer
> >   """
> > }
> > transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
> > kmer_ch = channel.of(21, 27, 31)
> >
> > workflow {
> >   INDEX(transcriptome_ch, kmer_ch)
> >   INDEX.out.view()
> > }
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

### Grouped inputs and outputs

So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the `tuple` qualifiers.

In tuples the first item is the grouping key and the second item is the list.

~~~
[group_key,[file1,file2,...]]
~~~


When using channel containing a tuple, such a one created with `.filesFromPairs` factory method, the corresponding input declaration must be declared with a `tuple` qualifier, followed by definition of each item in the tuple.


~~~
//process_tuple_input.nf
nextflow.enable.dsl=2

process TUPLEINPUT{
  input:
  tuple val(sample_id), path(reads)
  
  script:
  """
  echo $sample_id
  echo $reads
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  TUPLEINPUT(reads_ch)
}

~~~
{: .language-groovy }

outputs

~~~
ref1
ref1_1.fq.gz ref1_2.fq.gz
~~~
{: .output }

In the same manner an output channel containing tuple of values can be declared using the `tuple` qualifier following by the definition of each tuple element in the tuple.

In the code snippet below the output channel would contain a tuple with the grouping key value as the Nextflow variable `sample_id` and a list containing the files matching the following pattern `"*FP*.fq.gz"`.


~~~
output:
  tuple val(sample_id), path("*FP*.fq.gz")
~~~
{: .language-groovy }

An example can be seen in this script below.


~~~
//process_tuple_io_fastp.nf
nextflow.enable.dsl=2

process FASTP {
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("*FP*.fq.gz")
  
  script:
  """
  fastp \
   -i ${reads[0]} \
   -I ${reads[1]} \
   -o ${sample_id}_FP_R1.fq.gz \
   -O ${sample_id}_FP_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  FASTP(reads_ch)
  FASTP.out.view()
}

~~~
{: .language-groovy }

~~~
nextflow run process_tuple_io_fastp.nf
~~~
{: .language-bash }

The output is now a tuple containing the sample id and the two processed fastq files.

~~~
[ref1, [work/9b/fcca75db83718a15f7a95caabfbc15/ref1_FP_R1.fq.gz, work/9b/fcca75db83718a15f7a95caabfbc15/ref1_FP_R2.fq.gz]]
~~~
{: .output }

> ## Composite inputs and outputs
> Fill in the blank ___ input and output qualifers for `process_exercise_tuple.nf`.
> **Note:** the output for the FASTQC process is a directory named `fastq_out`.
> ~~~
> //process_exercise_tuple.nf
> nextflow.enable.dsl=2
>
> process FASTQC {
>   input:
>   tuple ___(sample_id), ___(reads)
> 
>   output:
>   tuple ___(sample_id), ___("fastqc_out")
> 
>   script:
>   """
>   mkdir fastqc_out
>   fastqc $reads -o fastqc_out -t 1
>   """
> }
> 
> reads_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')
>
> workflow{
>   FASTQC(reads_ch)
>   FASTQC.out.view()
> }
> ~~~
> {: .language-groovy }
> > ## Solution
> > ~~~
> > //process_exercise_tuple_answer.nf
> > nextflow.enable.dsl=2
> >
> >process FASTQC {
> >   input:
> >   tuple val(sample_id), path(reads)
> >   
> >   output:
> >   tuple val(sample_id), path("fastqc_out")
> >  
> >   script:
> >   """
> >   mkdir fastqc_out
> >   fastqc $reads -o fastqc_out -t 1
> >   """
> > }
> > 
> > reads_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')
> > 
> > workflow{
> >   FASTQC(reads_ch)
> >   FASTQC.out.view()
> >}
> > ~~~
> > {: .language-groovy }
> > ~~~
> > N E X T F L O W  ~  version 21.04.0
> > Launching `process_exercise_tuple.nf` [spontaneous_coulomb] - revision: 06ff22f1a9
> > executor >  local (3)
> > [75/f4a44d] process > FASTQC (3) [100%] 3 of 3 ✔
> > [ref3, work/99/a7d9176e332fdc0988973dbb89df63/fastqc_out]
> > [ref2, /work/53/e3cbd39afa9f0f84a3d9cd060d991a/fastqc_out]
> > [ref1, work/75/f4a44d0bc761fa4774c2f23a465766/fastqc_out]
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

## Conditional execution of a process

The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value; `true` or `false`.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters.

In the example below the process `CONDITIONAL` will only execute when the value of the `chr` variable is less than or equal to 5:

~~~
//process_when.nf
nextflow.enable.dsl=2

process CONDITIONAL {
  input:
  val chr

  when:
  chr <= 5

  script:
  """
  echo $chr
  """
}

chr_ch = channel.of(1..22)

workflow {
  CONDITIONAL(chr_ch)
}
~~~
{: .language-groovy }

~~~
4

5

2

3

1
~~~
{: .output }

## Directives

Directive declarations allow the definition of optional settings, like the number of `cpus` and amount of `memory`, that affect the execution of the current process without affecting the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e. `input`, `output`, etc).


**Note:** You do not use `=` when assigning a value to a directive.


Directives are commonly used to define the amount of computing resources to be used or extra information for configuration or logging purpose.

For example:

~~~
//process_directive.nf
nextflow.enable.dsl=2

process PRINTCHR {
  tag "tagging with chr$chr"
  cpus 1
  echo true

  input:
  val chr

  script:
  """
  echo processing chromosome: $chr
  echo number of cpus $task.cpus
  """
}

chr_ch = channel.of(1..22, 'X', 'Y')

workflow {
  PRINTCHR(chr_ch)
}
~~~
{: .language-groovy }

~~~
processing chromosome: 1
number of cpus 1

processing chromosome: 2
number of cpus 1

processing chromosome: 6
number of cpus 1
[..truncated..]
~~~
{: .output }

The above process uses the three directives, `tag`, `cpus` and `echo`.

The `tag` directive to allow you to give a custom tag to each process execution. This tag makes it easier to identify a particular task (executed instance of a process) in a log file or in the execution report.

The second directive `cpus`  allows you to define the number of CPUs required for each task.

The third directive `echo true` prints the stdout to the terminal.

We use the Nextflow `task.cpus` variable to capture the number of cpus assigned to a task. This is frequently used to specify the number of threads in a multi-threaded command in the script block.

Another commonly used directive is memory specification: `memory`.

A complete list of directives is available at this [link](https://www.nextflow.io/docs/latest/process.html#directives).

> ## Adding directives
> Modify the Nextflow script `process_exercise_directives.nf`
>
> 1. Add a `tag` directive logging the sample_id in the execution output.
> 1. Add a `cpus` directive to specify the number of cpus as 2.
> 1. Change the fastqc `-t` option value to `$task.cpus` in the script directive.
>
> ~~~
> //process_exercise_directives.nf
> nextflow.enable.dsl=2
>
> process FASTQC {
>   //add tag directive
>   //add cpu directive
>  
>   input:
>   tuple val(sample_id), path(reads)
>   
>   output:
>   tuple val(sample_id), path("fastqc_out")
>   
>   script:
>   """
>   mkdir fastqc_out
>   fastqc $reads -o fastqc_out -t 1
>   """
> }
>
> read_pairs_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')
> 
> workflow {
>   FASTQC(read_pairs_ch)
>   FASTQC.out.view()
> }
>  ~~~
>  {: .language-groovy }
> > ## solution
> > ~~~
> > //process_directives_answer.nf
> > nextflow.enable.dsl=2
> >
> > process FASTQC {
> >   tag "$sample_id"
> >   cpus 2
> >   
> >   input:
> >   tuple val(sample_id), path(reads)
> >   
> >   output:
> >   tuple val(sample_id), path("fastqc_out")
> >  
> >   script:
> >   """
> >   mkdir fastqc_out
> >   fastqc $reads -o fastqc_out -t 1
> >   """
> > }
> >
> > read_pairs_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')
> >
> > workflow {
> >   FASTQC(read_pairs_ch)
> >   FASTQC.out.view()
> > }
> >  ~~~
> > {: .language-groovy }
> > ~~~
> > N E X T F L O W  ~  version 21.04.0
> > Launching `process_exercise_directives.nf` [sad_rosalind] - revision: 2ccbfa4937
> > executor >  local (3)
> > [90/de1125] process > FASTQC (ref1) [100%] 3 of 3 ✔
> > [ref2, work/ea/9e6a341b88caf8879e8d18b77049c8/fastqc_out]
> > [ref3, work/94/d059b816a9ec3d868f2924c26813e7/fastqc_out]
> > [ref1, work/90/de11251d362f494d6650789d9f8c1d/fastqc_out]
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

## Organising outputs

### PublishDir directive

Nextflow manages intermediate results from the pipeline's expected outputs independently.

Files created by a `process` are stored in a task specific working directory which is considered as temporary. Normally this is under the `work` directory, which can be deleted upon completion.

The files you want the workflow to return as results need to be defined in the `process` `output` block and then the output directory specified using the `directive` `publishDir`. More information [here](https://www.nextflow.io/docs/latest/process.html#publishdir).

**Note:** A common mistake is to specify an output directory in the `publishDir` directive while forgetting to specify the files you want to include in the `output` block.

~~~
publishDir <directory>, parameter: value, parameter2: value ...
~~~

For example if we want to capture the results of the `QUANT` process in a `results/quant` output directory we
need to define the files in the `output` and  specify the location of the results directory in the `publishDir` directive:

~~~
//process_publishDir.nf
nextflow.enable.dsl=2

process QUANT {
  publishDir "results/quant"
  
  input:
  tuple val(sample_id), path(reads)
  each index
  
  output:
  tuple val(sample_id), path("${sample_id}_salmon_output")
  
  script:
  """
  salmon quant -i $index \
   -l A \
   -1 ${reads[0]} \
   -2 ${reads[1]} \
   -o ${sample_id}_salmon_output
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

workflow {
  QUANT(reads_ch, index_ch)
  QUANT.out.view()
}

~~~
{: .language-groovy }

~~~
$ nextflow run process_publishDir.nf
~~~
{: .language-bash }
~~~
N E X T F L O W  ~  version 21.04.0
Launching `process_publishDir.nf` [friendly_pauling] - revision: 9b5c315893
executor >  local (1)

[48/f97234] process > QUANT (1) [100%] 1 of 1 ✔
[ref1, work/48/f97234d7185cbfbd86e2f11c1afab5/ref1_salmon_output]
~~~
{: .output }

We can use the UNIX command `tree` to examine the contents of the results directory.

~~~
tree results
~~~
{: .language-bash }

~~~
results/
└── quant
    └── ref1_salmon_output -> work/48/f97234d7185cbfbd86e2f11c1afab5/ref1_salmon_output
~~~
{: .output }


In the above example, the `publishDir "results/quant"`,  creates a symbolic link `->` to the output files specified by the process `salmon_quant` to the directory path `results/quant`.

> ## publishDir
> The publishDir output is relative to the path the pipeline run has been launched. Hence, it is a good practice to use [implicit variables](https://www.nextflow.io/docs/latest/script.html?highlight=projectdir#script-implicit-variables) like `projectDir` to specify publishDir value.
{: .callout }

### publishDir parameters

The `publishDir` directive can take optional parameters, for example the `mode` parameter can take the value `"copy"` to specify that you wish to copy the file to output directory rather than just a symbolic link to the files in the working directory. Since the working directory is generally deleted on completion of a pipeline, it is safest to use `mode: "copy"` for results files. The default mode (symlink) is helpful for checking intermediate files which are not needed in the long term.
~~~
publishDir "results/quant", mode: "copy"
~~~
{: .language-groovy }

Full list [here](https://www.nextflow.io/docs/latest/process.html#publishdir).


###  Manage semantic sub-directories

You can use more than one `publishDir` to keep different outputs in separate directories. To specify which files to put in which output directory use the parameter `pattern` with the a glob pattern that selects which files to publish from the overall set of output files.

In the example below we will create an output folder structure in the directory results, which contains a separate sub-directory for bam files, `pattern:"*.bam"` ,  and a salmon output directory, `${sample_id}_salmon_output"`. Remember, we need to specify the files we want to copy as outputs.

~~~
//process_publishDir_semantic.nf
nextflow.enable.dsl=2

process QUANT {
  publishDir "results/bams", pattern: "*.bam", mode: "copy"
  publishDir "results/quant", pattern: "${sample_id}_salmon_output", mode: "copy"

  input:
  tuple val(sample_id), path(reads)
  path index
  
  output:
  tuple val(sample_id), path("${sample_id}.bam")
  path "${sample_id}_salmon_output"
  
  script:
  """
  salmon quant -i $index \
   -l A \
   -1 ${reads[0]} \
   -2 ${reads[1]} \
   -o ${sample_id}_salmon_output \
   --writeMappings | samtools sort | samtools view -bS -o ${sample_id}.bam
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

workflow {
  QUANT(reads_ch, index_ch)
}
~~~
{: .language-groovy }


~~~
$ nextflow run process_publishDir_semantic.nf
~~~
{: .language-bash}

~~~
N E X T F L O W  ~  version 21.04.0
Launching `process_publishDir_semantic.nf` [golden_poisson] - revision: 421a604840

executor >  local (1)
[be/950786] process > QUANT (1) [100%] 1 of 1 ✔
~~~
{: .output}

We can now use the `tree` command to examine the results directory.

~~~
$ tree results
~~~
{: .language-bash }

~~~
results/
├── bams
│   └── ref1.bam
└── quant
    └── ref1_salmon_output
        ├── aux_info
        │   ├── ambig_info.tsv
        │   ├── expected_bias.gz
        │   ├── fld.gz
        │   ├── meta_info.json
        │   ├── observed_bias.gz
        │   └── observed_bias_3p.gz
        ├── cmd_info.json
        ├── libParams
        │   └── flenDist.txt
        ├── lib_format_counts.json
        ├── logs
        │   └── salmon_quant.log
        └── quant.sf

6 directories, 12 files
~~~
{: .output }




> ## Publishing results
>  Add a `publishDir` directive to the nextflow script `process_exercise_publishDir.nf` that copies the index directory to the results folder .
>  ~~~
> //process_exercise_pubilishDir.nf
> nextflow.enable.dsl=2
>
> process INDEX {
>    //add publishDir directive here
>    input:
>    path transcriptome
> 
>    output:
>    path "index"
>
>    script:
>    """
>    salmon index -t $transcriptome -i index
>    """
> }
>
> params.transcriptome = "data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
> transcriptome_ch = channel.fromPath(params.transcriptome, checkIfExists: true)
>
> workflow {
>   INDEX(transcriptome_ch)
> }
>  ~~~
>  {: .language-groovy }
> > ## Solution
> > ~~~
> > //process_exercise_publishDir_answer.nf
> > nextflow.enable.dsl=2
> >
> > process INDEX {
> >   publishDir "results", mode: "copy"
> >   
> >   input:
> >   path transcriptome
> >
> >   output:
> >   path "index"
> >
> >   script:
> >   """
> >   salmon index -t $transcriptome -i index
> >   """
> > }
> >
> > params.transcriptome = "data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
> > transcriptome_ch = channel.fromPath(params.transcriptome, checkIfExists: true)
> >
> > workflow{
> >   INDEX(transcriptome_ch)
> > }
> > ~~~
> > {: .language-groovy }
> >
> > ~~~
> > $ nextflow run process_exercise_publishDir.nf
> > ~~~~
> > {: .language-bash }
> > ~~~
> > N E X T F L O W  ~  version 21.04.0
> > Launching `process_exercise_publishDir.nf` [infallible_becquerel] - revision: 4d865241a8
> >
> > executor >  local (1)
> > [fe/79c042] process > INDEX (1) [100%] 1 of 1 ✔
> > ~~~
> > {: .output }
> {: .solution}
{: .challenge}

> ## Nextflow Patterns
> If you want to find out common structures of Nextflow processes, the [Nextflow Patterns page](http://nextflow-io.github.io/patterns/) collects some recurrent implementation patterns used in Nextflow applications.
{: .callout}



{: .output }
{% include links.md %}
