---
title: "Simple Rna-Seq pipeline"
teaching: 20
exercises: 40
questions:
- "How can I create a pipeline?"
objectives:
- "Create a simple RNA-Seq pipeline."

keypoints:
- ""
---

During this tutorial you will implement a proof of concept of a RNA-Seq pipeline which:

1. Indexes a transcriptome file.

1. Performs quality controls

1. Performs quantification.

1. Create a MultiqQC report.

## Define the pipeline parameters

The script `script1.nf` defines the pipeline input parameters.

~~~
params.reads = "$baseDir/data/yeast/reads/*_{1,2}.fq.gz"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"

println "reads: $params.reads"
~~~~
{: .language-groovy }

Run it by using the following command:

~~~
nextflow run script1.nf
~~~
{: language-bash}

Try to specify a different input parameter, for example:

~~~
nextflow run script1.nf --reads "sample1.fq"
~~~
{: .language-groovy }

> ## Add parameter
> Modify the `script1.nf` adding a fourth parameter named `outdir` and set it to a default path that will be used as the pipeline output directory.
> > ## Solution
> > ~~~
> > params.outdir = "results"
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

> # log.info
> Modify the `script1.nf` to print all the pipeline parameters by using a single `log.info` command and a multiline string statement.
> See an example [here](https://github.com/nextflow-io/rnaseq-nf/blob/3b5b49f/main.nf#L41-L48).
> > ## Solution
> > ~~~
> > log.info """\
> >         R N A S E Q - N F   P I P E L I N E    
> >         ===================================
> >         transcriptome: ${params.transcript}
> >         reads        : ${params.reads}
> >         outdir       : ${params.outdir}
> >         """
> >         .stripIndent()
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

### Recap

In this step you have learned:

* How to define parameters in your pipeline script

* How to pass parameters by using the command line

* The use of $var and ${var} variable placeholders

* How to use multiline strings

* How to use log.info to print information and save it in the log execution file

## Create a transcriptome index file

Nextflow allows the execution of any command or user script by using a process definition.

A process is defined by providing three main declarations: the process [inputs](https://www.nextflow.io/docs/latest/process.html#inputs), the process [outputs](https://www.nextflow.io/docs/latest/process.html#outputs) and finally the command [script](https://www.nextflow.io/docs/latest/process.html#script).

The second example adds the `index` process.

~~~
/*
 * pipeline input parameters
 */
params.reads = "$baseDir/data/yeast/reads/*_{1,2}.fq.gz"
params.transcriptome = "$baseDir/data/transcriptome/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

println """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


/*
 * define the `index` process that create a binary index
 * given the transcriptome file
 */
process index {

    input:
    path transcriptome from params.transcriptome

    output:
    path 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
~~~
{: .language-groovy }

It takes the transcriptome params file as input and creates the transcriptome index by using the `salmon` tool.

Note how the input declaration defines a transcriptome variable in the process context that it is used in the command script to reference that file in the Salmon command line.

Try to run it by using the command:

~~~
nextflow run script2.nf
~~~
{: .language-bash }

The execution will fail because Salmon is not installed in your environment.

Add the command line option `-with-conda` to launch the execution through a conda environment as shown below:

~~~
nextflow run script2.nf -with-conda
~~~
{: .language-bash }

This time it works because it uses the conda environment `nextflow/rnaseq-nf` defined in the `nextflow.config` file.

In add the option `-with-conda` add the following line in the `nextflow.config` file:

`docker.enabled = true`

> ## Enable docker
> Enable the Docker execution by default adding the above setting in the nextflow.config file.
> > ## Solution
> {: .solution}
{: .challenge}


> ## Print output of the index_ch
> Print the output of the index_ch channel by using the view.
> > ## Solution
> {: .solution}
{: .challenge}

> ## Examine work directory
> Use the command `tree work` to see how Nextflow organises the process work directory.
> > ## Solution
> > ~~~
> > tree work
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}

### Recap

In this step you have learned:

* How to define a process executing a custom command

* How process inputs are declared

* How process outputs are declared

* How to access the number of available CPUs

* How to print the content of a channel

## Collect read files by pairs

This step shows how to match read files into pairs, so they can be mapped by Salmon.

Edit the script `script3.nf` and add the following statement as the last line:
~~~
read_pairs_ch.view()
~~~
{: .language-groovy }

Save it and execute it with the following command:

~~~
nextflow run script3.nf
~~~
{: .language-bash }

It will print an output similar to the one shown below:

~~~
[ggal_gut, [/.../data/ggal/gut_1.fq, /.../data/yeast/reads/ref1_2.fq]]
~~~
{: .output }

The above example shows how the `read_pairs_ch` channel emits tuples composed by two elements, where the first is the read pair prefix and the second is a list representing the actual files.

Try it again specifying different read files by using a glob pattern:

~~~
nextflow run script3.nf --reads 'data/ggal/*_{1,2}.fq'
~~~
{: .language-bash }

File paths including one or more wildcards ie. `*`, `?`, etc. MUST be wrapped in single-quoted characters to avoid Bash expands the glob.

> ## `set` operator
> Use the set operator in place of = assignment to define the read_pairs_ch channel.
> > ## Solution
> > ~~~
> > Channel .fromFilePairs(params.reads)
> > .set{read_pairs_ch}
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

> ## checkIfExists
> Use the `checkIfExists` option for the `fromFilePairs` method to check if the specified path contains at least file pairs.
> > ## Solution
> > ~~~
> > Channel .fromFilePairs(params.reads, checkIfExists: true)
> > .set{read_pairs_ch}
> > {: .language-groovy }
> {: .solution}
{: .challenge}

### Recap

In this step you have learned:

* How to use `fromFilePairs` to handle read pair files

* How to use the `checkIfExists` option to check input file existence

* How to use the `set` operator to define a new channel variable

## Perform expression quantification

The script `script4.nf` adds the quantification process.

In this script note as the `index_ch` channel, declared as output in the index process, is now used as a channel in the input section.

Also note as the second input is declared as a tuple composed by two elements: the pair_id and the reads in order to match the structure of the items emitted by the read_pairs_ch channel.

Execute it by using the following command:

~~~
nextflow run script4.nf -resume
~~~
{: .source}
You will see the execution of the quantification process.

The `-resume` option cause the execution of any step that has been already processed to be skipped.

Try to execute it with more read files as shown below:

~~~
nextflow run script4.nf -resume --reads 'data/ggal/*_{1,2}.fq'
~~~~
{: .source}

You will notice that the quantification process is executed more than one time.

Nextflow parallelizes the execution of your pipeline simply by providing multiple input data to your script.

> # Add a tag directive
> Add a tag directive to the quantification process to provide a more readable execution log.
> > ## Solution
> >
> {: .solution}
{: .challenge}

> # Add a publishDir directive
Add a publishDir directive to the quantification process to store the process results into a directory of your choice.
> > ## Solution
> >
> {: .solution}
{: .challenge}

### Recap
In this step you have learned:

* How to connect two processes by using the channel declarations

* How to resume the script execution skipping already already computed steps

* How to use the tag directive to provide a more readable execution output

* How to use the publishDir to store a process results in a path of your choice

## Quality control

This step implements a quality control of your input reads. The inputs are the same read pairs which are provided to the quantification steps

You can run it by using the following command:

~~~
nextflow run script5.nf -resume
~~~
{: .language-bash}

The script will report the following error message:

~~~
Channel `read_pairs_ch` has been used twice as an input by process `fastqc` and process `quantification`
~~~
{: .output}

> Exercise
> Modify the creation of the read_pairs_ch channel by using a [into](https://www.nextflow.io/docs/latest/operator.html#into) operator in place of a set.
> > ## Solution
> >
> {: .solution}
{: .challenge}

see an example [here](https://github.com/nextflow-io/rnaseq-nf/blob/3b5b49f/main.nf#L58).

### Recap

In this step you have learned:

* How to use the `into` operator to create multiple copies of the same channel

## MultiQC report

This step collect the outputs from the quantification and fastqc steps to create a final report by using the MultiQC tool.

Execute the script with the following command:
~~~~
nextflow run script6.nf -resume --reads 'data/yeast/reads/*_{1,2}.fq.gz'
~~~~
{: .language-bash}

It creates the final report in the results folder in the current work directory.

In this script note the use of the `mix` and `collect` operators chained together to get all the outputs of the quantification and fastqc process as a single input.

### Recap

In this step you have learned:

* How to collect many outputs to a single input with the collect operator

* How to mix two channels in a single channel

* How to chain two or more operators togethers

## Handle completion event

This step shows how to execute an action when the pipeline completes the execution.

Note that Nextflow processes define the execution of asynchronous tasks i.e. they are not executed one after another as they are written in the pipeline script as it would happen in a common imperative programming language.

The script uses the `workflow.onComplete` event handler to print a confirmation message when the script completes.

Try to run it by using the following command:

~~~
nextflow run script7.nf -resume --reads 'data/yeast/reads/*_{1,2}.fq.gz'
~~~
{: .language-bash}
