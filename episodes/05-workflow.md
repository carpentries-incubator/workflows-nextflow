---
title: "Workflow"
teaching: 30
exercises: 15
questions:
- "How do I combined channels and processes to create a workflow?"
objectives:
- "Create a Nextflow workflow joining multiple processes."
keypoints:
- "A Nextflow workflow is define by invoking `processes` inside the `workflow` scope."
- "Process outputs can be accessed using the `out` attribute for the respective `process`."
---

## Workflows

We now know how to define tasks using processes. How do we then combine multiple process into a workflow?


## Workflow scope

We can combined process to create our pipeline inside our `workflow` scope.

### Process definition

As seen previously, a `process` is invoked as a function in the `workflow` scope, passing the expected input channels as parameters as it if were a function `<process_name>(<input_ch1>,<input_ch2>,...)`. To combined multiple process in a workflow invoke the processes in the order they would appear in a workflow.

For example:

~~~
nextflow.enable.dsl=2

process INDEX {
    input:
      path transcriptome
    output:
      path 'index'
    script:
      """
      salmon index -t $transcriptome -i index
      """
}

 process QUANT {
    input:
      each  path(index)
      tuple(val(pair_id), path(reads))
    output:
      path pair_id
    script:
      """
      salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
      """
}

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz',checkIfExists: true)
    reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz',checkIfExists: true)

    //index process takes 1 input channel as a parameter
    index_obj = INDEX(transcriptome_ch)

    //quant channel takes 2 input channels as parameters
    QUANT(index_out_ch,reads).view()
}
~~~
{: .language-groovy }

In this example, the `INDEX` process is invoked first and the `QUANT` process second.
The `INDEX` object, `index_obj`, is passed as the first parameter to the `QUANT` process. The `reads` channel is passed as the second parameter.


### Process composition

Processes having matching input-output declaration can be composed so that the output of the first process is passed as input to the following process.

Taking in consideration the previous process definition, it’s possible to write it as the following:

~~~
[..truncated..]

workflow {
  transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
  reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')

  // pass index process as a parameter to quant process
  quant(index(transcriptome_ch),reads ).view()
}
~~~
{: .language-groovy }

### Process outputs

A process output can also be accessed using the `out` attribute for the respective `process object`.

For example:

~~~
[..truncated..]

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
    reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
    index(transcriptome_ch)

    // process output  accessed using the `out` attribute of the index object index_out
    quant(index.out,reads)
    quant.out.view()
}
~~~
{: .language-groovy }

When a process defines two or more output channels, each of them can be accessed using the list element operator e.g. `out[0]`, `out[1]`, or using named outputs.

### Process named output

The process output definition allows the use of the `emit:` option to define a name identifier that can be used to reference the channel in the external scope.

For example:

~~~
nextflow.enable.dsl=2

process index {

  input:
  path transcriptome

  output:
  path 'index', emit: salmon_index

  script:
  """
  salmon index -t $transcriptome -i index
  """
}

process quant {
   input:
     each  path(index)
     tuple(val(pair_id), path(reads))
   output:
     path pair_id
   script:
     """
     salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
     """
}

workflow {
  transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
  reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
  index(transcriptome_ch)
  quant(index.out.salmon_index,reads).view()
}
~~~
{: .language-groovy }


### Workflow parameters

A workflow component can access any variable and parameter defined in the outer scope:

For example:

~~~
[..truncated..]

params.transcriptome = 'data/yeast/transcriptome/*.fa.gz'
params.reads = 'data/yeast/reads/*_{1,2}.fq.gz'

workflow {

  transcriptome_ch = channel.fromPath(params.transcriptome)
  reads = channel.fromFilePairs(params.reads)
  index(transcriptome_ch)
  quant(index.out.salmon_index,reads).view()
}
~~~
{: .language-groovy }

In this example the `params.transcriptome` and `params.reads` be accessed inside the `workflow` scope.


### Implicit workflow

A workflow definition which does not declare any name is assumed to be the main workflow, and it is implicitly executed. Therefore it’s the entry point of the workflow application.


> ## Fixme
> Add exercise to check  learner knows how to invoke a workflow and connect processes.
> >
> > ## Solution
> >
> {: .solution}
{: .challenge}

{% include links.md %}
