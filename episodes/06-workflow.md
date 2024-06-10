---
title: Workflow
teaching: 20
exercises: 20
---

::::::::::::::::::::::::::::::::::::::: objectives

- Create a Nextflow workflow joining multiple processes.
- Understand how to to connect processes via their inputs and outputs within a workflow.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I connect channels and processes to create a workflow?
- How do I invoke a process inside a workflow?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Workflow

Our previous episodes have shown us how to parameterise workflows using `params`, move data around a workflow using `channels` and define individual tasks using `processes`. In this episode we will cover how connect multiple processes to create a workflow.

## Workflow definition

We can connect processes to create our pipeline inside a `workflow` scope.
The  workflow scope starts with the keyword `workflow`, followed by an optional name and finally the workflow body delimited by curly brackets `{}`.

::::::::::::::::::::::::::::::::::::::::  callout

## Implicit workflow

In contrast to processes, the workflow definition in Nextflow does not require a name. In Nextflow, if you don't give a name to a workflow, it's considered the main/implicit starting point of your workflow program.

A named workflow is a `subworkflow` that can be invoked from other workflows, subworkflows are not covered in this lesson, more information can be found in the official documentation [here](https://www.nextflow.io/docs/latest/workflow.html).

::::::::::::::::::::::::::::::::::::::::::::::::::

### Invoking processes with a workflow

As seen previously, a `process` is invoked as a function in the `workflow` scope, passing the expected input channels as arguments as it if were.

```
 <process_name>(<input_ch1>,<input_ch2>,...)
```

To combined multiple processes invoke them in the order they would appear in a workflow. When invoking a process with multiple inputs, provide them in the same order in which they are declared in the `input` block of the process.

For example:

```groovy 
//workflow_01.nf



 process FASTQC {
    input:
      tuple(val(sample_id), path(reads))
    output:
      path "fastqc_${sample_id}_logs"
    script:
      """
      mkdir fastqc_${sample_id}_logs
      fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
      """
}

process MULTIQC {
    publishDir "results/mqc"
    input:
      path transcriptome
    output:
      path "*"
    script:
      """
      multiqc .
      """
}

workflow {
    read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz',checkIfExists: true)

    //index process takes 1 input channel as a argument
    //assign process output to Nextflow variable fastqc_obj
    fastqc_obj = FASTQC(read_pairs_ch)

    //quant channel takes 1 input channel as an argument
    //We use the collect operator to gather multiple channel items into a single item
    MULTIQC(fastqc_obj.collect()).view()
}
```

### Process outputs

In the previous example we assigned the process output to a Nextflow variable `fastqc_obj`.

A process output can also be accessed directly using the `out` attribute for the respective `process object`.

For example:

```groovy 
[..truncated..]

workflow {
  read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz',checkIfExists: true)

  FASTQC(read_pairs_ch)

  // process output  accessed using the `out` attribute of the process object
  MULTIQC(FASTQC.out.collect()).view()
  MULTIQC.out.view()

}
```

When a process defines two or more output channels, each of them can be accessed using the list element operator e.g. `out[0]`, `out[1]`, or using named outputs.

### Process named output

It can be useful to name the output of a process, especially if there are multiple outputs.

The process `output` definition allows the use of the `emit:` option to define a named identifier that can be used to reference the channel in the external scope.

For example in the script below we name the output from the `FASTQC` process as `fastqc_results` using the `emit:` option. We can then reference the output as
`FASTQC.out.fastqc_results` in the workflow scope.

```groovy 
//workflow_02.nf


 process FASTQC {
    input:
      tuple val(sample_id), path(reads)
    output:
      path "fastqc_${sample_id}_logs", emit: fastqc_results
    script:
      """
      mkdir fastqc_${sample_id}_logs
      fastqc -o fastqc_${sample_id}_logs ${reads}
      """
}

process MULTIQC {
    publishDir "results/mqc"
    input:
      path fastqc_results
    output:
      path "*"
    script:
      """
      multiqc .
      """
}

workflow {
    read_pairs_ch = channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz',checkIfExists: true)
    
    //FASTQC process takes 1 input channel as a argument
    FASTQC(read_pairs_ch)

    //MULTIQC channel takes 1 input channels as arguments
    MULTIQC(FASTQC.out.fastqc_results.collect()).view()
}
```

### Accessing script parameters

A workflow component can access any variable and parameter defined in the outer scope:

For example:

```groovy 
//workflow_03.nf
[..truncated..]

params.reads = 'data/yeast/reads/*_{1,2}.fq.gz'

workflow {

  reads_ch_ = channel.fromFilePairs(params.reads)
  FASTQC(reads_ch_)
  MULTIQC(FASTQC.out.fastqc_results.collect()).view()
}
```

In this example `params.reads`, defined outside the workflow scope, can be accessed inside the `workflow` scope.

:::::::::::::::::::::::::::::::::::::::  challenge

## Workflow

Connect the output of the process `FASTQC` to `PARSEZIP` in the Nextflow script `workflow_exercise.nf`.

**Note:** You will need to pass the `read_pairs_ch` as an argument to FASTQC and you will need to use the `collect` operator to gather the items in the FASTQC channel output to a single List item.

```groovy 
//workflow_exercise.nf

params.reads = 'data/yeast/reads/*_{1,2}.fq.gz'

process FASTQC {
 input:
 tuple val(sample_id), path(reads)

 output:
 path "fastqc_${sample_id}_logs/*.zip"

 script:
 """
 mkdir fastqc_${sample_id}_logs
 fastqc -o fastqc_${sample_id}_logs  ${reads}
 """
}

process PARSEZIP {
 publishDir "results/fqpass", mode:"copy"
 input:
 path fastqc_logs

 output:
 path 'pass_basic.txt'

 script:
 """
 for zip in *.zip; do zipgrep 'Basic Statistics' \$zip|grep 'summary.txt'; done > pass_basic.txt
 """
}
read_pairs_ch = channel.fromFilePairs(params.reads,checkIfExists: true)

workflow {
//connect process FASTQC and PARSEZIP
// remember to use the collect operator on the FASTQC output
}
```

:::::::::::::::  solution

## Solution

```groovy 
//workflow_exercise.nf



params.reads = 'data/yeast/reads/*_{1,2}.fq.gz'

process FASTQC {
  input:
  tuple val(sample_id), path(reads)

  output:
  path "fastqc_${sample_id}_logs/*.zip"

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs  ${reads}
  """
}

process PARSEZIP {
  publishDir "results/fqpass", mode:"copy"
  input:
  path fastqc_logs

  output:
  path 'pass_basic.txt'

  script:
  """
  for zip in *.zip; do zipgrep 'Basic Statistics' \$zip|grep 'summary.txt'; done > pass_basic.txt
  """
}

read_pairs_ch = channel.fromFilePairs(params.reads,checkIfExists: true)

workflow {
  PARSEZIP(FASTQC(read_pairs_ch).collect())
}
```

```bash 
$ nextflow run workflow_exercise.nf
```

```bash 
$ wc -l  results/fqpass/pass_basic.txt
```

```output 
18
```

The file `results/fqpass/pass_basic.txt` should have 18 lines.
If you only have two lines it might mean that you did not use `collect()` operator on the FASTC output channel.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::: keypoints

- A Nextflow workflow is defined by invoking `processes` inside the `workflow` scope.
- A process is invoked like a function inside the `workflow` scope passing any required input parameters as arguments. e.g. `FASTQC(reads_ch)`.
- Process outputs can be accessed using the `out` attribute for the respective `process` object or assigning the output to a Nextflow variable. 
- Multiple outputs from a single process can be accessed using the list syntax `[]` and it's index or by referencing the a named process output .

::::::::::::::::::::::::::::::::::::::::::::::::::


