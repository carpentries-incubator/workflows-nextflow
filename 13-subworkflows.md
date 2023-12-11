---
title: Sub-workflows
teaching: 20
exercises: 0
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand how to create a sub-workflow.
- Understand how to run part of a workflow.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I reuse a workflow as part of a larger workflow?
- How do I run only a part of a workflow?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Sub-workflows

We have seen previously the Nextflow DSL2 syntax allows for the definition of reusable processes (modules). Nextflow DSL2 also allow the definition reusable  sub-workflow libraries.

## Workflow definition

The `workflow` keyword allows the definition of workflow components that enclose the invocation of one or more `processes` and `operators`.

For example,:

```groovy 
nextflow.enable.dsl=2

include {QUANT;INDEX} from './modules/module.nf'

workflow RNASEQ_QUANT_PIPE {
  read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
  transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa.gz')
  QUANT(INDEX(transcriptome_ch),read_pairs_ch)
}
```

The above snippet defines a workflow component, named `RNASEQ_QUANT_PIPE`, that can be invoked from another workflow component definition in the same way as any other function or `process` i.e. `RNASEQ_QUANT()`.

```groovy 
nextflow.enable.dsl=2

include {QUANT;INDEX} from './modules/module.nf'

workflow RNASEQ_QUANT_PIPE {
  read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
  transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa.gz')
  QUANT(INDEX(transcriptome_ch),read_pairs_ch)
}

// Implicit workflow
workflow  {
  /*
  * Call sub-workflow using <WORKFLOWNAME>() syntax
  */
  RNASEQ_QUANT_PIPE()
}
```

::::::::::::::::::::::::::::::::::::::::  callout

## Implicit workflow

A workflow definition which does not declare any name is assumed to be the main workflow, and it is implicitly executed. Therefore it's the entry point of the workflow application.


::::::::::::::::::::::::::::::::::::::::::::::::::

### Workflow parameters

A workflow component can access any variable and parameter defined in the outer scope.

For Example:

```groovy 
nextflow.enable.dsl=2

include {QUANT;INDEX} from './modules/module.nf'

params.transcriptome = '/some/data/file'
read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')

workflow RNASEQ_QUANT_PIPE {

  transcriptome_ch = channel.fromPath(params.transcriptome)
  QUANT(INDEX(transcriptome_ch),read_pairs_ch)
}
```

### Workflow inputs

A workflow component can declare one or more input channels using the `take` keyword.

For example:

```groovy 
nextflow.enable.dsl=2

include {QUANT;INDEX} from './modules/module.nf'

params.transcriptome = '/some/data/file'
read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')

workflow RNASEQ_QUANT_PIPE {
    take:
      transcriptome_ch
      read_pairs_ch
    main:
      transcriptome_ch = channel.fromPath(params.transcriptome)
      INDEX(transcriptome_ch)
      QUANT(INDEX.out,read_pairs_ch)
}
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Warning

When the `take` keyword is used, the beginning of the workflow body needs to be identified with the `main` keyword.
Then, the input can be specified as an argument in the workflow invocation statement:


::::::::::::::::::::::::::::::::::::::::::::::::::

These input channels can then be passed to the workflow as parameters inside the `()`. Multiple parameters are separated by a comma `,` and must be specified in the order they appear under `take`:

```groovy 
nextflow.enable.dsl=2

include {QUANT;INDEX} from './modules/module.nf'

params.transcriptome = '/some/data/file'
read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')

workflow RNASEQ_QUANT_PIPE {
    take:
      transcriptome_ch
      read_pairs_ch
    main:
      transcriptome_ch = channel.fromPath(params.transcriptome)
      INDEX(transcriptome_ch)
      QUANT(INDEX.out,read_pairs_ch)
}

workflow {
    RNASEQ_QUANT_PIPE(transcriptome_ch,read_pairs_ch )
}
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Note

Workflow inputs are by definition channel data structures. If a basic data type is provided instead, ie. number, string, list, etc. it's implicitly converted to a channel value (ie. non-consumable).


::::::::::::::::::::::::::::::::::::::::::::::::::

### Workflow outputs

A workflow component can declare one or more output channels using the `emit` keyword.

For example:

```groovy 
nextflow.enable.dsl=2

include {QUANT;INDEX} from './modules/module.nf'

params.transcriptome = '/some/data/file'
read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')

workflow RNASEQ_QUANT_PIPE {
    take:
     transcriptome_ch
     read_pairs_ch
    emit:
      QUANT.out
    main:
      transcriptome_ch = channel.fromPath(params.transcriptome)
      INDEX(transcriptome_ch)
      QUANT(INDEX.out,read_pairs_ch)
}
```

The above script declares one output, `QUANT.out`.

The result of the `RNASEQ_QUANT_PIPE` execution can be accessed using the `out` property ie. `RNASEQ_QUANT_PIPE.out`.

When there are multiple output channels declared, use the array bracket notation to access each output component as described for the Process outputs definition.

```groovy 
RNASEQ_QUANT_PIPE.out[0]
RNASEQ_QUANT_PIPE.out[1]
```

Alternatively, the output channel can be accessed using a name which it's assigned to in the emit declaration:

For example:

```groovy 
nextflow.enable.dsl=2

workflow RNASEQ_QUANT_PIPE {
   main:
     INDEX(transcriptome_ch)
     QUANT(INDEX.out,read_pairs_ch)
   emit:
     read_quant = QUANT.out
}
```

The output `QUANT.out` is assigned the name `read_quant`
The the result of the above snippet can accessed using:

```groovy 
RNASEQ_QUANT_PIPE.out.read_quant`.
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Note

Implicit workflow definition is ignored when a script is included as module. This allows the writing a workflow script that can be used either as a library module and as application script.


::::::::::::::::::::::::::::::::::::::::::::::::::

### Workflow composition

As with `modules` workflows components can be defined within your script or imported by a `include` statment. After which thet can then be invoked and composed as any other `workflow component` or process in your script.

```groovy
nextflow.enable.dsl=2

// file modules/qc.nf
include {FASTQC} from './modules.nf'

workflow READ_QC_PIPE {
    take:
      read_pairs_ch
      quant_out_ch
    main:
        FASTQC(read_pairs_ch)
    emit:
        FASTQC.out
}
```

```source
nextflow.enable.dsl=2

include { READ_QC_PIPE } from './modules/qc.nf'

workflow RNASEQ_QUANT_PIPE {
    take:
      transcriptome_ch
      read_pairs_ch
    main:
        INDEX(transcriptome)
        QUANT(INDEX.out)
    emit:
        QUANT.out
}

params.transcriptome = '/some/data/file'
read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
transcriptome_ch = channel.fromPath(params.transcriptome)

workflow {
    take:
      transcriptome_ch
      read_pairs_ch
    main:
      RNASEQ_QUANT(transcriptome_ch,read_pairs_ch)
      READ_QC(read_pairs_ch,RNASEQ_QUANT.out)
      MULTIQC(RNASEQ_QUANT.out.mix(READ_QC).collect())
}
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Nested workflow execution

Nested workflow execution determines an implicit scope. Therefore the same process can be invoked in two different workflow scopes, like for example in the above snippet `INDEX` could be used either in `RNASEQ_QUANT` and `RNASEQ_QC`. The workflow execution path along with the process names defines the process fully qualified name that is used to distinguish the two different process invocations i.e. `RNASEQ_QUANT:INDEX` and `RNASEQ_QC:INDEX` in the above example.


::::::::::::::::::::::::::::::::::::::::::::::::::

### Specific workflow entry points

By default, the unnamed workflow is assumed to be the main entry point for the script. Using named workflows, the entry point can be customised by using the `entry` option of the `run` command. This allows users to run a specific sub-workflow or a section of their entire workflow script.

For example:

```bash
$ nextflow run main.nf -entry RNASEQ_QUANT_PIPE
```

The above command would run the `RNASEQ_QUANT_PIPE` sub-workflow.



:::::::::::::::::::::::::::::::::::::::: keypoints

- Nextflow allows for definition of reusable sub-workflow libraries.
- Sub-workflow allows the definition of workflow processes that can be included from any other script and invoked as a custom function within the new workflow scope. This enables reuse of workflow components
- The `entry` option of the nextflow `run` command specifies the workflow name to be executed

::::::::::::::::::::::::::::::::::::::::::::::::::


