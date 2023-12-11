---
title: Modules
teaching: 30
exercises: 15
---

::::::::::::::::::::::::::::::::::::::: objectives

- Add modules to a Nextflow script.
- Create a Nextflow modules.
- Understand how to use parameters in a module.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How can I reuse a Nextflow `process` in different workflows?
- How do I use parameters in a module?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Modules

In most programming languages there is the concept of creating code blocks/modules that can be reused.

Nextflow (DSL2) allows the definition of `module` scripts that can be included and shared across workflow pipelines.

A module file is nothing more than a Nextflow script containing one or more `process` definitions that can be imported from another Nextflow script.

A module can contain the definition of a `function`, `process` and `workflow` definitions.

For example:

```groovy 
process INDEX {
  input:
    path transcriptome
  output:
    path 'index'
  script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
```

The Nextflow process `INDEX` above could be saved in a file `modules/rnaseq-tasks.nf` as a Module script.

### Importing module components

A component defined in a module script can be imported into another Nextflow script using the `include` statement.

For example:

```groovy 
nextflow.enable.dsl=2

include { INDEX } from './modules/rnaseq-tasks'

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
    //
    INDEX(transcriptome_ch)
}
```

The above snippets includes a process with name `INDEX` defined in the module script `rnaseq-tasks.nf` in the main execution context, as such it can be invoked in the workflow scope.

Nextflow implicitly looks for the script file `./modules/rnaseq-tasks.nf` resolving the path against the including script location.

**Note:** Relative paths must begin with the `./` prefix.

::::::::::::::::::::::::::::::::::::::::  callout

## Remote

You can not include a script from a remote URL in the `from` statement.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Add module

Add the Nextflow module `FASTQC` from the Nextflow script `./modules/rnaseq-tasks.nf`
to the following workflow.

```groovy 
nextflow.enable.dsl=2


params.reads = "data/yeast/reads/ref1_{1,2}.fq.gz"
read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists:true )

workflow {
    FASTQC(read_pairs_ch)
}
```

> ## Solution
> 
> ```
> nextflow.enable.dsl=2
> include { FASTQC } from './modules/rnaseq-tasks'
> 
> params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
> read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists:true )
> 
> workflow {
>     FASTQC(read_pairs_ch)
> }
> ```

```
{: .language-groovy }
{: .solution}
```

::::::::::::::::::::::::::::::::::::::::::::::::::

### Multiple inclusions

A Nextflow script allows the inclusion of any number of modules.
When multiple components need to be included from the some module script,
the component names can be specified in the same inclusion using the curly brackets `{}`.

**Note** Component names are separated by a semi-colon `;` as shown below:

```groovy 
nextflow.enable.dsl=2

include { INDEX; QUANT } from './modules/rnaseq-tasks'

workflow {
    reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
    INDEX(transcriptome_ch)
    QUANT(index.out,reads)
}
```

### Module aliases

A process component, such as `INDEX`, can be invoked only once in the same workflow context.

However, when including a module component it's possible to specify a name alias using the keyword `as` in the `include` statement. This allows the inclusion and the invocation of the same component multiple times in your script using different names.

For example:

```groovy 
nextflow.enable.dsl=2

include { INDEX } from './modules/rnaseq-tasks'
include { INDEX as SALMON_INDEX } from './modules/rnaseq-tasks'

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
    INDEX(transcriptome_ch)
    SALMON_INDEX(transcriptome_ch)
}
```

In the above script the `INDEX` process is imported as `INDEX` and an alias `SALMON_INDEX`.

The same is possible when including multiple components from the same module script as shown below:

```groovy 
nextflow.enable.dsl=2

include { INDEX; INDEX as SALMON_INDEX } from './modules/rnaseq-tasks'

workflow {
  transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa.gz)'
  INDEX(transcriptome)
  SALMON_INDEX(transcriptome)
}
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Add multiple modules

Add the Nextflow modules `FASTQC` and `MULTIQC` from the Nextflow script `modules/rnaseq-tasks.nf`
to the following workflow.

```groovy 
nextflow.enable.dsl=2
params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists:true )

workflow {
   FASTQC(read_pairs_ch)
   MULTIQC(fastqc.out.collect())
}
```

:::::::::::::::  solution

## Solution

```groovy 
nextflow.enable.dsl=2
include { FASTQC; MULTIQC } from './modules/rnaseq-tasks'

params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists:true )

workflow {
   FASTQC(read_pairs_ch)
   MULTIQC(fastqc.out.collect())
}
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Module parameters

A module script can define one or more parameters using the same syntax of a Nextflow workflow script:

```groovy 
//functions.nf file
params.message = 'parameter from module script'

//The def keyword allows use to define a function that we can use in the code
def sayMessage() {
    println "$params.message"
}
```

Then, parameters are inherited from the including context. For example:

```groovy 
nextflow.enable.dsl=2

params.message = 'parameter from workflow script'

include {sayMessage} from './modules/functions'

workflow {
    sayMessage()
}
```

The above snippet prints:

```output 
parameter from workflow script
```

The module uses the parameters define before the include statement, therefore any further parameter set later is ignored.

**Tip:** Define all pipeline parameters at the beginning of the script before any include declaration.

The option `addParams` can be used to extend the module parameters without affecting the parameters set before the `include` statement.

For example:

```groovy 
nextflow.enable.dsl=2

params.message = 'parameter from workflow script'

include {sayMessage} from './modules/module.nf' addParams(message: 'using addParams')

workflow {
    sayMessage()
}
```

The above code snippet prints:

```output
using addParams
```



:::::::::::::::::::::::::::::::::::::::: keypoints

- A module file is a Nextflow script containing one or more `process` definitions that can be imported from another Nextflow script.
- To import a module into a workflow use the `include` keyword.
- A module script can define one or more parameters using the same syntax of a Nextflow workflow script.
- The module inherits the parameters define before the include statement, therefore any further parameter set later is ignored.

::::::::::::::::::::::::::::::::::::::::::::::::::


