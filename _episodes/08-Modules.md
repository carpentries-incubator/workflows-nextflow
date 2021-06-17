---
title: "Modules"
teaching: 30
exercises: 15
questions:
- "How do I modularise my pipeline?"
- "How can I reuse a workflow as part of another larger workflow?"
objectives:
- "Create Nextflow modules."
- "Create a sub-workflows."
keypoints:
- "A module file is a Nextflow script containing one or more process definitions that can be imported from another Nextflow script."
- "To import a module into a workflow use the `include` keyword."
---

## Modules

In most programming languages there is the concept of creating code blocks/modules that can be reused.

Nextflow (DSL2) allows the definition of `module` scripts that can be included and shared across workflow applications.

A module file is nothing more than a Nextflow script containing one or more process definitions that can be imported from another Nextflow script.  

A module can contain the definition of a `function`, `process` and `workflow` definitions as described in the above sections.

For example:

~~~
process index {
  input:
    path transcriptome
  output:
    path 'index'
  script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
~~~
{: .language-groovy }

The Nextflow process `index` above could be saved in a file `rnaseq.nf` as a Module script.

### Importing module components

A component defined in a module script can be imported into another Nextflow script using the `include` keyword.

For example:

~~~
nextflow.enable.dsl=2

include { index } from './modules/rnaseq.nf'

workflow {
    transcriptome = channel.fromPath('/some/data/*.txt')
    index(data)
}
~~~
{: .language-groovy }

The above snippets includes a process with name `index` defined in the module script `rnaseq.nf` in the main execution context, as such it can be invoked in the workflow scope.

Nextflow implicitly looks for the script file `./modules/rnaseq.nf` resolving the path against the including script location.

**Note:** Relative paths must begin with the `./` prefix.

> ## Add module
> Add the Nextflow module `fastqc` from the Nextflow script `modules/rnaseq.nf`
> to the following workflow.
> ~~~
> nextflow.enable.dsl=2
>
>
> params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
> reads_ch = channel.fromFilePairs( params.reads, checkIfExists:true )
>
> workflow {
>     fastqc(reads_ch)
> }
> ~~~~
> {: .language-groovy }
> > ## Solution
> > ~~~
> > nextflow.enable.dsl=2
> > include { fastqc } from './modules/rnaseq.nf'
> >
> > params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
> > reads_ch = channel.fromFilePairs( params.reads, checkIfExists:true )
> >
> > workflow {
> >     fastqc(reads_ch)
> > }
> ~~~~
> {: .language-groovy }
> {: .solution}
{: .challenge}

### Multiple inclusions

A Nextflow script allows the inclusion of any number of modules.
When multiple components need to be included from the some module script,
the component names can be specified in the same inclusion using the curly brackets `{}`.
Component names are separated by a semi-colon `;` as shown below:

~~~
nextflow.enable.dsl=2

include { index; quant } from './modules/rnaseq.nf'

workflow {
    reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
    transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa'
    index(transcriptome_ch)
    quant(index.out,reads)
}
~~~
{: .language-groovy }

### Module aliases

A process component, such as `index`, can be invoked only once in the same workflow context.

However, when including a module component it’s possible to specify a name alias using the keyword `as` in the `include` statement. This allows the inclusion and the invocation of the same component multiple times in your script using different names.

For example:

~~~
nextflow.enable.dsl=2

include { index } from './modules/rnaseq.nf'
include { index as salmon_index } from './modules/rnaseq.nf'

workflow {
    transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa'
    index(transcriptome)
    salmon_index(transcriptome)
}
~~~
{: .language-groovy }

In the above script the `index` process is imported as `index` and an alias `salmon_index`.

The same is possible when including multiple components from the same module script as shown below:

~~~
nextflow.enable.dsl=2

include { index; index as salmon_index } from './modules/rnaseq.nf'

workflow {
  transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa'
  index(transcriptome)
  salmon_index(transcriptome)
}
~~~
{: .language-groovy }


> ## Add multiple modules
> Add the Nextflow modules `fastqc` and `multiqc` from the Nextflow script `modules/rnaseq.nf`
> to the following workflow.
> ~~~
> nextflow.enable.dsl=2
> params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
> reads_ch = channel.fromFilePairs( params.reads, checkIfExists:true )
>
> workflow {
>    fq_out= fastqc(reads_ch)
>    multiqc(fq_out.collect())
> }
> ~~~~
> {: .language-groovy }
> > ## Solution
> > ~~~
> > nextflow.enable.dsl=2
> > include { fastqc; multiqc } from './modules/rnaseq.nf'
> >
> > params.reads = "$baseDir/data/yeast/reads/ref1_{1,2}.fq.gz"
> > reads_ch = channel.fromFilePairs( params.reads, checkIfExists:true )
> >
> > workflow {
> >    fq_out = fastqc(reads_ch)
> >    multiqc(fq_out.collect())
> > }
> ~~~~
> {: .language-groovy }
> {: .solution}
{: .challenge}

### Module parameters

A module script can define one or more parameters using the same syntax of a Nextflow workflow script:

~~~
params.message = 'parameter from module script'

def sayMessage() {
    println "$params.message $params.bar"
}
~~~
{: .language-groovy }

Then, parameters are inherited from the including context. For example:

~~~
nextflow.enable.dsl=2


params.message = 'parameter from workflow script'

include {sayMessage} from './modules/module.nf'

workflow {
    sayMessage()
}
~~~
{: .language-groovy }

The above snippet prints:

~~~
parameter from workflow script
~~~
{: .output }

The module inherits the parameters define before the include statement, therefore any further parameter set later is ignored.

**Tip** Define all pipeline parameters at the beginning of the script before any include declaration.

The option `addParams` can be used to extend the module parameters without affecting the external scope. For example:

~~~
nextflow.enable.dsl=2

params.message = 'parameter from workflow script'

include {sayMessage} from './modules/module.nf' addParams(message: 'using addParams')

workflow {
    sayMessage()
}
~~~
{: .source}

The above snippet prints:
~~~
using addParams
~~~
{: .output}

## Sub-workflows

The DSL2 syntax also allows for the definition of sub-workflow libraries. The only requirement is to provide a `workflow` name that will be used to reference and declare the corresponding inputs and outputs using the new `take` and `emit` keywords.

For example:

~~~
nextflow.enable.dsl=2

workflow rnaseq {
  take:
    transcriptome
    read_pairs_ch

  main:
    index(transcriptome)
    fastqc(read_pairs_ch)
    quant(INDEX.out, read_pairs_ch)

  emit:
     quant.out.mix(fastqc.out).collect()
}
~~~
{: .language-groovy }

Now named sub-workflows can be used in the same way as processes, allowing you to easily include and reuse multi-step workflows as part of larger workflows.

~~~
nextflow.enable.dsl=2

process index {
    index:
      path transcriptome
    output:
      path 'index'
    script:
      """
      salmon index -t $transcriptome -i index
      """
}

 process quant {
    input:
      tuple pair_id, path(reads)
      path index
    output:
      path pair_id
    script:
      """
      salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
      """
}
`
workflow {
    transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa)
    reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
    index(transcriptome_ch)
    quant(index.out,reads)
}
~~~
{: .language-groovy }



### Workflow definition

The `workflow` keyword allows the definition of sub-workflow components that enclose the invocation of one or more processes and operators:

~~~
workflow salmon_quant {
  reads = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
  transcriptome_ch = channel.fromPath('/data/yeast/transcriptome/*.fa')
  quant(reads, index(transcriptome_ch))
}
~~~
{: .language-groovy }

For example, the above snippet defines a workflow component, named `salmon_quant`, that can be invoked from another workflow component definition as any other function or process i.e. `my_rnaseq_pipeline`


~~~
workflow my_rnaseq_pipeline {
  salmon_quant()
}
~~~
{: .language-groovy }


### Implicit workflow

A workflow definition which does not declare any name is assumed to be the main workflow, and it is implicitly executed. Therefore it’s the entry point of the workflow application.

> ## Note
> Implicit workflow definition is ignored when a script is included as module. This allows the writing a workflow script that can be used either as a library module and as application script.
{: .callout}


An alternative workflow entry can be specified using the `-entry` command line option.

### Workflow composition

Workflows defined in your script or imported by a module inclusion can be invoked and composed as any other process in your application.

~~~
workflow flow1 {
    take: data
    main:
        foo(data)
        bar(foo.out)
    emit:
        bar.out
}

workflow flow2 {
    take: data
    main:
        foo(data)
        baz(foo.out)
    emit:
        baz.out
}

workflow {
    take: data
    main:
      flow1(data)
      flow2(flow1.out)
}
~~~
{: .language-groovy }

> ## Nested workflow execution
> Nested workflow execution determines an implicit scope. Therefore the same process can be invoked in two different workflow scopes, like for example foo in the above snippet that is used either in `flow1` and `flow2`. The workflow execution path along with the process names defines the process fully qualified name that is used to distinguish the two different process invocations i.e. `flow1:foo` and `flow2:foo` in the above example.
{: .callout}  
