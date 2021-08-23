---
title: "Nextflow configuration"
teaching: 30
exercises: 15
questions:
- "How can I configure how Nextflow runs?"
- "How can I write a Nextflow configuration file?"
- "How can I control `process` settings using the Nextflow configuration file?"
- "How can I assign different resources to different process?"
- "How can I assign different profiles to different computational systems "
objectives:
- "Understand how Nextflow is configured."
- "Create a Nextflow configuration file."
- "Understand how to use the `process` scope to define process settings."
- "Understand how to specify custom resources with process selectors."
- "Understand how to group configurations using the `profile` scope. "
keypoints:
- "Nextflow configuration can be managed using a Nextflow configuration file."
- "Nextflow configuration are simple text files containing a set of properties."
- "You can define process setting, such as cpus and memory, within the `process` scope."
- "You can assign different resources to different process using the process selectors `withNames` or `withLabel`. "
- "You can define a profile for different configurations using the `profiles` scope. These profiles can be selected when launching a pipeline execution by using the `-profile` command line option"
---

## Nextflow configuration

A key Nextflow feature is the ability to decouple the workflow implementation from the configuration settings required by the underlying execution platform.

This enable portable deployment of the workflow on different computational platforms such as an institutional HPC or AWS without the need to modify the workflow code.

Nextflow uses configuration files to achieve this separation.

## Configuration file

We have seen in previous episodes how to configure how a workflow runs using parameters specified on the command line (`--variable_name`). You can also specify workflow parameters and settings using a Nextflow configuration file.

When a pipeline script is launched Nextflow looks for a file named `nextflow.config` in the current directory and in the script base directory (if it is not the same as the current directory). Finally it checks for the file `$HOME/.nextflow/config` (Note this file is not called nextflow.config, just config).

> ## $HOME/.nextflow/config
> `$HOME/.nextflow/config` setting are alway imported into your runs. This file is normally where you save the common setting for all your runs.
{: .callout }

The default config file search mechanism can be extended proving an extra configuration file by using the command line option `-c <config file>`.

When more than one on the above files exist they are merged, so that the settings in the first override the same ones that may appear in the second one, and so on. The configuration sources are ranked to decide which settings to apply to the workflow.

Possible configuration sources are reported below, listed in order of priority:

1. Parameters specified on the command line (`--something value`)
1. Config file specified using the `-c` my_config option
1. The config file named `nextflow.config` in the current directory
1. The config file named `nextflow.config` in the workflow project directory
1. The config file `$HOME/.nextflow/config`
1. Values defined within the pipeline script itself.

## Config syntax

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax:

`name = value`

String values need to be wrapped in quotation characters while numbers and boolean values (`true`, `false`) do not.

**Note:** that values are typed, meaning for example that, `1` is different from `'1'`, since the first is interpreted as the number one, while the latter is interpreted as a string value.

### Config variables

Configuration properties can be used as variables in the configuration file itself, by using the usual `$propertyName` or `${expression}` syntax.

**Note:** Importantly, these variables are not available in the Nextflow script.

For example:

~~~
//nextflow.config
kmer = 27
kmer_message = "kmer size is  ${kmer}"
~~~
{: .language-groovy}

You can use the `nextflow config` command to  to print the resolved configuration of the  `nextflow.config`  and is especially useful for understanding the resolved profiles and parameters that Nextflow will use run a pipeline.

~~~
$ nextflow config
~~~
{: .language-bash}

Would output:

~~~
kmer = 27
kmer_message = 'kmer size is  27'
~~~
{: .output}

In the configuration file it’s possible to access any variable defined in the host environment such as `$PATH`, `$HOME`, `$PWD`, etc.

~~~
//nextflow.config
my_home_dir = "$HOME"
~~~
{: .source}

### Config comments

Configuration files use the same conventions for comments used in the Nextflow script:

~~~
// comment a single config file

/*
   a comment spanning
   multiple lines
 */
~~~
{: .language-groovy }

### Config scopes

Configuration settings can be grouped together in a scope. We can create scopes by dot prefixing `.` the property names with a scope identifier or grouping the properties in the same scope using the scope identifier and curly brackets notation `{}`.

Both way of grouping a properties via a scope is shown in the following example:

~~~
//nextflow.config
aligner.name  = "salmon"
aligner.kmer  = 27

index {
    kmer = 27
    outdir = 'results/index'
}
~~~
{: .source }

As mentioned previously, these properties/variable in user defined scopes are not accessible in the main Nextflow script.

However, there are a number of special config scopes that are accessible and can enable you to alter various aspect of your workflow's run such as `params`, `process`, `env`, `executor` and `conda`.

### Parameter scope

The scope `params` allows the definition of workflow parameters that overrides the values defined in the main workflow script.

This is useful to consolidate one or more execution parameters in a separate file.

~~~
//nextflow.config
params.kmer= 27
params.transcriptome = "data/yeast/transcriptome/*.fa.gz"
~~~
{: .source }

~~~
// workflow script
params.kmer = 31

println "$params.kmer"
~~~~
{: .language-groovy }

A variable defined in the `nextflow.config` file have priority over those in the main pipeline script e.g. `main.nf`,  this would output:

~~~
$ nextflow run main.nf
~~~
{: .language-bash }

~~~
27
~~~
{: .output }

**Note:** A param defined on the command line has the highest priority.

For example:

~~~
$ nextflow run main.nf --kmer 21
~~~
{: .language-bash }

~~~
21
~~~
{: .output }

> ## Configuration parameters
> Save the code below in the file `nextflow.config`
> ~~~
> // config file
> params.genome = 'GRCh38'
> params.aligner = 'salmon'
> ~~~
> {: .source }
> and then save the code below in the Nextflow script `params.nf`.
> ~~~
> // workflow script
> params.genome = 'GRCh37'
> params.aligner = 'kallisto'
> println "$params.genome $params.aligner"
> ~~~
> {: .language-groovy }
> Then run:
> ~~~
> nextflow run params.nf
> ~~~
> {: .language-bash}
> Execute is again specifying the `genome` parameter on the command line:
>
> ~~~
> nextflow run params.nf --genome hg38
> ~~~
> {: .language-bash}
>
> Compare the result of the two executions.
> > ## Solution
> > The first command will print
> > ~~~
> > GRCh38 salmon
> > ~~~
> > {: .output}
> > The second command will print
> > ~~~
> > hg38 salmon
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}


### Environment scope

The `env` scope allows the definition one or more environmental variable that will be accessible in the environment where the workflow tasks will be executed e.g. .

Simply prefix your variable names with the `env` scope or surround them by curly brackets, as shown below:

~~~~
// configuration file
env.genome = 'hg38'
env {
   aligner = 'salmon'  
}
~~~~
{: .source }

You can now access the environment variables using the variable name without the `env` prefix.

~~~
// workflow script

println "aligner is ${aligner}"

process ENV {
  script:
  """
  env |grep 'genome'
  """
}

workflow {
  ENV()
}
~~~
{: .language-groovy }

~~~
aligner is salmon
hg38
~~~
{: .output }

### Process scope

In the process episode we saw that the `process` directives allow the specification of settings for the task execution such as `cpus`, `memory`, `conda` and other resources in the pipeline script.

This is useful when prototyping a small workflow script.

However it’s always a good practice to decouple the workflow execution logic from the process configuration settings, i.e. it’s strongly suggested to define the process settings in the workflow configuration `nextflow.config` file instead of the workflow script.

The `process` configuration scope allows the setting of any process directives in the Nextflow configuration file.

For example:

~~~
//nextflow.config
process {
    cpus = 2
    memory = 8.GB
}
~~~
{: .source }

> ## process directives require `=` in configuration file.
> The syntax for setting process directives in the configuration file requires `=` ie. assignment  operator.
{: .callout}




The above config snippet defines the `cpus`, `memory` for **all** processes in your workflow script.


The `process` selector can be used to apply the configuration to a specific process or group of processes (discussed later).

Memory and time duration unit can be specified either using a string based notation in which the digit(s) and the unit can be separated by a blank or
by using the numeric notation in which the digit(s) and the unit are separated by a dot character and it’s not enclosed by quote characters.

  | String syntax   | Numeric syntax | Value                 |
  |-----------------|----------------|-----------------------|
  | '10 KB'         | 10.KB          | 10240 bytes           |
  | '500 MB'        | 500.MB         | 524288000 bytes       |
  | '1 min'         | 1.min          | 60 seconds            |
  | '1 hour 25 sec' | -              | 1 hour and 25 seconds |


### Dynamic expressions

We can define a config setting using a dynamic expression using a closure.

For example we can specify the amount of `memory` required as a multiple of them
number of `cpus`. To access the number of cpus for a task use `task.cpus`.

For example:

~~~
//nextflow.config
process {
    cpus = 2
    memory = { 2.GB * task.cpus }
}
~~~
{: .source }


> ## Configure process scope
> 1. Create a Netxflow config file `wc-params.config`
> 1. Add a  process scope specifying the process run time as `time = '5s'``
> 1. Then run:
> ~~~
> $ nextflow run wc-params.nf --sleep 10 -c wc-params.config
> ~~~
> {: .language-bash}
> What output do you get?
> > ## Solution
> > ~~~
> > //wc-params.config
> > process {
> >
> > time = '5s'
> >
> > }
> > ~~~
> > {: .source}
> >
> > You will get a runtime error:
> > ~~~
> > Process exceeded running time limit (5s)
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}



### Process executor

A key Nextflow feature is the ability to decouple the workflow implementation from the actual execution platform. This allows the deployment of a workflow on any executing platform support by the Nextflow.

![nf-executors](https://seqera.io/training/img/nf-executors.png)

The `executor` defines the underlying system where processes are executed. By default a process uses the `executor` defined globally in the `nextflow.config` file.

The `process.executor` directive allows you to configure what executor has to be used by the process, overriding the default configuration.

For example:

|Name| Executor|
|----|---------|
|`local`|The process is executed in the computer where Nextflow is launched.|
|`sge`|The process is executed using the Sun Grid Engine / Open Grid Engine.|
|`uge`|The process is executed using the Univa Grid Engine job scheduler.|
|`lsf`|The process is executed using the Platform LSF job scheduler.|
|`slurm`|The process is executed using the SLURM job scheduler.|
|`pbs`|The process is executed using the PBS/Torque job scheduler.|

> ## Executors
> See full list of executors here [here](nextflow.io/docs/latest/config.html#scope-executor)
{: .callout }


To run your pipeline with a batch scheduler modify the `nextflow.config` file specifying the target executor and the required computing resources if needed. For example:

~~~
//nextflow.config
process.executor = 'sge'
~~~
{: .source }

#### Executor scope

The executor has it's own config scope which allows you to set optional executor settings.

For example, in the config below we specify the executor as Sun Grid Engine, `sge` and the number of tasks the executor will handle in a parallel manner  `queueSize` to 10.

~~~
//nextflow.config
executor {
    name = 'sge'
    queueSize = 10
}
~~~

More info [here](https://www.nextflow.io/docs/latest/config.html#scope-executor).

### Process selectors

In real world application different tasks need different amount of computing resources.

#### Process name

It is possible to define the resources for a specific task using  `withName:` followed by the process name and the directives within curly braces.

For example for the process `INDEX` and `FASTQC` we can specify different `cpus` and `memory` resources.:

~~~
//nextflow.config
process {
    withName: INDEX {
        cpus = 4
        memory = 8.GB
    }
    withName: FASTQC {
        cpus = 2
        memory = 4.GB
    }
}
~~~
{: .source }


#### Process labels

When a workflow is composed of many processes it can be overkill listing all process names in the configuration file to specifies the resources for each of them.

A better strategy consist to annotate the processes with a `label` directive. Then specify the resources in the configuration file using for all processes having the same label.

The `withLabel` selectors allow the configuration of all processes annotated with a label directive as shown below:

~~~
nextflow.enable.dsl=2

process P1 {
label "bigmem"
script:
  """
  echo P1: Using $task.cpus cpus and $task.memory memory.
  """
}

process P2 {
label "bigmem"

script:
  """
  echo P2: Using $task.cpus cpus and $task.memory memory.
  """
}

workflow {
 P1()
 P2()
}
~~~
{: .language-groovy }

~~~
//nextflow.config
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
    }
}
~~~
{: .source }

The above config example assigns 16 cpus, and 64 Gb of memory to all processes annotated with the label `big_mem`.


### Selectors priority

When mixing generic process configuration and selectors the following priority rules are applied (from lower to higher):

1. Process generic `process` configuration.
1. Process specific directive defined in the workflow script.
1. `withLabel` selector definition.
1. `withName` selector definition.

> ## Process selectors
> Create a Nextflow config, `nextflow.config`, file specifying different
> 1. `cpus`
> 1. `memory`
>
> resources for the two process `P1` (cpus 1 and memory 2.GB) and `P2` (cpus2 and memory 1.GB) for the Nextflow script below.
>
> ~~~
> nextflow.enable.dsl=2
>
> process P1 {
>
> script:
>   """
>   echo P1: Using $task.cpus cpus and $task.memory memory.
>   """
> }
>
> process P2 {
>
> script:
>   """
>   echo P2: Using $task.cpus cpus and $task.memory memory.
>   """
> }
>
> workflow {
>  P1()
>  P2()
> }
> ~~~
> > ## Solution
> > ~~~
> > //config file
> > process {
> >  echo = true
> >    withName: P1 {
> >    cpus = 1
> >    memory = 2.GB
> >   }
> >    withName: P2 {
> >    cpus = 2
> >    memory = 1.GB
> >   }
> > }
> > ~~~
> {: .solution}
 {: .challenge}

## Defining software requirements


A key feature of Nextflow is the ability to use different technologies such as the Conda package management system or container engines such as Docker, Singularity, Podman, Charliecloud, Shifter to manage the software requirements. This also facilitate portable and reproducible workflows.

We can specify the technology used for the entire workflow or a single or group of process using the technology specific scopes.

## Config Conda execution

To use  a Conda environment you can either specify the path of an existing Conda environment directory or the path of Conda environment `YAML` file process scope of the config file:

~~~
process.conda = "/home/user/miniconda3/envs/my_conda_env"
//or
process {
  conda = "/home/user/miniconda3/envs/my_conda_env"
}
~~~
{: .source }

Or specify the path of Conda environment `YAML` file.

For example;

~~~
process.conda = "environment.yml"
~~~
{: .source }

### conda scope

There is an optional `conda` scope allows you to control the creation of a Conda environment by the Conda package manager.

For example, `cacheDir` specifies the path  where the Conda environments are stored. By default this is in `conda` folder of the work directory.

**Note:** When using a compute cluster make sure to provide a shared file system path accessible from all computing nodes.


## Config Docker execution

To use docker we specify the container image to be used in the `process` scope
and set `docker.enabled = true` in the `docker` scope.

For example:

~~~
process.container = 'nextflow/rnaseq-nf'
docker.enabled = true
~~~
{: .source }

### Config Singularity execution

Similar to docker to use singularity container we provide the container image file path in `process.container` and set `singularity.enabled = true`:

~~~
process.container = '/some/singularity/image.sif'
singularity.enabled = true
~~~
{: .language-groovy }

**Note::** The container image file must be an absolute path i.e. it must start with a `/`.

> ## Container protocols
> The following protocols are supported:
>* `library://`` download the container image from the Singularity Library service.
>* `shub://`` download the container image from the Singularity Hub.
>* `docker://`` download the container image from the Docker Hub and convert it to the Singularity format.
>* `docker-daemon://` pull the container image from a local Docker installation and convert it to a Singularity image file.
{: .callout}

> Docker to Singularity image
Specifying a plain Docker container image name, Nextflow implicitly download and converts it to a Singularity image when the Singularity execution is enabled.
For example:
> ~~~
> process.container = 'nextflow/rnaseq-nf'
> singularity.enabled = true
> ~~~
> {: .source}
The above configuation instructs Nextflow to use Singularity engine to run your script processes.
>The container is pulled from the Dockr registry and cached in the current directory to be used for further runs.
>Alternatively if you have a Singularity image file, its location absolute path can be specified as the container name either using the `-with-singularity` option on the command line or the `process.container` setting in the config file.
{: .callout}



## Configuration profiles

One of the most powerful features of the configuration file is to define multiple different configuration or `profiles` , for different computational environments e.g. local computer vs. HPC.

A profile is a set of configuration attributes that can be activated/chosen when launching a pipeline execution by using the `-profile` command line option.

Configuration profiles are defined by using the special scope `profiles` which group the attributes that belong to the same profile using a common prefix.


For example:

~~~
profiles {

    local {
        params.genome = '/local/path/ref.fasta'
        process.executor = 'local'
    }

    cluster {
        params.genome = '/data/stared/ref.fasta'
        process.executor = 'sge'
        process.queue = 'long'
        process.memory = '10GB'
        process.conda = '/some/path/env.yml'
    }

    cloud {
        params.genome = '/data/stared/ref.fasta'
        process.executor = 'awsbatch'
        process.container = 'cbcrg/imagex'
        docker.enabled = true
    }

}
~~~
{: .source}

This configuration defines three different profiles: `standard`, `cluster` and `cloud` that set different process configuration strategies depending on the target runtime platform. By convention the standard profile is implicitly used when no other profile is specified by the user.

To enable a specific profile use `-profile` option followed by the profile name:

~~~
nextflow run <your script> -profile cluster
~~~
{: .language-bash}

{% include links.md %}
