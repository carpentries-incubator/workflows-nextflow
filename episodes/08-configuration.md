---
title: "Nextflow configuration"
teaching: 30
exercises: 15
questions:
- "What is the difference between the workflow implementation and the workflow configuration?"
- "How do I configure a Nextflow workflow?"
- "How do I assign different resources to different processes?"
- "How do I separate and provide configuration for different computational systems?"
- "How do I change configuration settings from the default settings provided by the workflow?"
objectives:
- "Understand the difference between workflow implementation and configuration."
- "Understand the difference between configuring Nextflow and a Nextflow script."
- "Create a Nextflow configuration file."
- "Understand what a configuration scope is."
- "Be able to assign resources to a process."
- "Be able to refine configuration settings using process selectors."
- "Be able to group configurations into profiles for use with different computer infrastructures."
- "Be able to override existing settings."
- "Be able to inspect configuration settings before running a workflow."
keypoints:
- "Nextflow configuration can be managed using a Nextflow configuration file."
- "Nextflow configuration files are plain text files containing a set of properties."
- "You can define process setting, such as cpus and memory, within the `process` scope."
- "You can assign different resources to different process using the process selectors `withNames` or `withLabel`. "
- "You can define a profile for different configurations using the `profiles` scope. These profiles can be selected when launching a pipeline execution by using the `-profile` command line option"
- "The workflow configuration settings can be inspected using `nextflow config <script> [options]`."
---

## Nextflow configuration

A key Nextflow feature is the ability to decouple the workflow implementation, which describes
the flow of data and operations to perform on that data, from the configuration settings required by the underlying execution platform. This enables the workflow to be portable, allowing it to run on different computational platforms such as an institutional HPC or cloud infrastructure, without needing to modify the workflow implementation.

We have seen earlier that it is possible to provide a `process` with
directives. These directives are process specific configuration settings.
Similarly, we have also provided parameters to our workflow which
are parameter configuration settings. These configuration settings
can be separated from the workflow implementation, into a
configuration file.

## Configuration files

Settings in a configuration file are sets of name-value pairs
(`name = value`). The `name` is a specific property to set,
while the `value` can be anything you can assign to a variable (see
[nextflow scripting](02-nextflow_scripting)), for example, strings,
booleans, or other variables.
It is also possible to access any variable defined in the
host environment such as `$PATH`, `$HOME`, `$PWD`, etc.

~~~
// nextflow.config
my_home_dir = "$HOME"
~~~
{: .language-groovy}

> ## Accessing variables in your configuration file
>
> Generally, variables and functions defined in a
> configuration file are not accessible from the
> workflow script. Only variables defined using the
> `params` scope and the `env` scope can
> be accessed from the workflow script.
>
> ~~~
> workflow {
>     MY_PROCESS( params.input )
> }
> ~~~
> {: .language-groovy}
{: .callout}

Settings are also partitioned into
scopes, which govern the behaviour of different elements of the
workflow. For example, workflow parameters are governed from the
`params` scope, while process directives are governed from the `process` scope. A full list of the available scopes can be found in the
[documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes). It is also possible to define your own scope.

Configuration settings for a workflow are often stored in the file `nextflow.config` which is in the same directory as the workflow script. Configuration can be written in either of two ways. The first is using
dot notation, and the second is using brace notation. Both forms
of notation can be used in the same configuration file.

An example of dot notation:
~~~
params.input = ''             // The workflow parameter "input" is assigned an empty string to use as a default value
params.outdir = './results'   // The workflow parameter "outdir" is assigned the value './results' to use by default.
~~~
{: .language-groovy }

An example of brace notation:
~~~
params {
    input  = ''
    outdir = './results'
}
~~~
{: .language-groovy }

Configuration files can also be separated into multiple files and
included into another using the `includeConfig` statement.

~~~
// nextflow.config
params {
    input  = ''
    outdir = './results'
}

includeConfig 'system_resources.config'
~~~
{: .language-groovy}
~~~
// system_resources.config
process {
    cpus = 1    // default cpu usage
    time = '1h' // default time limit
}
~~~
{: .language-groovy}

## How configuration files are combined

Configuration settings can be spread across several files. This also
allows settings to be overridden by other configuration files. The
priority of a setting is determined by the following order,
ranked from highest to lowest.

1. Parameters specified on the command line (`--something value`).
1. Parameters provided using the `-params-file` option.
1. Config file specified using the `-c` my_config option.
1. The config file named `nextflow.config` in the current directory.
1. The config file named `nextflow.config` in the workflow project directory.
1. The config file `$HOME/.nextflow/config`.
1. Values defined within the workflow script itself (e.g., `main.nf`).

If configuration is provided by more than one of these methods,
configuration is merged giving higher priority to configuration
provided higher in the list.

Existing configuration can be completely ignored by using `-C <custom.config>` to use only configuration provided in the `custom.config` file.

> ## Configuring Nextflow vs Configuring a Nextflow workflow
>
> Parameters starting with a single dash `-` (e.g., `-c my_config.config`) are configuration
> options for `nextflow`, while parameters starting with a double
> dash `--` (e.g., `--outdir`) are workflow parameters defined in the `params` scope.
>
> The majority of Nextflow configuration settings must be provided
> on the command-line, however a handful of settings can also
> be provided within a configuration file, such as
> `workdir = '/path/to/work/dir' (`-w /path/to/work/dir`) or
> `resume = true` (`-resume`), and do not
> belong to a configuration scope.
{: .callout}

FIXME: Add exercise to test understanding of configuration priority.

## Configuring process behaviour

Earlier we saw that `process` directives allow the specification of
settings for the task execution such as `cpus`, `memory`, `conda`
and other resources in the pipeline script. This is useful when
prototyping a small workflow script, however this ties the configuration
to the workflow, making it less portable. A good practice is to
separate the process configuration settings into another file.

The `process` configuration scope allows the setting of any process directives in the Nextflow configuration file.

For example:

~~~
//nextflow.config
process {
    cpus = 2
    memory = 8.GB
    time = '1 hour'
    publishDir = [ path: params.outdir, mode: 'copy' ]
}
~~~
{: .language-groovy }

> ## Unit values
>
> Memory and time duration unit can be specified either using a string
> based notation in which the digit(s) and the unit can be separated by a > blank or
> by using the numeric notation in which the digit(s) and the unit are
> separated by a dot character and it’s not enclosed by quote characters.
>
>  | String syntax   | Numeric syntax | Value                 |
>  |-----------------|----------------|-----------------------|
>  | '10 KB'         | 10.KB          | 10240 bytes           |
>  | '500 MB'        | 500.MB         | 524288000 bytes       |
>  | '1 min'         | 1.min          | 60 seconds            |
>  | '1 hour 25 sec' | -              | 1 hour and 25 seconds |
>
{: .callout}

These settings are applied to all processes in the workflow. A
process selector can be used to apply the configuration to a
specific process or group of processes.

### Process selectors

The resources for a specific process can be defined using `withName:`
followed by the process name ( either the simple name e.g., `'FASTQC'`,
or the fully qualified name e.g., `'NFCORE_RNASEQ:RNA_SEQ:SAMTOOLS_SORT'`), and the directives within curly braces.
For example, we can specify different `cpus` and `memory` resources
for the processes `INDEX` and `FASTQC` as follows:

~~~
// process_resources.config
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
{: .language-groovy }

When a workflow has many processes, it is inconvenient to specify
directives for all processes individually, especially if directives
are repeated for groups of processes. A helpful strategy is to annotate
the processes using the `label` directive (processes can have multiple
labels). The `withLabel` selector then allows the configuration of all
processes annotated with a specific label, as shown below:

~~~
//configuration_process_labels.nf
nextflow.enable.dsl=2

process P1 {

    label "big_mem"

    script:
    """
    echo P1: Using $task.cpus cpus and $task.memory memory.
    """
}

process P2 {

    label "big_mem"

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
{: .language-groovy}

~~~
//configuration_process-labels.config
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
    }
}
~~~
{: .language-groovy}

Another strategy is to use process selector expressions. Both
`withName:` and `withLabel:` allow the use of regular expressions
to apply the same configuration to all processes matching a pattern.
Regular expressions must be quoted, unlike individual process names
or labels.

- The `|` matches either-or, e.g., `withName: 'INDEX|FASTQC'`
applies the configuration to any process matching the name `INDEX`
or `FASTQC`.
- The `!` inverts a selector, e.g., `withLabel: '!small_mem'` applies
the configuration to any process without the `small_mem` label.
- The `.*` matches any number of characters, e.g.,
`withName: 'NFCORE_RNASEQ:RNA_SEQ:BAM_SORT:.*'` matches all processes
of the workflow `NFCORE_RNASEQ:RNA_SEQ:BAM_SORT`.

A regular expression cheat-sheet can be found
[here](https://www.jrebel.com/system/files/regular-expressions-cheat-sheet.pdf) if you would like to
write more expressive expressions.

#### Selector priority

When mixing generic process configuration and selectors the following
priority rules are applied (from highest to lowest):

1. `withName` selector definition.
1. `withLabel` selector definition.
1. Process specific directive defined in the workflow script.
1. Process generic `process` configuration.

> ## Process selectors
> Create a Nextflow config, `process-selector.config`, file specifying different.
> 1. `cpus`
> 1. `memory`
>
> resources for the two process
> * `P1` (cpus 1 and memory 2.GB) and
> * `P2` (cpus 2 and memory 1.GB)
>
> for the Nextflow script `process-selector.nf`.
>
> ~~~
> //process-selector.nf
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
> {: .language-groovy }
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
> > ~~~
> > $ nextflow run process-selector.nf -c process-selector.config
> > ~~~
> > {: .language-bash }
> {: .solution}
 {: .challenge}

#### Dynamic expressions

A common scenario is that configuration settings may depend on the
data being processed. Such settings can be dynamically expressed
using a closure. For example, we can specify the `memory` required
as a multiple of the number of `cpus`. Similarly, we can publish
results to a subfolder based on the sample name.

~~~
process FASTQC {

    input:
    tuple val(sample), path(reads)

    script:
    """
    fastqc -t $task.cpus $reads
    """
}
~~~
{: .language-groovy }

~~~
//nextflow.config
process {
    withName: FASTQC {
        cpus = 2
        memory = { 2.GB * task.cpus }
        publishDir = { "fastqc/$sample" }
    }
}
~~~
{: .language-groovy }

> ## Configure process scope
> 1. Create a Nextflow config file `wc-params.config`
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

## Configuring execution platforms

Nextflow supports a wide range of execution platforms, from
running locally, to running on HPC clusters or cloud infrastructures.
See https://www.nextflow.io/docs/latest/executor.html for the
full list of supported executors.

![nf-executors](https://seqera.io/training/img/nf-executors.png)

The default executor configuration is defined within the `executor`
scope (https://www.nextflow.io/docs/latest/config.html#scope-executor).
For example, in the config below we specify the executor as
Sun Grid Engine, `sge` and the number of tasks the executor will
handle in a parallel manner (`queueSize`) to 10.

~~~
//nextflow.config
executor {
    name = 'sge'
    queueSize = 10
}
~~~
{: .language-groovy }

The `process.executor` directive allows you to override
the executor to be used by a specific process. This can be
useful, for example, when there are short running tasks
that can be run locally, and are unsuitable for submission
to HPC executors (check for guidelines on best practice use
of your execution system).

~~~
//nextflow.config
executor {
    name = 'sge'
    queueSize = 10
}
process {
    withLabel: 'short' {
        executor = 'local'
    }
}
~~~
{: .language-groovy }

## Defining software requirements


A key feature of Nextflow is the ability to use different technologies such as the Conda package management system or container engines such as Docker, Singularity, Podman, Charliecloud, Shifter to manage the software requirements. This also facilitate portable and reproducible workflows.

We can specify the technology used for the entire workflow or a single or group of process using the technology specific scopes.

## Config Conda execution

To use  a Conda environment you can either specify the path of an existing Conda environment directory.

~~~
process.conda = "/home/user/miniconda3/envs/my_conda_env"
//or
process {
  conda = "/home/user/miniconda3/envs/my_conda_env"
}
~~~
{: .source }

Specify the path of Conda environment `YAML` file.

For example;

~~~
process.conda = "environment.yml"
//or
process {
  conda = "environment.yml"  
}
~~~
{: .source }

Or, specify the an individual package name using the syntax.

~~~
<channel>::<package_name>=<version>
~~~


For example,

~~~
process.conda = "bioconda::salmon"
//or
process {
  conda = "bioconda::salmon"
}
~~~

### Conda scope

There is an optional `conda` scope allows you to control the creation of a Conda environment by the Conda package manager.

For example, `cacheDir` specifies the path  where the Conda environments are stored. By default this is in `conda` folder of the `work` directory.

**Note:** When using a compute cluster make sure to provide a shared file system path accessible from all computing nodes.

>## Define a software requirement in the configuration file using conda
>  Create a config file for the Nextflow script `configuration_fastp.nf` that:
>  1. Add a conda scope for the process name `NUM_LINES` that includes the bioconda package `fastp`.
>
>  Run the Nextflow script `configuration_fastp.nf` with the configuration file using the `-c` option.
>
> ~~~
> nextflow.enable.dsl=2
>
> params.input = "data/yeast/reads/ref1_1.fq.gz"
>
> process NUM_LINES {
>
>    input:
>    path read
>
>    output:
>    stdout
>
>    script:
>    """
>    printf '${read} '
>    gunzip -c ${read} | wc -l
>    fastp -i ${read}  -o out.fq 2>&1
>    """
> }
>
>
> input_ch = Channel.fromPath(params.input)
>
>
> workflow {
>    NUM_LINES(input_ch).out.view()
> }
>~~~
> {: .language-groovy }
> > ## Solution
> > ~~~
> > //fastp.config
> > process {
> >   withName: "NUM_LINES" {
> >     conda = "bioconda::fastp"
> >  }
> > }
> > ~~~
> > ~~~
> > $ nextflow run configuration_fastp.nf -c fastp.config
> > ~~~
> > {: .language-bash }
> {: .solution}
{: .challenge}

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
//configuration_profiles.config
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

## Inspecting the Nextflow configuration

You can use the command `nextflow config` to print the resolved
configuration of a workflow. This allows you to see what settings Nextflow will use to run a workflow.

~~~
$ nextflow config workflow_02.nf -profile test
~~~
{: .language-bash}

Would output:

~~~
FIXME: fill in
~~~
{: .output}

{% include links.md %}
