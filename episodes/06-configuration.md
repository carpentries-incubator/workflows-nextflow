---
title: "Nextflow configuration"
teaching: 30
exercises: 15
questions:
- "How can I control how Nextflow runs?"
- "How can I write a Nextflow configuration file?"
objectives:
- "Understand how Nextflow is configured."
- "Create a Nextflow configuration file."
keypoints:
- "Nextflow configuration can be managed using a `nextflow.config` file."
- "Nextflow configuration are simple text files containing a set of properties defined using the syntax."
---


## Nextflow configuration

A key Nextflow feature is the ability to decouple the workflow implementation by the configuration setting required by the underlying execution platform.

This enable portable deployment without the need to modify the application code.

## Configuration file

When a pipeline script is launched Nextflow looks for a file named `nextflow.config` in the current directory and in the script base directory
(if it is not the same as the current directory). Finally it checks for the file `$HOME/.nextflow/config`.

When more than one on the above files exist they are merged, so that the settings in the first override the same ones that may appear in the second one,
and so on.

The default config file search mechanism can be extended proving an extra configuration file by using the command line option `-c <config file>`.

## Config syntax

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax:

`name = value`

String values need to be wrapped in quotation characters while numbers and boolean values (`true`, `false`) do not.
Also note that values are typed, meaning for example that, `1` is different from `'1'`,
since the first is interpreted as the number one, while the latter is interpreted as a string value.

### Config variables

Configuration properties can be used as variables in the configuration file itself, by using the usual `$propertyName` or `${expression}` syntax. These variables are not available in the Nextflow script.

~~~
kmer = 27
kmer_message = "kmer size is  ${kmer}"
~~~
{: .language-groovy }

In the configuration file it’s possible to access any variable defined in the host environment such as `$PATH`, `$HOME`, `$PWD`, etc.

~~~
my_home_dir = "$HOME"
~~~
{: .language-groovy }

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

In programming Scope is a concept that refers to where values and functions can be accessed. When you define a variable in a Nextflow config file it is accessible via the main script.

Configuration settings can be organised in different scopes by dot prefixing `.` the property names with a scope identifier or grouping the properties in the same scope using the curly brackets notation `{}`. This is shown in the following example:

~~~
aligner.name  = "salmon"
aligner.kmer  = 27

index {
    kmer = 27
    outdir = 'results/index'
}
~~~
{: .language-groovy }

### Config params

The scope `param`s allows the definition of workflow parameters that overrides the values defined in the main workflow script.

This is useful to consolidate one or more execution parameters in a separate file.

~~~
// config file nextflow.config
params.kmer= 27


// workflow script
params.kmer = 31

println "$params.kmer"
~~~~
{: .language-groovy }

This would output

~~~
27
~~~
{: .output }

> ## Configuration parameters
> Save the snippet below in the file `nextflow.config`
> ~~~
> // config file
> params.genome = 'GRCh38'
> params.aligner = 'salmon'
> ~~~
> and then save the code below in the Nextflow script `params.nf`.
> ~~~
> // workflow script
> params.genome = 'GRCh37'
> params.aligner = 'kallisto'
> println "$params.genome $params.aligner"
> ~~~
> Then run:
> ~~~
> nextflow run params.nf
> ~~~
> {: .bash-language}
> Execute is again specifying the foo parameter on the command line:
>
> ~~~
> nextflow run params.nf --genome hg38
> ~~~
{: .bash-language}
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


### Config env

The `env` scope allows the definition one or more variable that will be exported in the environment where the workflow tasks will be executed.

Simply prefix your variable names with the env scope or surround them by curly brackets, as shown below:

~~~~
// configuration file
env.genome = 'hg38'
env {
   aligner = 'salmon'  
}
~~~~
{: .language-groovy }

~~~
// workflow script
script:
"""
echo \$genome
"""
~~~
{: .language-groovy }

> env scope
> ~~~
> env.kmer = '21'
> env.genome = "hg38"
> ~~~
> {: .language-groovy }
> Save the above snippet a file named `my-env.config`. Then save the snippet below in a file named `my-env.nf`:
> ~~~
> process envtest {
>  echo true
>  '''
>  echo  \$kmer
>  env | egrep 'genome'
>  '''
> }
>~~~
> {: .language-groovy }
>
> Finally executed the following command:
>
> ~~~
> nextflow run my-env.nf -c my-env.config
> ~~~~
> {: .bash-language}
> > Solution
> > This will print
> > ~~~~
> > 21
> > genome=hg38
> > ~~~
> > > {: .output}
> {: .solution}
> {: .challenge}


### Config process

The `process` directives allow the specification of specific settings for the task execution such as `cpus`, `memory`, `conda` and other resources in the pipeline script.

This is useful specially when prototyping a small workflow script.

However it’s always a good practice to decouple the workflow execution logic from the process configuration settings,
i.e. it’s strongly suggested to define the process settings in the workflow configuration `nextflow.config` file instead of the workflow script.
The process configuration scope allows the setting of any process directives in the Nextflow configuration file. For example:

~~~
// configuration file
process {
    cpus = 2
    memory = 8.GB
    container = 'biocontainers/bamtools:v2.4.0_cv3'
}
~~~
{: .language-groovy }

The above config snippet defines the `cpus`, `memory` and `container` directives for **all** processes in your workflow script.

The `process` selector can be used to apply the configuration to a specific process or group of processes (discussed later).


> Unit
> Memory and time duration unit can be specified either using a string based notation in which the digit(s) and the unit can be separated by a blank or
by using the numeric notation in which the digit(s) and the unit are separated by a dot character and it’s not enclosed by quote characters.
{: .callout}

| String syntax   | Numeric syntax | Value                 |
|-----------------|----------------|-----------------------|
| '10 KB'         | 10.KB          | 10240 bytes           |
| '500 MB'        | 500.MB         | 524288000 bytes       |
| '1 min'         | 1.min          | 60 seconds            |
| '1 hour 25 sec' | -              | 1 hour and 25 seconds |

> process directives require `=` in configuration file.
>The syntax for setting process directives in the configuration file requires = ie. assignment  operator, instead it should not be used when setting process directives in the workflow script.
{: .callout}

This is important especially when you want to define a config setting using a dynamic expression using a closure. For example:

~~~
process {
    memory = { 4.GB * task.cpus }
}
~~~
{: .language-groovy }


Finally directives that allows to be repeated in the process definition, in the configuration files need to be defined as a list object. For example:

~~~
process {
    pod = [ [env: 'FOO', value: '123'],
            [env: 'BAR', value: '456'] ]
}
~~~
{: .language-groovy }



### Config Docker execution

The container image to be used for the process execution can be specified in the `nextflow.config` file:

~~~
process.container = 'nextflow/rnaseq-nf'
docker.enabled = true
~~~
{: .language-groovy }

The use of the unique SHA256 image ID guarantees that the image content do not change over time

~~~
process.container = 'nextflow/rnaseq-nf@sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266'
docker.enabled = true
~~~

{: .source}

### Config Singularity execution

The run the workflow execution with a Singularity container provide the container image file path in the Nextflow config file using the container directive:

~~~
process.container = '/some/singularity/image.sif'
singularity.enabled = true
~~~
{: .language-groovy }

The container image file must be an absolute path i.e. it must start with a `/`.

The following protocols are supported:

* `library://`` download the container image from the Singularity Library service.

* `shub://`` download the container image from the Singularity Hub.

* `docker://`` download the container image from the Docker Hub and convert it to the Singularity format.

* `docker-daemon://` pull the container image from a local Docker installation and convert it to a Singularity image file.

> Docker to Singularity image
Specifying a plain Docker container image name, Nextflow implicitly download and converts it to a Singularity image when the Singularity execution is enabled.
For example:
> ~~~
> process.container = 'nextflow/rnaseq-nf'
> singularity.enabled = true
> ~~~
> {: .source}
{: .callout}

The above configuration instructs Nextflow to use Singularity engine to run your script processes.
The container is pulled from the Docker registry and cached in the current directory to be used for further runs.

Alternatively if you have a Singularity image file, its location absolute path can be specified as the container name either using the `-with-singularity` option or the `process.container` setting in the config file.

Try to run the script as shown below:

~~~
nextflow run script7.nf
~~~
{: .bash-language}

Note: Nextflow will pull the container image automatically, it will require a few seconds depending the network connection speed.

### Config Conda execution

The use of a Conda environment can also be provided in the configuration file adding the following setting in the `nextflow.config` file:

~~~
process.conda = "/home/ubuntu/miniconda3/envs/nf-tutorial"
~~~
{: .language-groovy }

You can either specify the path of an existing Conda environment directory or the path of Conda environment `YAML` file.



{% include links.md %}
