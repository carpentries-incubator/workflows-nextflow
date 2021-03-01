---
title: "Getting Started with Nextflow"
teaching: 20
exercises: 10
questions:
- "What is Nextflow?"
- "Why should I use a workflow system as part of my research workflow?"
objectives:
- "Understand why you should use a workflow system as part of your workflow."
- "Explain the benefits of using Nextflow as part of your science workflow."
keypoints:
- "Nextflow is a "workflow orchestration engine and a programming domain specific language."
- "Using a workflow system facilitates portability and reproducibility of workflows."
---


## Basic concepts

Nextflow is a workflow management system and a programming domain specific language (DSL) that eases the writing of data-intensive computational pipelines.
It is designed around the idea that the Linux platform is the *lingua franca* of data science. Linux provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

Nextflow extends this approach, adding the ability to define complex program interactions and a high-level parallel computational environment based on the [dataflow programming model](https://devopedia.org/dataflow-programming). Nextflow core features are:

1. Enable workflows portability & reproducibility.

1. Simplify parallelization and large scale deployment.

1. Easily integrate existing tools, systems & industry standards

### Processes and channels

In practice a Nextflow pipeline script is made by joining together different processes. Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state. The only way they can communicate is via asynchronous FIFO queues, called channels in Nextflow.

Any process can define one or more channels as input and output. The interaction between these processes, and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.


> <p align="center">
>   <img alt="Processes and channels" src="../fig/channel-process.png" width="500">
> </p>
> 

### Execution abstraction

While a `process` defines what command or script has to be executed, the `executor` determines how that script is actually run in the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, Nextflow provides an abstraction between the pipeline’s functional logic and the underlying execution system. Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud, without modifying it, by simply defining the target execution platform in the configuration file.

It provides out-of-the-box support for major batch schedulers and cloud platforms:

* Grid engine (Open/Sun/Univa)

* IBM Platform LSF

* Linux SLURM

* PBS Works

* Torque

* Moab

* HTCondor

* Amazon Batch

* Google Life Sciences

* Kubernetes

### Scripting language

Nextflow implements declarative domain specific language (DSL) simplifies the writing of writing complex data analysis workflows as an extension of a general purpose programming language.

This approach makes Nextflow very flexible because allows to have in the same computing environment the benefit of concise DSL that allow the handling of recurrent use cases with ease and the flexibility and power of a general purpose programming language to handle corner cases, which may be difficult to implement using a declarative approach.

In practical terms Nextflow scripting is an extension of the [Groovy programming language](https://groovy-lang.org/), which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java in that is simplifies the writing of code and is more approachable.

## Your first script

Copy the following example into your favourite text editor and save it to a file named `hello.nf` :

~~~
#!/usr/bin/env nextflow

params.greeting  = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

process splitLetters {

    input:
    val x from greeting_ch

    output:
    file 'chunk_*' into letters

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process convertToUpper {

    input:
    file y from letters.flatten()

    output:
    stdout into result
    
    script:
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

result.view{ it.trim() }

~~~~
{: .source}

This script defines two processes. The first splits a string into files containing chunks of 6 characters. The second receives these files and transforms their contents to uppercase letters. The resulting strings are emitted on the result channel and the final output is printed by the view operator.

Execute the script by entering the following command in your terminal:

~~~
nextflow run hello.nf
~~~
{: .source}

It will output something similar to the text shown below:

~~~
N E X T F L O W  ~  version 20.01.0
Launching `hello.nf` [marvelous_plateau] - revision: 63f8ad7155
[warm up] executor > local
executor >  local (3)
[19/c2f873] process > splitLetters   [100%] 1 of 1 ✔
[05/5ff9f6] process > convertToUpper [100%] 2 of 2 ✔
HELLO
WORLD!
~~~
{: .output}

You can see that the first process is executed once, and the second twice. Finally the result string is printed.

It’s worth noting that the process `convertToUpper` is executed in parallel, so there’s no guarantee that the instance processing the first split (the chunk Hello) will be executed before before the one processing the second split (the chunk world!).

Thus, it is perfectly possible that you will get the final result printed out in a different order:

~~~

WORLD!
HELLO
~~~
{: .output}

> ## Process identification
> The hexadecimal numbers, like 22/7548fa, identify the unique process execution. 
> These numbers are also the prefix of the directories where each process is executed. 
> You can inspect the files produced by them changing to the directory `$PWD/work` and 
> using these numbers to find the process-specific execution path.
{: .callout}

## Modify and resume

Nextflow keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the convertToUpper process in the previous example, replacing the process script with the string `rev $x`, so that the process looks like this:

~~~
process convertToUpper {

    input:
    file y from letters.flatten()

    output:
    stdout into result

    """
    rev $y
    """
}
~~~
{: .source}

Then save the file with the same name, and execute it by adding the `-resume` option to the command line:

~~~
nextflow run hello.nf -resume
~~~
{: .source}

It will print output similar to this:

~~~
N E X T F L O W  ~  version 20.01.0
Launching `hello.nf` [naughty_tuckerman] - revision: 22eaa07be4
[warm up] executor > local
executor >  local (2)
[19/c2f873] process > splitLetters   [100%] 1 of 1, cached: 1 ✔
[a7/a410d3] process > convertToUpper [100%] 2 of 2 ✔
olleH
!dlrow
~~~
{: .output}

You will see that the execution of the process splitLetters is actually skipped (the process ID is the same), and its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.


> ## work directory
> The pipeline results are cached by default in the directory $PWD/work. 
> Depending on your script, this folder can take of lot of disk space. 
> If your are sure you won’t resume your pipeline execution, clean this folder periodically.
{: .callout}

## Pipeline parameters

Pipeline parameters are simply declared by prepending to a variable name the prefix params, separated by dot character. Their value can be specified on the command line by prefixing the parameter name with a double dash character, i.e. `--paramName`

For the sake of this tutorial, you can try to execute the previous example specifying a different input string parameter, as shown below:


~~~
nextflow run hello.nf --greeting 'Bonjour le monde!'
~~~
{: .bash}

The string specified on the command line will override the default value of the parameter. The output will look like this:

~~~
N E X T F L O W  ~  version 20.01.0
Launching `hello.nf` [wise_stallman] - revision: 22eaa07be4
[warm up] executor > local
executor >  local (4)
[48/e8315b] process > splitLetters   [100%] 1 of 1 ✔
[01/840ca7] process > convertToUpper [100%] 3 of 3 ✔
uojnoB
m el r
!edno

~~~
{: .output}




{% include links.md %}
