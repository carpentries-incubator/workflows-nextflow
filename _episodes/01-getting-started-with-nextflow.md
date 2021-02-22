---
title: "Getting Started with Nextflow"
teaching: 15
exercises: 5
questions:
- "What is Nextflow?"
- "Why should I use a workflow system as part of my research workflow?"
- "Why use Nextflow ?"
objectives:
- "Understand why you should use a workflow system as part of your 
  (data) science workflow." 
- "Explain the benefits of using Nextflow as part of your (data) science workflow."
keypoints:
- "Nextflow is a "workflow orchestration engine and a programming domain specific language (DSL)."
- "Using a workflow system facilitates portability and reproducibility 
  of (data) science workflows."
- "Nextflow is the workflow language, processes can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.)"
---

## Basic concepts

Nextflow is workflow orchestration engine and a programming domain specific language (DSL) that eases the writing of data-intensive computational pipelines.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

Nextflow extends this approach, adding the ability to define complex program interactions and a high-level parallel computational environment based on the dataflow programming model. Nextflow core features are:

* enable workflows portability & reproducibility

* simplify parallelization and large scale deployment

* easily integrate existing tools, systems & industry standards

### Processes and channels

In practice a Nextflow pipeline script is made by joining together different processes. Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state. The only way they can communicate is via asynchronous FIFO queues, called channels in Nextflow.

Any process can define one or more channels as input and output. The interaction between these processes, and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.


> <p align="center">
>   <img alt="Processes and channels" src="../fig/channel-process.png" width="500">
> </p>
> 

### Execution abstraction

While a process defines what command or script has to be executed, the executor determines how that script is actually run in the target system.

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

Copy the following example into your favourite text editor and save it to a file named hello.nf :

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

    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

result.view{ it.trim() }

~~~~

This script defines two processes. The first splits a string into files containing chunks of 6 characters. The second receives these files and transforms their contents to uppercase letters. The resulting strings are emitted on the result channel and the final output is printed by the view operator.

Execute the script by entering the following command in your terminal:

~~~
nextflow run hello.nf
~~~

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


{% include links.md %}



