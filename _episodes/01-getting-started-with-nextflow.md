---
title: "Getting Started with Nextflow"
teaching: 30
exercises: 10
questions:
- "What is workflow?"
- "What is workflow management system?"
- "Why should I use a workflow management system?"
- "What is Nextflow?"
- "What are the main features of Nextflow?"
- "What are the main components of a Nextflow script?"
- "How do you run a Nextflow script?"
objectives:
- "Understand what a workflow management system is."
- "Understand the benefits of using a workflow management system."
- "Explain the benefits of using Nextflow as part of your science workflow."
- "Explain the components of a Nextflow script."
- "Run a simple Nextflow script."
keypoints:
- "A workflow is sequence of tasks that processes a set of data. "
- "A workflow management system is a computational platform that manages the execution of a workflow."
- "Nextflow is a workflow management system that comprises both a runtime environment and a domain specific language (DSL)."
- "Using a workflow system facilitates portability and reproducibility of workflows."
- "Nextflow scripts comprises of channels for controlling inputs and outputs and processes for defining workflow tasks."
- "You run a Nextflow script using the `nextflow run` command."
---


## Workflows and Workflow management systems

Analysing data involves a sequence of tasks, including, gathering, cleaning and processing data. These sequence of tasks is  generally called a workflow or a pipeline. These workflows typically requires multiple software and packages, sometimes in running in different computer environments. Traditionally these steps have been joined together in scripts using general purpose programming languages such as Bash or Python.

 *Workflow management  systems* (WfMS)  have  emerged specifically catering to computational data-analysis  in field such as Bioinformatics, Medical Imaging, Astronomy, Physics, and Chemistry.  

These *Workflow management systems* contain multiple features that aid in the development, monitoring,  execution and sharing of pipelines.

Key features include;

* Run time management: Management of program execution and task parallelisation.
* Software management: Use of container and environment managers to manage software dependencies.
* Portability & Interoperability: Bioinformatic workflow written on one system can be run on another computing infrastructure e.g. local vs HPC.
* Reproducibility: Workflow produces the same results when run on different platforms.


## Nextflow Basic concepts

Nextflow is a *workflow management system* that combines a runtime environment (software that is designed to run other software) and a *programming domain specific language (DSL)* that eases the writing of data-intensive computational pipelines.

It is designed around the idea that the Linux platform is the *lingua franca* of data science. Linux provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

Nextflow extends this approach, adding the ability to define complex program interactions and a high-level parallel computational environment based on the [dataflow programming model](https://devopedia.org/dataflow-programming) whereby the processes are connected via their `outputs` and `inputs` to other `processes`, and processes run as soon as they receive an input.


Nextflow core features are:

1. A custom domain specific language (DSL) for writing pipelines that enables fast prototyping.

1. Enable workflows portability & reproducibility: Nextflow's syntax is separated out from where   the pipeline is deployed e.g. local compute vs HPC.  

1. Simplify parallelisation and large scale deployment: The dataflow programming model enables implicit parallelism.

1. Easily integrate existing tools, systems & industry standards.

note <https://www.youtube.com/watch?v=8_i8Tn335X0&ab_channel=Nextflow> minute 20:15

### Processes and Channels

In practice a Nextflow pipeline script is made by joining together different processes (workflow steps). Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state. The only way they can communicate is via asynchronous FIFO queues, called `channels` in Nextflow.

Any process can define one or more channels as input and output. The interaction between these processes, and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.


Here we have a channel containing three elements, e.g. 3 data files, . We have a process that takes the channel as input. The fact that the channels has three elements would mean that three independent instances of that process are being run in parallel. The processes would generate three outputs.

<p align="center">
   <img alt="Processes and channels" src="../fig/channel-process.png" width="500">
</p>


### Execution abstraction

While a `process` defines what command or script has to be executed, the `executor` determines how that script is actually run in the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, Nextflow provides an abstraction between the pipeline’s functional logic and the underlying execution system. Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud, without modifying it, by simply defining the target execution platform in the configuration file.

It provides out-of-the-box support for major batch schedulers and cloud platforms:

|Name|Executor|
|----|--------|
|local|The process is executed in the computer where Nextflow is launched.|
|sge|The process is executed using the Sun Grid Engine / Open Grid Engine.|
|uge|The process is executed using the Univa Grid Engine job scheduler.|
|lsf|The process is executed using the Platform LSF job scheduler.|
|slurm|The process is executed using the SLURM job scheduler.|
|pbs|The process is executed using the PBS/Torque job scheduler.|
|pbspro|The process is executed using the PBS Pro job scheduler.|
|moab|The process is executed using the Moab job scheduler.|
|condor|The process is executed using the HTCondor job scheduler.|
|nqsii|The process is executed using the NQSII job scheduler.|
|ignite|The process is executed using the Apache Ignite cluster.|
|k8s|The process is executed using the Kubernetes cluster.|
|awsbatch|The process is executed using the AWS Batch service.|
|google-pipelines|The process is executed using the Google Genomics Pipelines service.|

### Scripting language

Nextflow implements declarative domain specific language (DSL) that simplifies the syntax of writing complex data analysis workflows as an extension of a general purpose programming language.

This approach makes Nextflow very flexible because allows to have in the same computing environment the benefit of concise DSL that allow the handling of recurrent use cases with ease and the flexibility and power of a general purpose programming language to handle corner cases, which may be difficult to implement using a declarative approach.

In practical terms Nextflow scripting is an extension of the [Groovy programming language](https://groovy-lang.org/), which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java in that is simplifies the writing of code and is more approachable.

## Your first script

Copy the following example into your favourite text editor and save it to a file named `wc.nf` :

~~~
#!/usr/bin/env nextflow

//pipeline parameter
params.samples  = "data/ggal/gut_1.fq"
samples_ch = Channel.fromPath(params.samples)


process numLines {

    input:
    path read from samples_ch

    output:
    stdout into read_out_ch

    """
    sleep 5
    wc -l ${read}
    """
}

read_out_ch.view()

~~~~
{: .source}

This script defines one `process` that counts the number of lines in a single fastq file. The output, as captured by stdout qualifier, is emitted on the result `channel` and the final output is printed by the `view` operator.


> ## Run a Nextflow  script
> Run the script by entering the following command in your terminal:
>
> ~~~
> nextflow run wc.nf
> ~~~
> {: .language-bash}
> > ## Solution
> > It will output something similar to the text shown below:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [fervent_babbage] - revision: c54a707593
> > executor >  local (1)
> > [21/b259be] process > numLines (1) [100%] 1 of 1 ✔
> >   11748 gut_1.fq
> >  ~~~
{: .challenge}

You can see that the  process `numLines` is executed once and the result string is printed.


## Pipeline parameters

The Nextflow `wc.nf` script defines a pipeline parameter `params.samples`.
Pipeline parameters are simply declared by prepending to a variable name the prefix `params`, separated by dot character e.g. `params.reads`. Their value can be specified on the command line by prefixing the parameter name with a **double dash** character, i.e. `--paramName` e.g. `--reads`

We can changed the input using the `params` variable on the command line.

> ## Add a pipeline parameter
> Run the Netxflow script by entering the following command in your terminal:
>
> ~~~
> nextflow run wc.nf --samples 'data/ggal/*.fq'
> ~~~
> {: .language-bash}
> > ## Solution
> > The string specified on the command line will override the default value of the parameter. The output will look like this:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [soggy_miescher] - revision: c54a707593
> > executor >  local (6)
> > [b3/c9f4ee] process > numLines (1) [100%] 6 of 6 ✔
> >    11748 gut_2.fq
> >
> >    11748 liver_2.fq
> >
> >    11748 gut_1.fq
> >
> >    11748 lung_1.fq
> >
> >    11748 liver_1.fq
> >
> >    11748 lung_2.fq
> >  ~~~
{: .challenge}

The pipeline has now executed the `numLines` process six time using the string `data/ggal/*.fq` to capture the six fastq files matching the glob pattern `data/ggal/*.fq`.


It’s worth noting that the process `wc` is executed in parallel, so there’s no guarantee on the output order.

Thus, it is perfectly possible that you will get the final result printed out in a different order:


> ## Process identification
> The hexadecimal numbers, like b3/c9f4ee, identify the unique process execution.
> These numbers are also the prefix of the directories where each process is executed.
> You can inspect the files produced by them changing to the directory `$PWD/work` and
> using these numbers to find the process-specific execution path.
{: .callout}

## Nextflow log

You can print Nextflow's execution history and log information using the  `nextflow log` command.

## Show Execution Log
> Listing the execution logs of previous invocations of all pipelines in a directory.
>
> ~~~
> nextflow log
> ~~~
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
> >TIMESTAMP          	DURATION	RUN NAME       	STATUS	REVISION ID	SESSION ID                          	COMMAND
> >2021-03-19 13:45:53	6.5s    	fervent_babbage	OK    	c54a707593 	15487395-443a-4835-9198-229f6ad7a7fd	nextflow run wc.nf
> > 2021-03-19 13:46:53	6.6s    	soggy_miescher 	OK    	c54a707593 	58da0ccf-63f9-42e4-ba4b-1c349348ece5	nextflow run wc.nf --samples 'data/ggal/*.fq'
> >  ~~~
{: .challenge}


## Modify and resume

Nextflow keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.


> ## Re-run the pipeline using -resume option
> Execute the script by entering the following command in your terminal:
>
> ~~~
> nextflow run wc.nf --samples 'data/ggal/*.fq' -resume
> ~~~
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [crazy_sax] - revision: c54a707593
> > [0a/0425b3] process > numLines (3) [100%] 6 of 6, cached: 6 ✔
> >    11748 lung_1.fq
> >
> >    11748 liver_2.fq
> >
> >    11748 liver_1.fq
> >
> >    11748 lung_2.fq
> >
> >    11748 gut_2.fq
> >
> >    11748 gut_1.fq
> >  ~~~
{: .challenge}

{: .source}


You will see that the execution of the process `numLines` is actually skipped (cached text appears), and its results are retrieved from the cache.


> ## Modify the wc.nf script and re-run the pipeline using -resume option
> Modify the wc.nf script changing the sleep time and execute the script by entering the following command in your terminal:
>
> ~~~
> nextflow run myfirst.nf --samples 'data/ggal.*.fq' -resume
> ~~~
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
N E X T F L O W  ~  version 20.10.0
Launching `wc.nf` [backstabbing_joliot] - revision: 714e17a273
executor >  local (6)
[68/1ba655] process > numLines (5) [100%] 6 of 6 ✔
   11748 liver_1.fq

   11748 lung_1.fq

   11748 gut_2.fq

   11748 liver_2.fq

   11748 lung_2.fq

   11748 gut_1.fq
> >  ~~~
{: .challenge}

As you have changed the script the pipeline will re-run and won't use the cached results for that process.


~~~
nextflow log
~~~
{: .language-bash}

~~~
IMESTAMP          	DURATION	RUN NAME           	STATUS	REVISION ID	SESSION ID                          	COMMAND
2021-03-19 13:45:53	6.5s    	fervent_babbage    	OK    	c54a707593 	15487395-443a-4835-9198-229f6ad7a7fd	nextflow run wc.nf
2021-03-19 13:46:53	6.6s    	soggy_miescher     	OK    	c54a707593 	58da0ccf-63f9-42e4-ba4b-1c349348ece5	nextflow run wc.nf --samples 'data/ggal/*.fq'
2021-03-19 13:51:40	2s      	crazy_sax          	OK    	c54a707593 	58da0ccf-63f9-42e4-ba4b-1c349348ece5	nextflow run wc.nf --samples 'data/ggal/*.fq' -resume
2021-03-19 13:52:15	4.5s    	backstabbing_joliot	OK    	714e17a273 	58da0ccf-63f9-42e4-ba4b-1c349348ece5	nextflow run wc.nf --samples 'data/ggal/*.fq' -resume
~~~
{: .output}
## work directory

The pipeline results are cached by default in the directory `work` where the pipeline is launched.

~~~
work/
├── 0a
│   └── 0425b3789460591b49ef8dae6ec44b
│       └── gut_2.fq -> /Users/ggrimes2/Documents/nf/data/ggal/gut_2.fq
├── 21
│   └── b259bed0f7dacf8b83bcc0da5a079f
│       └── gut_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/gut_1.fq
├── 48
│   └── fe7f1b6c838cfeeaf022b08b1f3fc8
│       └── liver_2.fq -> /Users/ggrimes2/Documents/nf/data/ggal/liver_2.fq
├── 4f
│   └── f773dcbb2045831171497637e31e88
│       └── gut_2.fq -> /Users/ggrimes2/Documents/nf/data/ggal/gut_2.fq
├── 54
│   └── 1cf9da3ac999d80b99b7dc4ac9b270
│       └── lung_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/lung_1.fq
├── 55
│   └── 3536aa774ef8e1a51e901ec9eb27bb
│       └── liver_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/liver_1.fq
├── 57
│   └── f851843afff11994536dfd0300090f
│       └── liver_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/liver_1.fq
├── 68
│   └── 1ba655e16066348420844c00c9c081
│       └── gut_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/gut_1.fq
├── 6d
│   ├── 211fd4c064bf5b7176865fa2ffeaca
│   │   └── lung_2.fq -> /Users/ggrimes2/Documents/nf/data/ggal/lung_2.fq
│   └── e0f493681b75fc2f9d6bfa35c6ec54
│       └── gut_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/gut_1.fq
├── 8b
│   └── b7e9351f3e7254249170294c3cd6eb
│       └── lung_1.fq -> /Users/ggrimes2/Documents/nf/data/ggal/lung_1.fq
├── b3
│   └── c9f4eeef85981d8569e4d8b7deb5a8
│       └── lung_2.fq -> /Users/ggrimes2/Documents/nf/data/ggal/lung_2.fq
└── de
    └── 62ee012039ea0af53e2898c4dc4bc8
        └── liver_2.fq -> /Users/ggrimes2/Documents/nf/data/ggal/liver_2.fq
~~~


Depending on your script, this folder can take of lot of disk space.
You can specify another work directory using the command line option `-w`
~~~
nextflow run <script> -w /some/scratch/dir
~~~
{: .language-bash}

If your are sure you won’t resume your pipeline execution, clean this folder periodically using the command `nextflow clean`.

~~~
nextflow clean [run_name|session_id] [options]
~~~
{: .language-bash}



{% include links.md %}
