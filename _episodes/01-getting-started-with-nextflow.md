---
title: "Getting Started with Nextflow"
teaching: 30
exercises: 10
questions:
- "What is a workflow and what are workflow management systems?"
- "Why should I use a workflow management system?"
- "What is Nextflow?"
- "What are the main features of Nextflow?"
- "What are the main components of a Nextflow script?"
- "How do I run a Nextflow script?"
- "How can I use the nextflow logs?"
objectives:
- "Understand what a workflow management system is."
- "Understand the benefits of using a workflow management system."
- "Explain the benefits of using Nextflow as part of your bioinformatics workflow."
- "Explain the components of a Nextflow script."
- "Run a Nextflow script."
- "Use the nextflow log command to show information about executed pipelines."
- "Use the `-resume` option to execute the script using the cached results."
keypoints:
- "A workflow is a sequence of tasks that process a set of data, and a workflow management system (WfMS) is a computational platform that provides an infrastructure for the set-up, execution and monitoring of workflows."
- "Nextflow is a workflow management system that comprises both a runtime environment and a domain specific language (DSL)."
- "Nextflow scripts comprise of channels for controlling inputs and outputs, and processes for defining workflow tasks."
- "Nextflow stores working files in the work directory."
- "You run a Nextflow script using the `nextflow run` command."
- "You can resume a workflow, skipping cached steps, using the `-resume` option"
- "The `nextflow log` command can be used to see information about executed pipelines."
---


## Workflows

Analysing data involves a sequence of tasks, including gathering, cleaning and processing data. These sequence of tasks are called a workflow or a pipeline. These workflows typically require executing multiple software packages, sometimes running on different computing environments, such as a desktop or a compute cluster. Traditionally these workflows have been joined together in scripts using general purpose programming languages such as Bash or Python.



<br>
<center>
<img src="../fig/sarek_workflow.png" width="200" height="400" >
<img src="../fig/flowchart.png" width="300" height="400" >
<br>
<em>Example bioinformatics variant calling workflow/pipeline diagram from [nf-core](https://nf-co.re/sarek) and simple RNA-Seq pipeline in DAG format.
</em>
</center>
<br>


However, as workflows become larger and more complex, the management of the programming logic and software becomes difficult.

##  Workflow management systems

Recently Workflow Management Systems (WfMS), such as Snakemake, Galaxy, and Nextflow have emerged specifically to manage computational data-analysis workflows in fields such as Bioinformatics, Imaging, Physics, and Chemistry.  

These *Workflow management systems* contain multiple features that simplify the development, monitoring, execution and sharing of pipelines.

Key features include;

* **Run time management**: Management of program execution on the operating system and splitting tasks and data up to run at the same time in a process called parallelisation.
* **Software management**: Use of software management technology like containers, such as docker or singularity, that packages up code and all its dependencies so the application runs reliably from one computing environment to another.
* **Portability & Interoperability**: Workflows written on one system can be run on another computing infrastructure e.g. local computer, compute cluster, or cloud infrastructure.
* **Reproducibility**: The use of Software management systems and a pipeline specification means that the workflow will produce the same results when re-run, including on different computing platforms.
* **Reentrancy**: Continuous checkpoints allow workflows to resume
from the last successfully executed steps.

## Nextflow Basic concepts

Nextflow is a workflow management system that combines a runtime environment, software that is designed to run other software, and a *programming domain specific language (DSL)* that eases the writing of computational pipelines.

Nextflow is built around the idea that Linux is the lingua franca of data science. Nextflow follows Linux "small pieces loosely joined" philosophy : in which many simple but powerful command-line and scripting tools that, when chained together, facilitate more complex data manipulations.

Nextflow extends this approach, adding the ability to define complex program interactions and an accessible (high-level) parallel computational environment based on the [dataflow programming model](https://devopedia.org/dataflow-programming), whereby the processes are connected via their `outputs` and `inputs` to other `processes`, and processes run as soon as they receive an input.

The diagram below illustrates the differences between a dataflow model and a simple linear program .



<br>
<center>
<img src="../fig/dataflow.png" width="420" height="300" >
<br>
<em>A simple program (a) and its dataflow equivalent (b) https://doi.org/10.1145/1013208.1013209.
</em>
</center>
<br>

In a simple program **(a)**, these statements would be executed sequentially. Thus, the program would execute in three units of time. In the dataflow programming model **(b)**, this program takes only two units of time. This is because the read quantitation and QC steps have no dependencies on each other and therefore can execute simultaneously in parallel.

### Nextflow core features are:

1. Fast prototyping: A simple syntax for writing pipelines that enables you to reuse existing scripts and tools for fast prototyping.

1. Reproducibility: Nextflow supports several container technologies, such as Docker and Singularity, as well as the package manager conda. This, along with the integration of the GitHub code sharing platform, allows you to write self-contained pipelines, manage versions and to reproduce any former configuration.

1. Portability: Nextflow's syntax separates the functional logic (the steps of the workflow) from the execution settings (how the workflow is executed). This allows the pipeline to be run on multiple platforms, e.g. local compute vs. a university compute cluster or a cloud service like AWS, without changing the steps of the workflow.  

1. Simple parallelism:  Nextflow is based on the dataflow programming model which greatly simplifies the splitting of tasks that can be run at the same time (parallelisation).

1. Continuous checkpoints: All the intermediate results produced during the pipeline execution are automatically tracked. This allows you to resume its execution, from the last successfully executed step, no matter what the reason was for it stopping.

### Processes and Channels

 In practice a Nextflow pipeline is a script made by joining together different commands tasks in process blocks. Each process can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.).

Processes create a task for each complete input set. Each task is executed independently, and can not interact with another task. The only way data can be passed between processes is via asynchronous queues, called `channels` in Nextflow.

Processes uses these channels to define inputs and outputs. The interaction between  processes, and ultimately the pipeline execution flow itself, is implicitly defined by these input and output declarations.


Here we have a channel containing three elements, e.g. 3 data files. We have a process that takes the channel as input. The fact that the channel has three elements would mean that three independent instances (tasks) of that process are being run in parallel. The tasks then generate three outputs, that are used as input for another process.

<p align="center">
   <img alt="Processes and channels" src="../fig/channel-process_fqc.png" width="500">
   <br>
   <em>Nextflow process flow diagram</em>
</p>


### Workflow Execution

While a `process` defines what command or script has to be executed, the `executor` determines how that script is actually run in the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

<p align="center">
   <img alt="Processes and channels" src="../fig/executor.png" width="250">
   <br>
   <em>Nextflow Executors</em>
</p>


In this way Nextflow provides a separation between the pipeline’s functional logic and the underlying execution platform. This makes it possible to write a pipeline once, and then run it on your computer, compute cluster, or the cloud, without modifying the workflow, by simply defining the target execution platform in the configuration file.

Nextflow provides out-of-the-box support for major batch schedulers and cloud platforms such as
Sun Grid Engine, SLURM job scheduler, AWS Batch service and Kubernetes. A full list can be found [here](https://www.nextflow.io/docs/latest/executor.html).


### Scripting language

Nextflow scripts are written using a scripting language that simplifies the writing of workflows. Languages that are written for a specific fields are called Domain Specific Languages (DSL), e.g., SQL is used to work with databases, and AWK is designed for text processing.

In practical terms the Nextflow scripting language is an extension of the [Groovy programming language](https://groovy-lang.org/), which in turn is a super-set of the Java programming language. Groovy simplifies the writing of code and is more approachable than Java. Groovy semantics
(syntax, control structures, etc) are documented [here](https://groovy-lang.org/semantics.html).

The approach of having a simple DSL built on top of a more powerful general purpose programming language makes Nextflow very flexible. The Nextflow syntax can handle most workflow use cases with ease, and then Groovy can be used to handle corner cases which may be difficult to implement using the DSL.


## Your first script

We are now going to look at a simple Nextflow script that counts the number of lines in a file.

Open the file `wc.nf` in the script directory with your favourite text editor.

This is a Nextflow script. It contains;

1. A Shebang line, specifying the location of the Nextflow interpreter
1. A multi-line Nextflow comment, written using C style block comments.
1. A pipeline parameter `params.samples` which contains the relative path to the location of a fastq file.
1. A Nextflow channel created using the variable `params.samples` as input.
1. A Nextflow process block named `numLines`.
1. An input definition block that takes the samples channel `samples_ch` as input.
1. A script block that contains the bash commands `sleep` and `wc -l`
1. An output definition block that uses the Linux/Unix standard output stream `stdout` from the script block to create an output channel, named `read_out_ch`.
1. Finally the script writes the result to the terminal using the view operator.

~~~
#!/usr/bin/env nextflow

/*
*  Comment that starts with a slash asterisk /* and finishes with an asterisk slash  and you can place it anywhere in your code, on the same line or several lines.
*/

/*
 * pipeline input parameters
 */
params.samples  = "data/ggal/gut_1.fq"

/*
* Sample input channel
*/
samples_ch = Channel.fromPath(params.samples)


/*
* Nextflow process block
*/
process numLines {

    input:
    path read from samples_ch

    output:
    stdout into read_out_ch

    script:
    """
    sleep 5
    wc -l ${read}
    """
}

/*
* Nextflow view operator for showing contents of a channel
*/
read_out_ch.view()

~~~~
{: .language-groovy}

To run a nextflow script use the command `nextflow run <script_name>`

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
> {: .solution}
{: .challenge}

1. The first line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
1. The second line tells you how the process has been executed `executor >  local`.
1. The third line shows the process id `21/b259be`, name, number of cpus, percentage task completion and how many it times it has been run.
1. The final line is the result of the view operator.


## Pipeline parameters

The Nextflow `wc.nf` script defines a pipeline parameter `params.samples`. Pipeline parameters enable you to change the input to the workflow at runtime via the command line or a config file so they are not hard-coded into the script. In this way you change the input data or pipeline execution e.g., change the fastq file for input.

Pipeline parameters are declared by prepending the prefix `params`, separated by dot character, to a variable name e.g., `params.samples`. Their value can be specified on the command line by prefixing the parameter name with a **double dash** character, i.e. `--paramName` e.g. `--samples`

We can changed the input using the `params` variable on the command line.

> ## Add a pipeline parameter
> Re-run the Nextflow script by entering the following command in your terminal:
>
> ~~~
> nextflow run wc.nf --samples 'data/ggal/*.fq'
> ~~~
> {: .language-bash}
> > ## Solution
> > The string specified on the command line will override the default value of the parameter in the script. The output will look like this:
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
> {: .solution}
{: .challenge}

The pipeline has now executed the `numLines` process six time using the string `data/ggal/*.fq` to capture the six fastq files matching the pattern `data/ggal/*.fq`.


It’s worth noting that the process `wc -l` is executed in parallel, so there’s no guarantee on the output order. So, it is perfectly possible that you will get the final result printed out in a different order:


> ## Process identification
> The hexadecimal numbers, like b3/c9f4ee, identify the unique process execution.
> These numbers are also the prefix of the directories where each process is executed.
> You can inspect the files produced by them changing to the directory `$PWD/work` and
> using these numbers to find the process-specific execution path.
{: .callout}

## Nextflow log


Once a script has run Nextflow stores a log of all the pipelines executed in the current folder.
Similar to an electronic lab book this means you have a have a record of all processing steps and command run.

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
> {: .solution}
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
> {: .solution}
{: .challenge}

{: .source}


You will see that the execution of the process `numLines` is actually skipped (cached text appears), and its results are retrieved from the cache.


> ## Modify the wc.nf script and re-run the pipeline using -resume option
> Modify the wc.nf script changing the sleep time and execute the script by entering the following command in your terminal:
>
> ~~~
> nextflow run wc.nf --samples 'data/ggal.*.fq' -resume
> ~~~
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [backstabbing_joliot] - revision: 714e17a273
> > executor >  local (6)
> > [68/1ba655] process > numLines (5) [100%] 6 of 6 ✔
> >    11748 liver_1.fq
> >
> >    11748 lung_1.fq
> >
> >    11748 gut_2.fq
> >  
> >    11748 liver_2.fq
> >
> >    11748 lung_2.fq
> >
> >    11748 gut_1.fq
> >  ~~~
> {: .solution}
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

By default the pipeline results are cached by default in the directory `work` where the pipeline is launched.

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
