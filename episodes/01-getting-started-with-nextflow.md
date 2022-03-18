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
objectives:
- "Understand what a workflow management system is."
- "Understand the benefits of using a workflow management system."
- "Explain the benefits of using Nextflow as part of your bioinformatics workflow."
- "Explain the components of a Nextflow script."
- "Run a Nextflow script."
keypoints:
- "A workflow is a sequence of tasks that process a set of data."
- "A workflow management system (WfMS) is a computational platform that provides an infrastructure for the set-up, execution and monitoring of workflows."
- "Nextflow is a workflow management system that comprises both a runtime environment and a domain specific language (DSL)."
- "Nextflow scripts comprise of channels for controlling inputs and outputs, and processes for defining workflow tasks."
- "You run a Nextflow script using the `nextflow run` command."
---

## Workflows

Analysing data involves a sequence of tasks, including gathering, cleaning, and processing data. These sequence of tasks are called a workflow or a pipeline. These workflows typically require executing multiple software packages, sometimes running on different computing environments, such as a desktop or a compute cluster. Traditionally these workflows have been joined together in scripts using general purpose programming languages such as Bash or Python.

<br>
<center>
<img src="https://sateeshperi.github.io/nextflow_varcal/nextflow/images/analysis_workflow.PNG" width="200">
<br>
<em> Example bioinformatics variant calling workflow/pipeline diagram from nf-core (https://nf-co.re/sarek) and simple RNA-Seq pipeline in DAG format. </em>
</center>
<br>

However, as workflows become larger and more complex, the management of the programming logic and software becomes difficult.

## Workflow management systems

*Workflow Management Systems* (WfMS) such as [Snakemake](https://snakemake.readthedocs.io/en/stable/),
[Galaxy](https://usegalaxy.org/), and [Nextflow](nextflow.io/) have been
developed specifically to manage computational data-analysis workflows in
fields such as bioinformatics, imaging, physics, and chemistry. These systems
contain multiple features that simplify the development, monitoring, execution
and sharing of pipelines:

* Run time management
* Software management
* Portability & Interoperability
* Reproducibility
* Re-entrancy

<br>
<center>
<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-021-01254-9/MediaObjects/41592_2021_1254_Fig1_HTML.png?as=webp" width="600">
<br>
<em>
An example of differences between running a specific analysis workflow using
a traditional pipeline or a WfMS-based pipeline. Source: Wratten, L., Wilm, A.
& Göke, J. Reproducible, scalable, and shareable analysis pipelines with
bioinformatics workflow managers. Nat Methods 18, 1161–1168 (2021).
https://doi.org/10.1038/s41592-021-01254-9
</em>
</center>
<br>

## Nextflow core features

<br>
<center>
<img src="https://training.seqera.io/img/execution_abstraction.png" width="600">
<br>
<em> Overview of Nextflow core features. </em>
</center>
<br>

* **Fast prototyping**: A simple syntax for writing pipelines that enables you
  to reuse existing scripts and tools for fast prototyping.

* **Reproducibility**: Nextflow supports several container technologies, such
  as [Docker](https://www.docker.com/) and [Singularity](https://sylabs.io/singularity),
  as well as the package manager [Conda](https://docs.conda.io). This, along
  with the integration of the [GitHub](https://www.github.com) code sharing
  platform, allows you to write self-contained pipelines, manage versions and
  to reproduce any previous result when re-run, including on different
  computing platforms.

* **Portability & interoperability**: Nextflow's syntax separates the
  functional logic (the steps of the workflow) from the execution settings (how
  the workflow is executed). This allows the pipeline to be run on multiple
  platforms, e.g. local compute vs. a university compute cluster or a cloud
  service like [AWS](https://aws.amazon.com/), without changing the steps of
  the workflow.

* **Simple parallelism**:  Nextflow is based on the dataflow programming model
  which greatly simplifies the splitting of tasks that can be run at the same
  time (parallelisation).

* **Continuous checkpoints & re-entrancy**: All the intermediate results
  produced during the pipeline execution are automatically tracked. This allows
  you to resume its execution from the last successfully executed step, no
  matter what the reason was for it stopping.

## Processes, channels, and workflows

Nextflow workflows have three main parts: *processes*, *channels*, and
*workflows*.

* *Processes* describe a task to be run. A process script can be written in any
  scripting language that can be executed by the Linux platform (Bash, Perl,
  Ruby, Python, R, etc.). Processes spawn a task for each complete input set.
  Each task is executed independently and cannot interact with other tasks. The
  only way data can be passed between process tasks is via asynchronous queues,
  called *channels*.

* Processes define inputs and outputs for a task. *Channels* are then used to
  manipulate the flow of data from one process to the next.

* The interaction between processes, and ultimately the pipeline execution flow
  itself, is then explicitly defined in a *workflow* section.

In the following example we have a channel containing three elements, e.g.,
three data files. We have a process that takes the channel as input. Since the
channel has three elements, three independent instances (tasks) of that process
are run in parallel. Each task generates an output, which is passed to another
channel and used as input for the next process.

<p align="center">
   <img alt="Processes and channels" src="../fig/channel-process_fqc.png" width="700">
   <br>
   <em> Nextflow process flow diagram. </em>
</p>

## Workflow execution

While a `process` defines what command or script has to be executed, the
`executor` determines how that script is actually run in the target system.

If not otherwise specified, processes are executed on the local computer. The
local executor is very useful for pipeline development, testing, and
small-scale workflows, but for large-scale computational pipelines, a High
Performance Cluster (HPC) or Cloud platform is often required.

<p align="center">
   <img alt="Processes and channels" src="../fig/executor.png" width="350">
   <br>
   <em>Nextflow Executors</em>
</p>

Nextflow provides a separation between the pipeline’s functional logic and the
underlying execution platform. This makes it possible to write a pipeline once,
and then run it on your computer, compute cluster, or the cloud, without
modifying the workflow, by defining the target execution platform in
a configuration file.

Nextflow provides out-of-the-box support for major batch schedulers and cloud
platforms such as Sun Grid Engine, SLURM job scheduler, AWS Batch service and
Kubernetes; a full list can be found [here](https://www.nextflow.io/docs/latest/executor.html).

## Your first script

We are now going to look at a sample Nextflow script that counts the number of
lines in a file. Create the file `word_count.nf` in the current directory using
your favourite text editor and copy-paste the following code:

~~~
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator.
========================================================================================
*/

params.input = "data/untrimmed_fastq/SRR2584863_1.fastq.gz"

/*
========================================================================================
    Input data is received through channels
========================================================================================
*/

input_ch = Channel.fromPath(params.input)

/*
========================================================================================
   Main Workflow
========================================================================================
*/

workflow {
    //  The script to execute is called by it's process name, and input is provided between brackets.

    NUM_LINES(input_ch)

    /*  Process output is accessed using the `out` channel.
        The channel operator view() is used to print process output to the terminal. */

    NUM_LINES.out.view()

}

/*
========================================================================================
    A Nextflow process block. Process names are written, by convention, in uppercase.
    This convention is used to enhance workflow readability.
========================================================================================
*/
process NUM_LINES {

    input:
    path read

    output:
    stdout

    script:
    """
    # Print reads
    printf '${read}\t'

    # Unzip file and count number of lines
    gunzip -c ${read} | wc -l
    """
}
~~~~
{: .language-groovy}

This is a Nextflow script, which contains the following:

1. An optional interpreter directive ("Shebang") line, specifying the location of the Nextflow interpreter.
1. `nextflow.enable.dsl=2` to enable DSL2 syntax.
1. A multi-line Nextflow comment, written using C style block comments, followed by a single line comment.
1. A pipeline parameter `params.input` which is given a default value, of the relative path to the location of a compressed fastq file, as a string.
1. An unnamed `workflow` execution block, which is the default workflow to run.
1. A Nextflow channel used to read in data to the workflow.
1. A call to the process `NUM_LINES`.
1. A Nextflow process block named `NUM_LINES`, which defines what the process does.
1. An `input` definition block that assigns the `input` to the variable `read`, and declares that it should be interpreted as a file path.
1. An `output` definition block that uses the Linux/Unix standard output stream `stdout` from the script block.
1. A script block that contains the bash commands `printf '${read}'` and `gunzip -c ${read} | wc -l`.
1. A Nextflow channel `input_ch` used to read in data to the workflow.
1. An unnamed `workflow` execution block, which is the default workflow to run.
1. A call to the process `NUM_LINES` with input channel `input_ch`.
1. An operation on the process output, using the channel operator `.view()`.

## Running Nextflow scripts

To run a Nextflow script use the command `nextflow run <script_name>`.

> ## Run a Nextflow  script
> Run the script by entering the following command in your terminal:
>
> ~~~
> $ nextflow run word_count.nf
> ~~~
> {: .language-bash}
> > ## Solution
> > You should see output similar to the text shown below:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `word_count.nf` [fervent_babbage] - revision: c54a707593
> > executor >  local (1)
> > [21/b259be] process > NUM_LINES (1) [100%] 1 of 1 ✔
> >
> >  ref1_1.fq.gz 58708
> > ~~~
> > {: .output}
> >
> > 1. The first line shows the Nextflow version number.
> > 1. The second line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
> > 1. The third line tells you the process has been executed locally (`executor >  local`).
> > 1. The next line shows the process id `21/b259be`, process name, number of cpus, percentage task completion, and how many instances of the process have been run.
> > 1. The final line is the output of the `.view()` operator.
> {: .solution}
{: .challenge}

{% include links.md %}
