---
title: Getting Started with Nextflow
teaching: 30
exercises: 10
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand what a workflow management system is.
- Understand the benefits of using a workflow management system.
- Explain the benefits of using Nextflow as part of your bioinformatics workflow.
- Explain the components of a Nextflow script.
- Run a Nextflow script.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What is a workflow and what are workflow management systems?
- Why should I use a workflow management system?
- What is Nextflow?
- What are the main features of Nextflow?
- What are the main components of a Nextflow script?
- How do I run a Nextflow script?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Workflows

Analysing data involves a sequence of tasks, including gathering, cleaning, and processing data. This sequence of tasks is called a workflow or a pipeline. These workflows typically require executing multiple software packages, sometimes running on different computing environments, such as a desktop or a compute cluster. Traditionally these workflows have been joined together in scripts using general purpose programming languages such as Bash or Python.

<br>
<center>
    <img src="https://sateeshperi.github.io/nextflow_varcal/nextflow/images/analysis_workflow.PNG" width="200" alt="Flowchart illustrating a simple bioinformatics RNA-Seq pipeline for transcript expression quantification. The process begins with inputs at the top: 'Fastq' files, a 'Reference sequence', and 'Grch38 Ensembl 91' annotations. The first step is 'quality control', using 'fastQC' software version 0.11.9. The second step is 'index creation' with 'Salmon' software version 1.3.0 using the '-i' option. The third step is 'quantification', again with 'Salmon' version 1.3.0, this time using the '-l A' option. There are two outputs from this process: 'Output 1' is a 'QC report' and 'Output 2' is 'Transcript expression' data.">
    <br>
    <em> An example of a simple bioinformatics RNA-Seq pipeline. </em>
</center>

<br>

However, as workflows become larger and more complex, the management of the programming logic and software becomes difficult.

## Workflow management systems

*Workflow Management Systems* (WfMS) such as [Snakemake](https://snakemake.readthedocs.io/en/stable/),
[Galaxy](https://usegalaxy.org/), and [Nextflow](nextflow.io) have been
developed specifically to manage computational data-analysis workflows in
fields such as bioinformatics, imaging, physics, and chemistry. These systems
contain multiple features that simplify the development, monitoring, execution
and sharing of pipelines, such as:

- Run time management
- Software management
- Portability \& Interoperability
- Reproducibility
- Re-entrancy

<br>
<center>
<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-021-01254-9/MediaObjects/41592_2021_1254_Fig1_HTML.png?as=webp" alt="A comparison of three bioinformatics pipeline diagrams. Panel A shows an 'Analysis workflow' for transcript expression quantification with three main steps: 1) quality control using fastQC v0.11.9, 2) index creation with Salmon v.1.3.0, and 3) quantification also using Salmon v.1.3.0. Inputs include Fastq files, a Reference sequence, and Grch38 Ensembl 91, leading to outputs of a QC report and transcript expression data. Panel B illustrates a 'Traditional pipeline' emphasizing platform-specific requirements and local execution with steps leading to two outputs. Panel C depicts a 'Workflow manager', highlighting platform-independent requirements, portability, local and cloud execution options, scalability, and containerized steps for automatic resource management, leading to an output and an execution report. The color-coding indicates input data (gray), output data (yellow), and software, versions, parameters (green and blue)" width="600">
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
<img src="fig/execution_abstraction.png" alt= "Infographic illustrating the components and supported platforms of a nextflow pipeline. The top section 'nextflow pipeline' is divided into three: writing code in any language (represented by R, Python, and Bash icons), orchestrating tasks with dataflow programming (represented by papers marked 'Data Flow' and 'Programming Model'), and defining software dependencies via containers (represented by Conda, Docker, and Singularity icons) and built-in version control with Git (represented by Git, GitHub, GitLab, and Bitbucket icons). Below, in the 'nextflow runtime' section, is 'Task orchestration and execution'. Arrows point downwards to the 'Supported Platforms' section, showcasing various platforms such as AWS, Google Cloud, Azure, Grid Engine, Slurm, HTCondor, Platform Computing, Kubernetes, and PBS Works." width="600">
<br>
<em> Overview of Nextflow core features. </em>
</center>
<br>

- **Fast prototyping**: A simple syntax for writing pipelines that enables you
  to reuse existing scripts and tools for fast prototyping.

- **Reproducibility**: Nextflow supports several container technologies, such
  as [Docker](https://www.docker.com/) and [Singularity](https://sylabs.io/singularity),
  as well as the package manager [Conda](https://docs.conda.io). This, along
  with the integration of the [GitHub](https://www.github.com) code sharing
  platform, allows you to write self-contained pipelines, manage versions and
  to reproduce any previous result when re-run, including on different
  computing platforms.

- **Portability \& interoperability**: Nextflow's syntax separates the
  functional logic (the steps of the workflow) from the execution settings (how
  the workflow is executed). This allows the pipeline to be run on multiple
  platforms, e.g. local compute vs. a university compute cluster or a cloud
  service like [AWS](https://aws.amazon.com/), without changing the steps of
  the workflow.

- **Simple parallelism**:  Nextflow is based on the dataflow programming model
  which greatly simplifies the splitting of tasks that can be run at the same
  time (parallelisation).

- **Continuous checkpoints \& re-entrancy**: All the intermediate results
  produced during the pipeline execution are automatically tracked. This allows
  you to resume its execution from the last successfully executed step, no
  matter what the reason was for it stopping.

## Processes, channels, and workflows

Nextflow workflows have three main parts: *processes*, *channels*, and
*workflows*.

- *Processes* describe a task to be run. A process script can be written in any
  scripting language that can be executed by the Linux platform (Bash, Perl,
  Ruby, Python, R, etc.). Processes spawn a task for each complete input set.
  Each task is executed independently and cannot interact with other tasks. The
  only way data can be passed between process tasks is via asynchronous queues,
  called *channels*.

- Processes define inputs and outputs for a task. *Channels* are then used to
  manipulate the flow of data from one process to the next.

- The interaction between processes, and ultimately the pipeline execution flow
  itself, is then explicitly defined in a *workflow* section.

In the following example we have a channel containing three elements, e.g.,
three data files. We have a process that takes the channel as input. Since the
channel has three elements, three independent instances (tasks) of that process
are run in parallel. Each task generates an output, which is passed to another
channel and used as input for the next process.

<p align="center">   <img src="fig/channel-process_fqc.png" alt="Diagram depicting part of a bioinformatics data processing workflow. On the left, there is a 'channel' labeled 'samples' containing three items: Fastq1, Fastq2, and Fastq3. This channel flows into a 'process' called 'fastqc' represented by a rounded rectangle containing the command 'fastqc -o out ${reads}'. The output of this process goes into a channel named 'out_ch', which lists 'outdir' three times as its contents. This channel then flows into a channel operator 'collect' and then into  another 'process' called 'multiqc', indicated by a rounded rectangle with the command 'multiqc -o mqc_res .'. The output of 'multiqc' goes into a channel called 'mqc_ch', which also lists 'outdir' one time." width="700">   <br>   <em> Nextflow process flow diagram. </em>
</p>

## Workflow execution

While a `process` defines what command or script has to be executed, the
`executor` determines how that script is actually run in the target system.

If not otherwise specified, processes are executed on the local computer. The
local executor is very useful for pipeline development, testing, and
small-scale workflows, but for large-scale computational pipelines, a High
Performance Cluster (HPC) or Cloud platform is often required.

<p align="center">   <img  src="fig/executor.png" alt="Diagram of a computational process within a bioinformatics workflow. The image features a large, central, rounded rectangle labeled 'process' with a smaller rectangle inside it labeled 'script', indicating the code or commands that are being executed. Above the script box, there is a smaller inset labeled 'Executors' with three icons: a desktop computer labeled 'Local', a stack of servers labeled 'High Performance Compute Cluster', and a cloud symbol labeled 'Cloud Compute'. These represent the different computing environments where the script can be executed. To the left of the process box is a green left-pointing arrowhead, suggesting input into the process, and to the right is a yellow right-pointing arrowhead, indicating the direction of output from the process." width="350">   <br>   <em>Nextflow Executors</em>
</p>

Nextflow provides a separation between the pipeline's functional logic and the
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

::::::::::::::::::::::::::::::::::::: instructor

The nextflow scripts for each lesson episode are available in the `scripts` directory created
during the course setup. You should copy the script into the current directory

For example, to copy the script for this lesson episode, run the following command:

```bash
$ cp scripts/introduction/word_count.nf .
```

:::::::::::::::::::::::::::::::::::::::::::::::::

```groovy
#!/usr/bin/env nextflow



/*
========================================================================================
    Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator.
========================================================================================
*/

params.input = "data/yeast/reads/ref1_1.fq.gz"

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
    # Print file name
    printf '${read}\\t'

    # Unzip file and count number of lines
    gunzip -c ${read} | wc -l
    """
}
```

This is a Nextflow script, which contains the following:

1. An optional interpreter directive ("Shebang") line, specifying the location of the Nextflow interpreter.
2. A multi-line Nextflow comment, written using C style block comments, there are more comments later in the file.
3. A pipeline parameter `params.input` which is given a default value, of the relative path to the location of a compressed fastq file, as a string.
4. A Nextflow channel `input_ch` used to read in data to the workflow.
5. An unnamed `workflow` execution block, which is the default workflow to run.
6. A call to the process `NUM_LINES`.
7. An operation on the process output, using the channel operator `.view()`.
8. A Nextflow process block named `NUM_LINES`, which defines what the process does.
9. An `input` definition block that assigns the `input` to the variable `read`, and declares that it should be interpreted as a file path.
10. An `output` definition block that uses the Linux/Unix standard output stream `stdout` from the script block.
11. A script block that contains the bash commands `printf '${read}'` and `gunzip -c ${read} | wc -l`.

## Running Nextflow scripts

To run a Nextflow script use the command `nextflow run <script_name>`.

:::::::::::::::::::::::::::::::::::::::  challenge

## Run a Nextflow  script

Run the script by entering the following command in your terminal:

```bash
$ nextflow run word_count.nf
```

:::::::::::::::  solution

## Solution

You should see output similar to the text shown below:

```output
N E X T F L O W  ~  version 21.04.3
Launching `word_count.nf` [fervent_babbage] - revision: c54a707593
executor >  local (1)
[21/b259be] process > NUM_LINES (1) [100%] 1 of 1 ✔

 ref1_1.fq.gz 58708
```

1. The first line shows the Nextflow version number.
2. The second line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
3. The third line tells you the process has been executed locally (`executor >  local`).
4. The next line shows the process id `21/b259be`, process name, number of cpus, percentage task completion, and how many instances of the process have been run.
5. The final line is the output of the `.view()` operator.
  
  

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::



## Quick recap

- A workflow is a sequence of tasks that process a set of data, and a workflow
  management system (WfMS) is a computational platform that provides an
  infrastructure for the set-up, execution and monitoring of workflows.
- Nextflow scripts comprise of *channels* for controlling inputs and outputs,
  and *processes* for defining workflow tasks.
- You run a Nextflow script using the `nextflow run` command.



:::::::::::::::::::::::::::::::::::::::: keypoints

- A workflow is a sequence of tasks that process a set of data.
- A workflow management system (WfMS) is a computational platform that provides an infrastructure for the set-up, execution and monitoring of workflows.
- Nextflow is a workflow management system that comprises both a runtime environment and a domain specific language (DSL).
- Nextflow scripts comprise of channels for controlling inputs and outputs, and processes for defining workflow tasks.
- You run a Nextflow script using the `nextflow run` command.

::::::::::::::::::::::::::::::::::::::::::::::::::


