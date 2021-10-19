---
title: "Workflow caching and checkpointing"
teaching: 20
exercises: 10
questions:
- "How can I restart a Nextflow workflow after an error?"
- "How can I add new data to a workflow?"
- "Where can I find intermediate data and results?"
objectives:
- "Resume a Nextflow workflow using the `-resume` option."
- "Restart a Nextflow workflow using new data."
keypoints:
- "Nextflow automatically keeping track of all the processes executed in your pipeline via  caching  and checkpointing."
- "You can restart a Nextflow workflow using new data using new data skipping steps that have been successfully executed."
- "Nextflow stores intermediate data in a working directory."
---

A key features of workflow management systems, like Nextflow, is the ability to restart a pipeline after an error from the last successful process. Nextflow achieves this by automatically keeping track of all the processes executed in your pipeline via  caching  and checkpointing.

## Resume

To restart from the last successfully executed process we add the command line option `-resume` to the Nextflow command.

For example, the command below would resume the `wc.nf` script from the last successful process.

~~~
$ nextflow run wc.nf --input 'data/yeast/reads/ref1*.fq.gz' -resume
~~~
{: .language-bash}

We can see in the output that the results from the process `NUM_LINES` has been retrieved from the cache.
~~~
Launching `wc.nf` [condescending_dalembert] - revision: fede04a544
[c9/2597d5] process > NUM_LINES (1) [100%] 2 of 2, cached: 2 ✔
ref1_1.fq.gz 58708

ref1_2.fq.gz 58708
~~~
{: .output}


> ## Resume a pipeline
> Resume the Nextflow script `wc.nf` by re-running the command and adding the parameter `-resume`
> and the parameter `--input 'data/yeast/reads/temp33*'`:
>
> {: .language-bash}
> > ## Solution
> > ~~~
> > $ nextflow run wc.nf --input 'data/yeast/reads/temp33*' -resume
> > ~~~
> > If your previous run was successful the output will look similar to this:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [nauseous_leavitt] - revision: fede04a544
> > [21/6116de] process > NUM_LINES (4) [100%] 6 of 6, cached: 6 ✔
> >temp33_3_2.fq.gz 88956
> >
> >temp33_3_1.fq.gz 88956
> >
> >temp33_1_1.fq.gz 82372
> >
> >temp33_2_2.fq.gz 63116
> >
> >temp33_1_2.fq.gz 82372
> >
> >temp33_2_1.fq.gz 63116
> >  ~~~
> > {: .output }
> > You will see that the execution of the process `NUMLINES` is actually skipped (cached text appears), and its results are retrieved from the cache.
> {: .solution}
{: .challenge}

## How does resume work?

The mechanism works by assigning a unique ID to each task. This unique ID is used to create a separate execution directory, within the `work` directory, where the tasks are executed and the results stored. A task’s unique ID is generated as a 128-bit hash number obtained from a composition of the task’s:

* Inputs values
* Input files
* Command line string
* Container ID
* Conda environment
* Environment modules
* Any executed scripts in the bin directory

When we resume a workflow Nextflow uses this unique ID to check if:

1. The working directory exists
1. It contains a valid command exit status
1. It contains the expected output files.

If these conditions are satisfied, the task execution is skipped and the previously computed outputs are applied. When a task requires recomputation, ie. the conditions above are not fulfilled, the downstream tasks are automatically invalidated.


Therefore, if you modify some parts of your script, or alter the input data using `-resume`, will only execute the processes that are actually changed.

The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

> ## Modify Nextflow script and re-run.
>  Alter the timestamp on the file temp33_3_2.fq.gz using the UNIX `touch` command.
> ~~~~
> $ touch data/yeast/reads/temp33_3_2.fq.gz
> ~~~~
> {: .language-bash}
> Run command below.
> ~~~
> $ nextflow run wc.nf --input 'data/yeast/reads/temp33*' -resume
> ~~~
> How many processes will be cached and how many will run ?
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [gigantic_minsky] - revision: fede04a544
> > executor >  local (1)
> > [20/cda0d5] process > NUM_LINES (5) [100%] 6 of 6, cached: 5 ✔
> > temp33_1_2.fq.gz 82372
> >
> > temp33_3_1.fq.gz 88956
> >
> > temp33_2_1.fq.gz 63116
> >
> > temp33_1_1.fq.gz 82372
> >
> > temp33_2_2.fq.gz 63116
> >
> > temp33_3_2.fq.gz 88956
> >  ~~~
> > {: .output }
> > As you changed the timestamp on one file it will only re-run that process.
> > The results from the other 5 processes are cached.
> {: .solution}
{: .challenge}


## The Work directory

By default the pipeline results are cached in the directory `work` where the pipeline is launched.

We can use the Bash `tree` command to list the contents of the work directory.
**Note:** By default tree does not print hidden files (those beginning with a dot `.`). Use the `-a`   to view  all files.

~~~
$ tree work
~~~
{: .language-bash}



~~~
work/
├── 02
│   └── 26541c31e597c2243b2489b06f51ef
│       └── ref1_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/ref1_2.fq.gz
├── 07
│   └── 54363a267def098c0544708d3f6dd3
│       └── temp33_1_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/temp33_1_2.fq.gz
├── 4a
│   └── aeed908acc5481ee887736386ac8b8
│       └── ref3_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/ref3_2.fq.gz
├── 57
│   ├── 6a47b10d1569561ae321720d0b8f15
│   │   └── etoh60_1_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/etoh60_1_2.fq.gz
│   └── 99ee9a312eb3f3b8afaa6cf881f972
│       └── etoh60_3_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/etoh60_3_1.fq.gz
├── 5e
│   └── 7ee2308079f20c96a7d5b70291e167
│       └── temp33_1_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/temp33_1_1.fq.gz
├── 7b
│   └── 963dc9c1592ade5978e5baa7619cef
│       └── etoh60_2_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/etoh60_2_1.fq.gz
├── 7e
│   ├── a706851ba2ad4f3ce796d6b55faadf
│   │   └── temp33_3_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/temp33_3_2.fq.gz
│   └── b00105a9510bc0182e0382e2a91710
│       └── temp33_2_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/temp33_2_2.fq.gz
├── 8d
│   └── 19a8e8da614510599fb4a8c1080176
│       └── temp33_2_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/temp33_2_1.fq.gz
├── 9f
│   └── edb44e296ae4710ade2da9c6a2cd13
│       └── temp33_3_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/temp33_3_1.fq.gz
├── b4
│   └── 81ccb20b2920b042c0916d7fb2c071
│       └── etoh60_3_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/etoh60_3_2.fq.gz
├── b5
│   └── 337cfa936e4670c386efcd93497dcd
│       └── etoh60_1_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/etoh60_1_1.fq.gz
├── c3
│   └── db85075e09b1743f73f4dc78657423
│       └── etoh60_2_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/etoh60_2_2.fq.gz
├── d7
│   └── c070aebb4232b6171043f1c29066e6
│       └── ref1_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/ref1_1.fq.gz
├── e3
│   └── 9c644d276e02018189862afedc6ab6
│       └── ref3_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/ref3_1.fq.gz
├── f0
│   └── 9421f2eb5b12ac3213e04e5682324e
│       └── ref2_1.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/ref2_1.fq.gz
└── fc
    └── 4022a522f2964dde97fa1484bd5742
        └── ref2_2.fq.gz -> /Users/ggrimes2/Downloads/nextflow_rnaseq_training_dataset/data/yeast/reads/ref2_2.fq.gz
~~~
{: .output }

You will see the input Fastq files are symbolically linked to their original location.

### Task execution directory

Within the `work` directory there are multiple task execution directories. 
There is one directory for each time a  process is executed. 
These task directories are identified by the process execution hash e.g. `fc/4022a5` .

The task execution directory contains:

* `.command.sh`: The command script.

* `.command.run`: The command wrapped used to run the job.

* `.command.out`: The complete job standard output.

* `.command.err`: The complete job standard error.

* `.command.log`: The wrapper execution output.

* `.command.begin`: Sentinel file created as soon as the job is launched.

* `.exitcode`: A file containing the task exit code.

* Task input files (symlinks)

* Task output files



### Specifying another work directory

Depending on your script, this work folder can take a lot of disk space.
You can specify another work directory using the command line option `-w`

~~~
$ nextflow run <script> -w /some/scratch/dir
~~~
{: .language-bash}

### Clean the work directory

If you are sure you won’t resume your pipeline execution, clean this folder periodically using the command `nextflow clean`.

~~~
$ nextflow clean [run_name|session_id] [options]
~~~
{: .language-bash}

Typically, results before the last successful result are cleaned:

~~~
$ nextflow clean -f -before [run_name|session_id]
~~~
{: .language-bash}

> ## Remove a Nextflow run.
>  Remove  the last nextflow run using the command `nextflow clean`. 
>  First use the option `-dry-run` to see which files would be deleted and then re-run removing the run and associated files.
> > ## Solution
> > An example nextflow clean command with `dry-run` .
> >
> > ~~~
> > $ nextflow clean nauseous_leavitt -dry-run
> > ~~~
> > {: .language-bash}
> > An example nextflow clean command removing the files.
> > ~~~
> > $ nextflow clean nauseous_leavitt -f
> > ~~~
> > {: .language-bash}
> {: .solution}
{: .challenge}


