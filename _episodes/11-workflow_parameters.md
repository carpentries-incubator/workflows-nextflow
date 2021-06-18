---
title: "Workflow parameterization"
teaching: 30
exercises: 10
questions:
- "How can I change the data a workflow uses?"
- "How can parameterize a workflow?"
objectives:
- "FIXME"
keypoints:
- "FIXME"
---

## Nextflow log

Once a script has run, Nextflow stores a log of all the workflows executed in the current folder.
Similar to an electronic lab book, this means you have a have a record of all processing steps and commands run.

You can print Nextflow's execution history and log information using the  `nextflow log` command.

> ## Show Execution Log
> Listing the execution logs of previous invocations of all pipelines in a directory.
>
> ~~~
> $ nextflow log
> ~~~
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
> >TIMESTAMP          	DURATION	RUN NAME       	STATUS	REVISION ID	SESSION ID                          	COMMAND
> >2021-03-19 13:45:53	6.5s    	fervent_babbage	OK    	c54a707593 	15487395-443a-4835-9198-229f6ad7a7fd	nextflow run wc.nf
> > 2021-03-19 13:46:53	6.6s    	soggy_miescher 	OK    	c54a707593 	58da0ccf-63f9-42e4-ba4b-1c349348ece5	nextflow run wc.nf --samples 'data/yeast/reads/*.fq.gz'
> >  ~~~
> > {: .output }
> {: .solution}
{: .challenge}


## Modify and resume

When Nextflow is run, it runs the entire workflow by default.
However, Nextflow keeps track of all the processes executed in your pipeline. By using the Nextflow specific parameter `-resume`, Nextflow
will start from the last successfully executed process. If you modify some parts of your script, only the processes that are actually changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.


> ## Re-run the pipeline using -resume option
> Execute the script by entering the following command in your terminal:
>
> ~~~
> $ nextflow run wc.nf --samples 'data/ggal/*.fq' -resume
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
> > {: .output }
> {: .solution}
{: .challenge}




You will see that the execution of the process `numLines` is actually skipped (cached text appears), and its results are retrieved from the cache.


> ## Modify the wc.nf script and re-run the pipeline using -resume option
> Modify the wc.nf script changing the sleep time and execute the script by entering the following command in your terminal:
>
> ~~~
> $ nextflow run wc.nf --samples 'data/yeast/reads/*.fq.gz' -resume
> ~~~
> {: .language-bash}
> > ## Solution
> > The output will look similar to this:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `wc.nf` [backstabbing_joliot] - revision: 714e17a273
[fc/deef86] process > numLines (18) [100%] 18 of 18, cached: 18 ✔
> > 3252 ref3_1.fq.gz
> >
> >5393 temp33_3_1.fq.gz
> >
> >3602 ref1_2.fq.gz
> >
> >4950 temp33_1_2.fq.gz
> >
> >5040 ref2_1.fq.gz
> >
> >5347 etoh60_1_1.fq.gz
> >
> >3823 temp33_2_1.fq.gz
> >
> >6320 etoh60_3_2.fq.gz
> >
> >6327 etoh60_3_1.fq.gz
> >
> >4904 temp33_1_1.fq.gz
> >
> >5038 ref2_2.fq.gz
> >
> >3628 ref1_1.fq.gz
> >
> >5434 etoh60_2_2.fq.gz
> >
> >3293 ref3_2.fq.gz
> >
> >3858 temp33_2_2.fq.gz
> >
> >5386 temp33_3_2.fq.gz
> >
> >5371 etoh60_1_2.fq.gz
> >
> >5468 etoh60_2_1.fq.gz
> >  ~~~
> > {: .output }
> {: .solution}
{: .challenge}

As you have changed the script the pipeline will re-run and won't use the cached results for that process.


~~~
$ nextflow log
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

> ## Nextflow specific and workflow specific parameters
> Command line parameters that start with a single dash e.g., `-resume`,
> are parameters specifically for Nextflow to interpret.
>
> Command line parameters that start with a double dash e.g., `--samples`,
> are parameters to your workflow script and can be accessed via
> the `params.<variable>` variable.
{: .callout}


## Work directory

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

### Specifying another work directory

Depending on your script, this work folder can take a lot of disk space.
You can specify another work directory using the command line option `-w`

~~~
$ nextflow run <script> -w /some/scratch/dir
~~~
{: .language-bash}

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
