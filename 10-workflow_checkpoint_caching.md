---
title: Workflow caching and checkpointing
teaching: 20
exercises: 10
---

::::::::::::::::::::::::::::::::::::::: objectives

- Resume a Nextflow workflow using the `-resume` option.
- Restart a Nextflow workflow using new data.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How can I restart a Nextflow workflow after an error?
- How can I add new data to a workflow without starting from the beginning?
- Where can I find intermediate data and results?

::::::::::::::::::::::::::::::::::::::::::::::::::

A key feature of workflow management systems, like Nextflow, is re-entrancy, which is the ability to restart a pipeline after an error from the last successfully executed process. Re-entrancy enables time consuming successfully completed steps, such as index creation, to be skipped when adding more data to a pipeline. This in turn leads to faster prototyping and development of workflows, and faster analyses of additional data.
Nextflow achieves re-entrancy by automatically keeping track of all the processes executed in your pipeline via  caching  and checkpointing.

## Resume

To restart from the last successfully executed process we add the command line option `-resume` to the Nextflow command.

For example, the command below would resume the `wc.nf` script from the last successful process.

```bash
$ nextflow run wc.nf --input 'data/yeast/reads/ref1*.fq.gz' -resume
```

We can see in the output that the results from the process `NUM_LINES` has been retrieved from the cache.

```output
Launching `wc.nf` [condescending_dalembert] - revision: fede04a544
[c9/2597d5] process > NUM_LINES (1) [100%] 2 of 2, cached: 2 ✔
ref1_1.fq.gz 58708

ref1_2.fq.gz 58708
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Resume a pipeline

Resume the Nextflow script `wc.nf` by re-running the command and adding the parameter `-resume`
and the parameter `--input 'data/yeast/reads/temp33*'`:

:::::::::::::::  solution

## Solution

```
$ nextflow run wc.nf --input 'data/yeast/reads/temp33*' -resume
```

If your previous run was successful the output will look similar to this:

```output 
N E X T F L O W  ~  version 20.10.0
Launching `wc.nf` [nauseous_leavitt] - revision: fede04a544
[21/6116de] process > NUM_LINES (4) [100%] 6 of 6, cached: 6 ✔
temp33_3_2.fq.gz 88956

temp33_3_1.fq.gz 88956

temp33_1_1.fq.gz 82372

temp33_2_2.fq.gz 63116

temp33_1_2.fq.gz 82372

temp33_2_1.fq.gz 63116
```

You will see that the execution of the process `NUMLINES` is actually skipped (cached text appears), and the results are retrieved from the cache.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## How does resume work?

The mechanism works by assigning a unique ID to each task. This unique ID is used to create a separate execution directory, within the `work` directory, where the tasks are executed and the results stored. A task's unique ID is generated as a 128-bit hash number obtained from a composition of the task's:

- Inputs values
- Input files
- Command line string
- Container ID
- Conda environment
- Environment modules
- Any executed scripts in the bin directory

When we resume a workflow Nextflow uses this unique ID to check if:

1. The working directory exists
2. It contains a valid command exit status
3. It contains the expected output files.

If these conditions are satisfied, the task execution is skipped and the previously computed outputs are applied. When a task requires recomputation, ie. the conditions above are not fulfilled, the downstream tasks are automatically invalidated.

Therefore, if you modify some parts of your script, or alter the input data using `-resume`, will only execute the processes that are actually changed.

The execution of the processes that are not changed will be skipped and the cached result used instead.

This helps a lot when testing or modifying part of your pipeline without having to re-execute it from scratch.

:::::::::::::::::::::::::::::::::::::::  challenge

## Modify Nextflow script and re-run.

Alter the timestamp on the file temp33\_3\_2.fq.gz using the UNIX `touch` command.

```bash
$ touch data/yeast/reads/temp33_3_2.fq.gz
```

Run command below.

```
$ nextflow run wc.nf --input 'data/yeast/reads/temp33*' -resume
```

How many processes will be cached and how many will run ?
{: .language-bash}

:::::::::::::::  solution

## Solution

The output will look similar to this:

```output 
N E X T F L O W  ~  version 20.10.0
Launching `wc.nf` [gigantic_minsky] - revision: fede04a544
executor >  local (1)
[20/cda0d5] process > NUM_LINES (5) [100%] 6 of 6, cached: 5 ✔
temp33_1_2.fq.gz 82372

temp33_3_1.fq.gz 88956

temp33_2_1.fq.gz 63116

temp33_1_1.fq.gz 82372

temp33_2_2.fq.gz 63116

temp33_3_2.fq.gz 88956
```

As you changed the timestamp on one file it will only re-run that process.
The results for the other 5 processes are retrieved from the cache.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## The Work directory

By default the pipeline results are cached in the directory `work` where the pipeline is launched.

We can use the Bash `tree` command to list the contents of the work directory.
**Note:** By default tree does not print hidden files (those beginning with a dot `.`). Use the `-a`   to view  all files.

```bash
$ tree -a work
```

Example output

```output 
work/
├── 12
│   └── 5489f3c7dbd521c0e43f43b4c1f352
│       ├── .command.begin
│       ├── .command.err
│       ├── .command.log
│       ├── .command.out
│       ├── .command.run
│       ├── .command.sh
│       ├── .exitcode
│       └── temp33_1_2.fq.gz -> /home/training/data/yeast/reads/temp33_1_2.fq.gz
├── 3b
│   └── a3fb24ad3242e4cc8e5aa0c24d174b
│       ├── .command.begin
│       ├── .command.err
│       ├── .command.log
│       ├── .command.out
│       ├── .command.run
│       ├── .command.sh
│       ├── .exitcode
│       └── temp33_2_1.fq.gz -> /home/training/data/yeast/reads/temp33_2_1.fq.gz
├── 4c
│   └── 125b5e5a5ee144fa25dd9bccd467e9
│       ├── .command.begin
│       ├── .command.err
│       ├── .command.log
│       ├── .command.out
│       ├── .command.run
│       ├── .command.sh
│       ├── .exitcode
│       └── temp33_3_1.fq.gz -> /home/training/data/yeast/reads/temp33_3_1.fq.gz
├── 54
│   └── eb9d72e9ac24af8183de569ab0b977
│       ├── .command.begin
│       ├── .command.err
│       ├── .command.log
│       ├── .command.out
│       ├── .command.run
│       ├── .command.sh
│       ├── .exitcode
│       └── temp33_2_2.fq.gz -> /home/training/data/yeast/reads/temp33_2_2.fq.gz
├── e9
│   └── 31f28c291481342cc45d4e176a200a
│       ├── .command.begin
│       ├── .command.err
│       ├── .command.log
│       ├── .command.out
│       ├── .command.run
│       ├── .command.sh
│       ├── .exitcode
│       └── temp33_1_1.fq.gz -> /home/training/data/yeast/reads/temp33_1_1.fq.gz
└── fa
    └── cd3e49b63eadd6248aa357083763c1
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── temp33_3_2.fq.gz -> /home/training/data/yeast/reads/temp33_3_2.fq.gz
```

### Task execution directory

Within the `work` directory there are multiple task execution directories.
There is one directory for each time a process is executed.
These task directories are identified by the process execution hash. For example
the task directory `fa/cd3e49b63eadd6248aa357083763c1` would be location for the process identified by the hash `fa/cd3e49` .

The task execution directory contains:

- `.command.sh`: The command script.

- `.command.run`: The command wrapped used to run the job.

- `.command.out`: The complete job standard output.

- `.command.err`: The complete job standard error.

- `.command.log`: The wrapper execution output.

- `.command.begin`: A file created as soon as the job is launched.

- `.exitcode`: A file containing the task exit code.

- Any task input files (symlinks)

- Any task output files

### Specifying another work directory

Depending on your script, this work folder can take a lot of disk space.
You can specify another work directory using the command line option `-w`.
**Note** Using a different work directory will mean that any jobs will need to re-run from the beginning.

```bash
$ nextflow run wc.nf --input 'data/yeast/reads/temp33*' -w second_work_dir -resume
```

```output
N E X T F L O W  ~  version 21.04.0
Launching `wc.nf` [deadly_easley] - revision: fede04a544
executor >  local (6)
[9d/0f5e89] process > NUM_LINES (5) [100%] 6 of 6 ✔
temp33_3_2.fq.gz 88956

temp33_1_1.fq.gz 82372

temp33_3_1.fq.gz 88956

temp33_1_2.fq.gz 82372

temp33_2_2.fq.gz 63116

temp33_2_1.fq.gz 63116
```

### Clean the work directory

If you are sure you won't resume your pipeline execution, clean this folder periodically using the command `nextflow clean`.

```bash
$ nextflow clean [run_name|session_id] [options]
```

Supply the option `-n` to print names of files to be removed without deleting them, or `-f` to force the removal of the files. If you only want to remove files from a run but retain execution log entries and metadata, add the option `-k`. Multiple runs can be cleaned with the options, `-before`, `-after` or `-but` before the run name.
For example, the command below would remove all the temporary files and log entries for runs before the run `gigantic_minsky`.

```bash
$ nextflow clean -f -before gigantic_minsky
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Remove a Nextflow run.

Remove the last Nextflow run using the command `nextflow clean`.
First use the option `-dry-run` to see which files would be deleted and then re-run removing the run and associated files.

:::::::::::::::  solution

## Solution

An example nextflow clean command with `dry-run` .

```bash
$ nextflow clean nauseous_leavitt -dry-run
```

An example nextflow clean command removing the files.

```bash
$ nextflow clean nauseous_leavitt -f
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: keypoints

- Nextflow automatically keeps track of all the processes executed in your pipeline via checkpointing.
- Nextflow caches intermediate data in task directories within the work directory.
- Nextflow caching and checkpointing allows re-entrancy into a workflow after a pipeline error or using new data, skipping steps that have been successfully executed. - Re-entrancy is enabled using the `-resume` option.

::::::::::::::::::::::::::::::::::::::::::::::::::


