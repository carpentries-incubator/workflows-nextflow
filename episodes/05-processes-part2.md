---
title: Processes Part 2
teaching: 30
exercises: 10
---

::::::::::::::::::::::::::::::::::::::: objectives

- Define outputs to a process.
- Understand how to handle grouped input and output using the tuple qualifier.
- Understand how to use conditionals to control process execution.
- Use process directives to control execution of a process.
- Use the `publishDir` directive to save result files to a directory.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I get data, files, and values,  out of processes?
- How do I handle grouped input and output?
- How can I control when a process is executed?
- How do I control resources, such as number of CPUs and memory, available to processes?
- How do I save output/results from a process?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Outputs

We have seen how to input data into a process; now we will see how to output files and values from a process.

The `output` declaration block allows us to define the channels used by the process to send out the files and values produced.

An output block is not required, but if it is present it can contain one or more output declarations.

The output block follows the syntax shown below:

```groovy 
output:
  <output qualifier> <output name>
  <output qualifier> <output name>
  ...
```

### Output values

Like the input, the type of output data is defined using type qualifiers.

The `val` qualifier allows us to output a value defined in the script.

Because Nextflow processes can only communicate through channels, if we want to share a value output of one process as input to another process, we would need to define that value in the output declaration block as shown in the following example:

```groovy 
//process_output_value.nf


params.transcriptome="${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

process COUNT_CHR_SEQS {
  input:
  val chr

  output:
  val chr

  script:
  """
  zgrep -c '^>Y'$chr $params.transcriptome
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
chr_ch = channel.of('A'..'P')

workflow {
  COUNT_CHR_SEQS(chr_ch)
  // use the view operator to display contents of the channel
  COUNT_CHR_SEQS.out.view()
}
```

```output 
N E X T F L O W  ~  version 21.10.6
Launching `p1.nf` [jovial_lavoisier] - revision: a652ef75d4
executor >  local (16)
executor >  local (16)
[6a/d82669] process > COUNT_CHR_SEQS (16) [100%] 16 of 16 ✔
B
456

A
118

C
186

[..truncated..]

```

### Output files

If we want to capture a file instead of a value as output we can use the
`path` qualifier that can capture one or more files produced by the process, over the specified channel.

```groovy 
//process_output_file.nf


params.transcriptome="${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

process COUNT_CHR_SEQS {
  input:
  val chr

  output:
  path "${chr}_seq_count.txt"

  script:
  """
  zgrep -c '^>Y'$chr $params.transcriptome > ${chr}_seq_count.txt
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
chr_ch = channel.of('A'..'P')

workflow {
  COUNT_CHR_SEQS(chr_ch)
  // use the view operator to display contents of the channel
  COUNT_CHR_SEQS.out.view()
}
```

```output 
N E X T F L O W  ~  version 21.10.6
Launching `process_output_file.nf` [angry_lichterman] - revision: 6a46c69413
executor >  local (16)
[95/ec5d62] process > COUNT_CHR_SEQS (13) [100%] 16 of 16 ✔
/Users/ggrimes2/Documents/process_wf/work/f2/6d5c44985a15feb0555b7b71c37a9c/J_seq_count.txt
executor >  local (16)
[95/ec5d62] process > COUNT_CHR_SEQS (13) [100%] 16 of 16 ✔
work/f2/6d5c44985a15feb0555b7b71c37a9c/J_seq_count.txt
work/4f/f810942341d003acc80c2603671177/B_seq_count.txt
work/23/883ccf187b5357137a9a87d98717c0/I_seq_count.txt
[..truncated..]
```

In the above example the process `COUNT_CHR_SEQS` creates a file named `<chr>_seq_count.txt` in the work directory containing the number of transcripts within that chromosome.

Since a file parameter using the same name, `<chr>_seq_count.txt`, is declared in the output block, when the task is completed that file is sent over the output channel.

A downstream `operator`, such as `.view` or a `process` declaring the same channel as input will be able to receive it.

### Multiple output files

When an output file name contains a `*` or `?` metacharacter it is interpreted as a pattern match.
This allows us to capture multiple files into a list and output them as a one item channel.

For example, here we will capture the files `sequence_ids.txt` and  `sequence.txt` produced as results from SPLIT\_FASTA in the output channel.

```groovy 
//process_output_multiple.nf


params.transcriptome="${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

process SPLIT_FASTA {
  input:
  path transcriptome

  output:
  path "*"

  script:
  """
  zgrep  '^>' $transcriptome > sequence_ids.txt
  zgrep -v '^>' $transcriptome > sequence.txt
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
transcriptome_ch = channel.fromPath(params.transcriptome)

workflow {
  SPLIT_FASTA(transcriptome_ch)
  // use the view operator to display contents of the channel
  SPLIT_FASTA.out.view()
}
```

```bash 
$ nextflow run process_output_multiple.nf
```

```output
N E X T F L O W  ~  version 21.10.6
Launching `process_output_multiple.nf` [goofy_meitner] - revision: 53cbf7e5a4
executor >  local (1)
[21/01e6ba] process > SPLIT_FASTA (1) [100%] 1 of 1 ✔
[/work/21/01e6baac41d2f37531f86dc7a57034/sequence.txt, work/21/01e6baac41d2f37531f86dc7a57034/sequence_ids.txt]

```

**Note:** There are some caveats on glob pattern behaviour:

- Input files are not included in the list of possible matches.
- Glob pattern matches against both files and directories path.
- When a two stars pattern `**` is used to recurse through subdirectories, only file paths are matched i.e. directories are not included in the result list.

:::::::::::::::::::::::::::::::::::::::  challenge

## Output channels

Modify the nextflow script `process_exercise_output.nf` to include an output block that captures the different output file `${chr}_seqids.txt`.

```groovy 
//process_exercise_output.nf


process EXTRACT_IDS {
  input:
  path transcriptome
  each chr

  //add output block here to capture the file "${chr}_seqids.txt"

  script:
  """
  zgrep '^>Y'$chr $transcriptome > ${chr}_seqids.txt
  """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
chr_ch = channel.of('A'..'P')

workflow {
  EXTRACT_IDS(transcriptome_ch, chr_ch)
  EXTRACT_IDS.out.view()
}
```

:::::::::::::::  solution

## Solution

```groovy 
//process_exercise_output_answer.nf


process EXTRACT_IDS {
  input:
  path transcriptome
  each chr

  //add output block here to capture the file "${chr}_seqids.txt"
  output:
  path "${chr}_seqids.txt"

  script:
  """
  zgrep '^>Y'$chr $transcriptome > ${chr}_seqids.txt
  """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
chr_ch = channel.of('A'..'P')

workflow {
  EXTRACT_IDS(transcriptome_ch, chr_ch)
  EXTRACT_IDS.out.view()
}
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Grouped inputs and outputs

So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the `tuple` qualifiers.

In tuples the first item is the grouping key and the second item is the list.

```
[group_key,[file1,file2,...]]
```

When using channel containing a tuple, such a one created with `.filesFromPairs` factory method, the corresponding input declaration must be declared with a `tuple` qualifier, followed by definition of each item in the tuple.

```groovy 
//process_tuple_input.nf


process TUPLEINPUT{
  input:
  tuple val(sample_id), path(reads)

  script:
  """
  echo $sample_id
  echo $reads
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  TUPLEINPUT(reads_ch)
}

```

outputs

```output 
ref1
ref1_1.fq.gz ref1_2.fq.gz
```

In the same manner an output channel containing tuple of values can be declared using the `tuple` qualifier following by the definition of each tuple element in the tuple.

In the code snippet below the output channel would contain a tuple with the grouping key value as the Nextflow variable `sample_id` and a list containing the files matching the following pattern `"${sample_id}.fq.gz"`.

```groovy 
output:
  tuple val(sample_id), path("${sample_id}.fq.gz")
```

An example can be seen in this script below.

```groovy 
//process_tuple_io.nf


process COMBINE_FQ {
  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.fq.gz")

  script:
  """
  cat $reads > ${sample_id}.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  COMBINE_FQ(reads_ch)
  COMBINE_FQ.out.view()
}

```

```bash 
nextflow run process_tuple_io.nf
```

The output is now a tuple containing the sample id and the combined fastq files.

```output 
[ref1, work/2d/a073d34b5b3231b1f57872599bd308/ref1.fq]
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Composite inputs and outputs

Fill in the blank \_\_\_ input and output qualifiers for `process_exercise_tuple.nf`.
**Note:** the output for the COMBINE\_REPS process.

```groovy 
//process_exercise_tuple.nf


process COMBINE_REPS {
  input:
  tuple ___(sample_id), ___(reads)

  output:
  tuple ___(sample_id), ___("*.fq.gz")

  script:
  """
  cat *_1.fq.gz > ${sample_id}_R1.fq.gz
  cat *_2.fq.gz > ${sample_id}_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref{1,2,3}*.fq.gz',size:-1)

workflow{
  COMBINE_REPS(reads_ch)
  COMBINE_REPS.out.view()
}
```

:::::::::::::::  solution

## Solution

```groovy 
//process_exercise_tuple_answer.nf


process COMBINE_REPS {
  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*.fq.gz")

  script:
  """
  cat *_1.fq.gz > ${sample_id}_R1.fq.gz
  cat *_2.fq.gz > ${sample_id}_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2,3}.fq.gz',size:-1)

workflow{
  COMBINE_REPS(reads_ch)
  COMBINE_REPS.out.view()
}
```

```output
N E X T F L O W  ~  version 21.04.0
Launching `process_exercise_tuple.nf` [spontaneous_coulomb] - revision: 06ff22f1a9
executor >  local (3)
[75/f4a44d] process > COMBINE_REPS (3) [100%] 3 of 3 ✔
[ref3, work/99/a7d9176e332fdc0988973dbb89df63/ref3_R1.fq.gz, work/99/a7d9176e332fdc0988973dbb89df63/ref3_R2.fq.gz]
[ref2, /work/53/e3cbd39afa9f0f84a3d9cd060d991a/ref2_R1.fq.gz, /work/53/e3cbd39afa9f0f84a3d9cd060d991a/ref2_R2.fq.gz]
[ref1, work/75/f4a44d0bc761fa4774c2f23a465766/ref1_R1.fq.gz, work/75/f4a44d0bc761fa4774c2f23a465766/ref1_R2.fq.gz]
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Conditional execution of a process

The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value; `true` or `false`.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters.

In the example below the process `CONDITIONAL` will only execute when the value of the `chr` variable is less than or equal to 5:

```groovy 
//process_when.nf


process CONDITIONAL {
  input:
  val chr

  when:
  chr <= 5

  script:
  """
  echo $chr
  """
}

chr_ch = channel.of(1..22)

workflow {
  CONDITIONAL(chr_ch)
}
```

```output 
4

5

2

3

1
```

## Directives

Directive declarations allow the definition of optional settings, like the number of `cpus` and amount of `memory`, that affect the execution of the current process without affecting the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e. `input`, `output`, etc).

**Note:** You do not use `=` when assigning a value to a directive.

Directives are commonly used to define the amount of computing resources to be used or extra information for configuration or logging purpose.

For example:

```groovy 
//process_directive.nf


process PRINTCHR {
  tag "tagging with chr$chr"
  cpus 1
  echo true

  input:
  val chr

  script:
  """
  echo processing chromosome: $chr
  echo number of cpus $task.cpus
  """
}

chr_ch = channel.of(1..22, 'X', 'Y')

workflow {
  PRINTCHR(chr_ch)
}
```

```output 
processing chromosome: 1
number of cpus 1

processing chromosome: 2
number of cpus 1

processing chromosome: 6
number of cpus 1
[..truncated..]
```

The above process uses the three directives, `tag`, `cpus` and `echo`.

The `tag` directive to allow you to give a custom tag to each process execution. This tag makes it easier to identify a particular task (executed instance of a process) in a log file or in the execution report.

The second directive `cpus`  allows you to define the number of CPUs required for each task.

The third directive `echo true` prints the stdout to the terminal.

We use the Nextflow `task.cpus` variable to capture the number of cpus assigned to a task. This is frequently used to specify the number of threads in a multi-threaded command in the script block.

Another commonly used directive is memory specification: `memory`.

A complete list of directives is available at this [link](https://www.nextflow.io/docs/latest/process.html#directives).

:::::::::::::::::::::::::::::::::::::::  challenge

## Adding directives


Many software tools allow users to configure the number of CPU threads used, optimizing performance for faster and more efficient data processing in high-throughput sequencing tasks.

In this exercise, we will use the bioinformatics tool [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of high-throughput sequencing read data. FastQC generates an HTML report along with a directory containing detailed analysis results. We can specify the number of CPU threads for FastQC to use with the -t option, followed by the desired number of threads.

Modify the Nextflow script `process_exercise_directives.nf`

1. Add a `tag` directive logging the sample_id in the execution output.
2. Add a `cpus` directive to specify the number of cpus as 2.
3. Change the fastqc `-t` option value to `$task.cpus` in the script directive.

```groovy 
//process_exercise_directives.nf


process FASTQC {
  //add tag directive
  //add cpu directive
 
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("fastqc_out")
  
  script:
  """
  mkdir fastqc_out
  fastqc $reads -o fastqc_out -t 1
  """
}

read_pairs_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')

workflow {
  FASTQC(read_pairs_ch)
  FASTQC.out.view()
}
```

:::::::::::::::  solution

## solution

```groovy 
//process_directives_answer.nf


process FASTQC {
  tag "$sample_id"
  cpus 2
  
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("fastqc_out")
 
  script:
  """
  mkdir fastqc_out
  fastqc $reads -o fastqc_out -t $task.cpus
  """
}

read_pairs_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')

workflow {
  FASTQC(read_pairs_ch)
  FASTQC.out.view()
}
```

```output
N E X T F L O W  ~  version 21.04.0
Launching `process_exercise_directives.nf` [sad_rosalind] - revision: 2ccbfa4937
executor >  local (3)
[90/de1125] process > FASTQC (ref1) [100%] 3 of 3 ✔
[ref2, work/ea/9e6a341b88caf8879e8d18b77049c8/fastqc_out]
[ref3, work/94/d059b816a9ec3d868f2924c26813e7/fastqc_out]
[ref1, work/90/de11251d362f494d6650789d9f8c1d/fastqc_out]
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Organising outputs

### PublishDir directive

Nextflow manages intermediate results from the pipeline's expected outputs independently.

Files created by a `process` are stored in a task specific working directory which is considered as temporary. Normally this is under the `work` directory, which can be deleted upon completion.

The files you want the workflow to return as results need to be defined in the `process` `output` block and then the output directory specified using the `directive` `publishDir`. More information [here](https://www.nextflow.io/docs/latest/process.html#publishdir).

**Note:** A common mistake is to specify an output directory in the `publishDir` directive while forgetting to specify the files you want to include in the `output` block.

```
publishDir <directory>, parameter: value, parameter2: value ...
```

For example if we want to capture the results of the `COMBINE_READS` process in a `results/merged_reads` output directory we
need to define the files in the `output` and  specify the location of the results directory in the `publishDir` directive:

```groovy 
//process_publishDir.nf


process COMBINE_READS {
  publishDir "results/merged_reads"

  input:
  tuple val(sample_id), path(reads)

  output:
  path("${sample_id}.merged.fq.gz")

  script:
  """
  cat ${reads} > ${sample_id}.merged.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')


workflow {
  COMBINE_READS(reads_ch)
}

```

```bash 
$ nextflow run process_publishDir.nf
```

```output 
N E X T F L O W  ~  version 21.04.0
Launching `process_publishDir.nf` [friendly_pauling] - revision: 9b5c315893
executor >  local (1)

[a1/5956bd] process > COMBINE_READS (1) [100%] 1 of 1 ✔
```

We can use the UNIX command `ls -l` to examine the contents of the results directory.

```bash 
ls -l results/merged_reads/ref1.merged.fq.gz
```

```output 
results/merged_reads/ref1.merged.fq.gz -> /Users/nf-user/nf-training/work/a1/5956bd9a92f13694b3ada1941f0d2d/ref1.merged.fq.gz
```

In the above example, the `publishDir "results/merged_reads"`,  creates a symbolic link `->` to the output files specified by the process `merged_reads` to the directory path `results/merged_reads`.

A symbolic link, often referred to as a symlink, is a type of file that serves as a reference or pointer to another file or directory, allowing multiple access paths to the same resource without duplicating its actual data

::::::::::::::::::::::::::::::::::::::::  callout

## publishDir

The publishDir output is relative to the path the pipeline run has been launched. Hence, it is a good practice to use [implicit variables](https://www.nextflow.io/docs/latest/script.html?highlight=projectdir#script-implicit-variables) like `projectDir` to specify publishDir value.


::::::::::::::::::::::::::::::::::::::::::::::::::

### publishDir parameters

The `publishDir` directive can take optional parameters, for example the `mode` parameter can take the value `"copy"` to specify that you wish to copy the file to output directory rather than just a symbolic link to the files in the working directory. Since the working directory is generally deleted on completion of a pipeline, it is safest to use `mode: "copy"` for results files. The default mode (symlink) is helpful for checking intermediate files which are not needed in the long term.

```groovy 
publishDir "results/merged_reads", mode: "copy"
```

Full list [here](https://www.nextflow.io/docs/latest/process.html#publishdir).

### Manage semantic sub-directories

You can use more than one `publishDir` to keep different outputs in separate directories. To specify which files to put in which output directory use the parameter `pattern` with the a glob pattern that selects which files to publish from the overall set of output files.

In the example below we will create an output folder structure in the directory results, which contains a separate sub-directory for sequence id file, `pattern:"*_ids.txt"` ,  and a sequence directory, `results/sequence"` for the `sequence.txt` file. Remember, we need to specify the files we want to copy as outputs.

```groovy 
//process_publishDir_semantic.nf


params.transcriptome="${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

process SPLIT_FASTA {
  publishDir "results/ids", pattern: "*_ids.txt", mode: "copy"
  publishDir "results/sequence", pattern: "sequence.txt", mode: "copy"


  input:
  path transcriptome

  output:
  path "*"

  script:
  """
  zgrep  '^>' $transcriptome > sequence_ids.txt
  zgrep -v '^>' $transcriptome > sequence.txt
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
transcriptome_ch = channel.fromPath(params.transcriptome)

workflow {
  SPLIT_FASTA(transcriptome_ch)
  // use the view operator to display contents of the channel
  SPLIT_FASTA.out.view()
}
```

```bash
$ nextflow run process_publishDir_semantic.nf
```

```output
N E X T F L O W  ~  version 21.04.0
Launching `process_publishDir_semantic.nf` [golden_poisson] - revision: 421a604840

executor >  local (1)
[be/950786] process > SPLIT_FASTA (1) [100%] 1 of 1 ✔
```

We can now use the `ls results/*` command to examine the results directory.

```bash 
$ls results/*
```

```output
results/ids:
sequence_ids.txt

results/sequence:
sequence.txt
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Publishing results

Add a `publishDir` directive to the nextflow script `process_exercise_publishDir.nf` that copies the merged reads  to the results folder merged\_reps.

```groovy 
//process_exercise_publishDir.nf


params.reads= "data/yeast/reads/ref{1,2,3}*{1,2}.fq.gz"

process MERGE_REPS {
 
 input:
 tuple val(sample_id), path(reads)
 
 output:
 path("*fq.gz")

 script:
 """
 cat *1.fq.gz > ${sample_id}.merged.R1.fq.gz
 cat *2.fq.gz > ${sample_id}.merged.R2.fq.gz
 """
}
reads_ch = Channel.fromFilePairs(params.reads,checkIfExists:true,size:6)

workflow {
 MERGE_REPS(reads_ch)
}
```

:::::::::::::::  solution

## Solution

```groovy 
//process_exercise_publishDir_answer.nf


params.reads= "data/yeast/reads/ref{1,2,3}*{1,2}.fq.gz"

process MERGE_REPS {
  publishDir "results/merged_reps"
  input:
  tuple val(sample_id), path(reads)
  output:
  path("*fq.gz")

  script:
  """
  cat *1.fq.gz > ${sample_id}.merged.R1.fq.gz
  cat *2.fq.gz > ${sample_id}.merged.R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs(params.reads,checkIfExists:true,size:6)

workflow {
  MERGE_REPS(reads_ch)
}
```

```bash 
$ nextflow run process_exercise_publishDir.nf
```

```output 
N E X T F L O W  ~  version 21.04.0
Launching `process_exercise_publishDir.nf` [infallible_becquerel] - revision: 4d865241a8

executor >  local (1)
[22/88aa22] process > MERGE_REPS (1) [100%] 1 of 1 ✔
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Nextflow Patterns

If you want to find out common structures of Nextflow processes, the [Nextflow Patterns page](https://nextflow-io.github.io/patterns/) collects some recurrent implementation patterns used in Nextflow applications.


::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::: keypoints

- Outputs to a process are defined using the output blocks.
- You can group input and output data from a process using the tuple qualifier.
- The execution of a process can be controlled using the `when` declaration and conditional statements.
- Files produced within a process and defined as `output` can be saved to a directory using the `publishDir` directive.

::::::::::::::::::::::::::::::::::::::::::::::::::


