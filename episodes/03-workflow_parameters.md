---
title: "Workflow parameterisation"
teaching: 20
exercises: 5
---

::::::::::::::::::::::::::::::::::::::: objectives

- "Use pipeline parameters to change the input to a workflow."
- "Add a pipeline parameters to a Nextflow script."
- "Understand how to create and use a parameter file."

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- "How can I change the data a workflow uses?"
- "How can I parameterise a workflow?"
- "How can I add my parameters to a file?"

::::::::::::::::::::::::::::::::::::::::::::::::::



In the first episode we ran the Nextflow script, `wc.nf`, from the
command line and it counted the number of lines in the file
`data/yeast/reads/ref1_1.fq.gz`. To change the input to script we can
make use of pipeline parameters.

## Pipeline parameters

The Nextflow `wc.nf` script defines a pipeline parameter `params.input`.
Pipeline parameters enable you to change the input to the workflow at
runtime, via the command line or a configuration file, so they are not
hard-coded into the script.

Pipeline parameters are declared in the workflow by prepending the
prefix `params`, separated by the dot character, to a variable name
e.g., `params.input`.

Their value can be specified on the command line by prefixing the
parameter name with a **double dash** character, e.g., `--input`.

In the script `wc.nf` the pipeline parameter `params.input` was
specified with a value of `"data/yeast/reads/ref1_1.fq.gz"`.

To process a different file, e.g. `data/yeast/reads/ref2_2.fq.gz`, in
the `wc.nf` script we would run:

``` bash
$ nextflow run wc.nf --input 'data/yeast/reads/ref2_2.fq.gz'
```

``` output
N E X T F L O W  ~  version 21.04.0
Launching `wc.nf` [gigantic_woese] - revision: 8acb5cb9b0
executor >  local (1)
[26/3cf986] process > NUM_LINES (1) [100%] 1 of 1 ✔
ref2_2.fq.gz 81720
```

We can also use wild cards to specify multiple input files (This will be
covered in the channels episode). In the example below we use the `*` to
match any sequence of characters between `ref2_` and `.fq.gz`. **Note:**
If you use wild card characters on the command line you must enclose the
value in quotes.

``` bash
$ nextflow run wc.nf --input 'data/yeast/reads/ref2_*.fq.gz'
```

This runs the process NUM_LINES twice, once for each file it matches.

``` output
N E X T F L O W  ~  version 21.04.0
Launching `wc.nf` [tender_lumiere] - revision: 8acb5cb9b0
executor >  local (2)
[cc/b6f793] process > NUM_LINES (1) [100%] 2 of 2 ✔
ref2_2.fq.gz 81720

ref2_1.fq.gz 81720
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Change a pipeline's input using a parameter

Re-run the Nextflow script `wc.nf` by changing the pipeline input to all
files in the directory `data/yeast/reads/` that begin with `ref` and end
with `.fq.gz`:

:::::::::::::::  solution

## Solution

``` bash
$ nextflow run wc.nf --input 'data/yeast/reads/ref*.fq.gz'
```

The string specified on the command line will override the default value
of the parameter in the script. The output will look like this:

``` output
N E X T F L O W  ~  version 20.10.0
Launching `wc.nf` [soggy_miescher] - revision: c54a707593
executor >  local (6)
[d3/9ca185] process > NUM_LINES (2) [100%] 6 of 6 ✔
ref3_2.fq.gz 52592

ref2_2.fq.gz 81720

ref1_1.fq.gz 58708

ref1_2.fq.gz 58708

ref3_1.fq.gz 52592

ref2_1.fq.gz 81720
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


## Adding a parameter to a script

To add a pipeline parameter to a script prepend the prefix `params`,
separated by a dot character `.`, to a variable name e.g.,
`params.input`.

Let's make a copy of the `wc.nf` script as `wc-params.nf` and add a new
input parameter.

``` bash
$ cp wc.nf wc-params.nf
```

To add a parameter `sleep` with the default value `2` to `wc-params.nf`
we add the line:

``` groovy
params.sleep = 2
```

**Note:** You should always add a sensible default value to the pipeline
parameter. We can use this parameter to add another step to our
`NUM_LINES` process.

``` groovy
script:
 """
 sleep ${params.sleep}
 printf '${read} '
 gunzip -c ${read} | wc -l
 """
```

This step, `sleep ${params.sleep}`, will add a delay for the amount of
time specified in the `params.sleep` variable, by default 2 seconds. To
access the value inside the script block we use `{variable_name}` syntax
e.g. `${params.sleep}`.

We can now change the sleep parameter from the command line, For
Example:

``` bash
nextflow run wc-params.nf --sleep 10
```



:::::::::::::::::::::::::::::::::::::::  challenge


## Add a pipeline parameter

If you haven't already make a copy of the `wc.nf` as `wc-params.nf`.

```bash
$ cp wc.nf wc-params.nf
```

Add the param `sleep` with a default value of 2 below the `params.input` line. 
Add the line `sleep ${params.sleep}` in the process `NUM_LINES` above the line
printf `${read}`.

Run the new script `wc-params.nf` changing the sleep input time.

What input file would it run and why?

How would you get it to process all `.fq.gz` files in the `data/yeast/reads` 
directory as well as changing the sleep input to 1 second?

:::::::::::::::  solution

## Solution

```groovy
params.sleep=2
```

```groovy 
script: 
"""
sleep ${params.sleep} > printf '${read} '
gunzip -c ${read} | wc -l 
"""
```

```bash
$ nextflow run wc-params.nf --sleep 1 
```

This would use 1 as a value of `sleep` parameter instead of default value (which is 2) and run the pipeline. 
The input file would be  `data/yeast/reads/ref1_1.fq.gz` as this is the default.
To run all input files we could add the param
`--input 'data/yeast/reads/*.fq.gz'` 
```bash
$ nextflow run wc-params.nf --sleep 1 --input 'data/yeast/reads/\*.fq.gz' 
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Parameter File

If we have many parameters to pass to a script it is best to create a
parameters file. Parameters are stored in JSON or YAML format. JSON and
YAML are data serialization languages, that are a way of storing data
objects and structures, such as the `params` object in a file.

The `-params-file` option is used to pass the parameters file to the
script.

For example the file `wc-params.json` contains the parameters `sleep`
and `input` in JSON format.

```json         
{
  "sleep": 5,
  "input": "data/yeast/reads/etoh60_1*.fq.gz"
}
```


To run the `wc-params.nf` script using these parameters we add the
option `-params-file` and pass the file `wc-params.json`:

```bash         
$ nextflow run wc-params.nf -params-file wc-params.json
```



```output         
N E X T F L O W  ~  version 21.04.0
Launching `wc-params.nf` [nostalgic_northcutt] - revision: 2f86c9ac7e
executor >  local (2)
[b4/747eaa] process > NUM_LINES (1) [100%] 2 of 2 ✔
etoh60_1_2.fq.gz 87348

etoh60_1_1.fq.gz 87348
```


:::::::::::::::::::::::::::::::::::::::  challenge

## Create and use a Parameter file.

Create a parameter file `params.json` for the Nextflow file
`wc-params.nf`, and run the Nextflow script using the created
parameter file, specifying:

-   sleep as 10
-   input as `data/yeast/reads/ref3_1.fq.gz`

:::::::::::::::  solution

## Solution

```json         
{
"sleep": 10,
"input": "data/yeast/reads/ref3_1.fq.gz"

```bash
$ nextflow run wc-params.nf -params-file params.json 
```

```output
N E X T F L O W 
 version 21.04.0 Launching `wc-params.nf` [small_wiles] - revision:
f5ef7b7a01 executor \> local (1) [f3/4fa480] process \> NUM_LINES
(1) [100%] 1 of 1 ✔ ref3_1.fq.gz 52592 
```
:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::::::::: keypoints
- "Pipeline parameters are specified by prepending the prefix `params` to a variable name, separated by dot character."
- "To specify a pipeline parameter on the command line for a Nextflow run use `--variable_name` syntax."
- "You can add parameters to a JSON or YAML formatted file and pass them to the script using option `-params-file`."
::::::::::::::::::::::::::::::::::::::::::::::::::
