---
title: Deploying nf-core pipelines
teaching: 30
exercises: 10
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand what nf-core is and how it relates to Nextflow.
- Use the nf-core helper tool to find nf-core pipelines.
- Understand how to configuration nf-core pipelines.
- Run a small demo nf-core pipeline using a test dataset.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- Where can I find best-practice Nextflow bioinformatic pipelines?
- How do I run nf-core pipelines?
- How do I configure nf-core pipelines to use my data?
- How do I reference nf-core pipelines?

::::::::::::::::::::::::::::::::::::::::::::::::::

### What is nf-core?

nf-core is a community-led project to develop a set of best-practice pipelines built using Nextflow workflow management system.
Pipelines are governed by a set of guidelines, enforced by community code reviews and automatic code testing.

![A diagram showcasing the key aspects of nf-core, a community effort to provide best-practice analysis pipelines. The diagram is divided into three sections: Deploy, Participate, and Develop. The Deploy section includes features like Stable pipelines, Centralized configs, List and update pipelines, and Download for offline use. The Participate section highlights Documentation, Slack workspace, Twitter updates, and Hackathons. The Develop section emphasizes the Starter template, Code guidelines, CI code linting and tests, and Helper tools.](fig/nf-core.png 'nf-core')

In this episode we will covering finding, deploying and configuring nf-core pipelines.

### What are nf-core pipelines?

nf-core pipelines are an organised collection of Nextflow scripts,  other non-nextflow scripts (written in any language), configuration files, software specifications, and documentation hosted on [GitHub](https://github.com/nf-core). There is generally a single pipeline for a given data and analysis type e.g. There is a single pipeline for bulk RNA-Seq. All nf-core pipelines are distributed under the, permissive free software, [MIT licences](https://en.wikipedia.org/wiki/MIT_License).

### What is nf-core tools?

nf-core provides a suite of helper tools aim to help people run and develop pipelines.
The [nf-core tools](https://nf-co.re/tools) package is written in Python and can run from the command line or imported and used within other packages.

#### nf-core tools sub-commands

You can use the `--help` option to see the range of nf-core tools sub-commands.
In this episode we will be covering the `list`, `launch` and `download` sub-commands which
aid in the finding and deployment of the nf-core pipelines.

```bash
nf-core --help
```

```output
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

 nf-core/tools version 2.14.1 - https://nf-co.re



 Usage: nf-core [OPTIONS] COMMAND [ARGS]...

 nf-core/tools provides a set of helper tools for use with nf-core Nextflow pipelines.
 It is designed for both end-users running pipelines and also developers creating new pipelines.

╭─ Options ────────────────────────────────────────────────────────────────────────────────────────╮
│ --version                        Show the version and exit.                                      │
│ --verbose        -v              Print verbose output to the console.                            │
│ --hide-progress                  Don't show progress bars.                                       │
│ --log-file       -l  <filename>  Save a verbose log to a file.                                   │
│ --help           -h              Show this message and exit.                                     │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands for users ─────────────────────────────────────────────────────────────────────────────╮
│ list                  List available nf-core pipelines with local info.                          │
│ launch                Launch a pipeline using a web GUI or command line prompts.                 │
│ create-params-file    Build a parameter file for a pipeline.                                     │
│ download              Download a pipeline, nf-core/configs and pipeline singularity images.      │
│ licences              List software licences for a given workflow (DSL1 only).                   │
│ tui                   Open Textual TUI.                                                          │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands for developers ────────────────────────────────────────────────────────────────────────╮
│ create            Create a new pipeline using the nf-core template.                              │
│ lint              Check pipeline code against nf-core guidelines.                                │
│ modules           Commands to manage Nextflow DSL2 modules (tool wrappers).                      │
│ subworkflows      Commands to manage Nextflow DSL2 subworkflows (tool wrappers).                 │
│ schema            Suite of tools for developers to manage pipeline schema.                       │
│ create-logo       Generate a logo with the nf-core logo template.                                │
│ bump-version      Update nf-core pipeline version number.                                        │
│ sync              Sync a pipeline TEMPLATE branch with the nf-core template.                     │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
```

### Listing available nf-core pipelines

The simplest sub-command is `nf-core list`, which lists all available nf-core pipelines in the nf-core Github repository.

The output shows the latest version number and when that was released.
If the pipeline has been pulled locally using Nextflow, it tells you when that was and whether you have the latest version.

Run the command below.

```bash
nf-core list
```

An example of the output from the command is as follows:

```output


                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.14.1 - https://nf-co.re


┏━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name             ┃ Stars ┃ Latest Release ┃      Released ┃ Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ mcmicro                   │     4 │            dev │    2 days ago │           - │ -                    │
│ fastquorum                │    13 │          1.0.0 │  2 months ago │           - │ -                    │
│ rnaseq                    │   821 │         3.14.0 │  6 months ago │           - │ -                    │
│ crisprseq                 │    22 │          2.2.0 │  1 months ago │           - │ -                    │
│ funcscan                  │    62 │          1.1.6 │   2 weeks ago │           - │ -                    │
│ pairgenomealign           │     0 │            dev │    3 days ago │           - │ -                    │
│ multiplesequencealign     │    11 │            dev │    3 days ago │           - │ -                    │
│ denovotranscript          │     0 │            dev │    3 days ago │           - │ -                    │
│ demo                      │     1 │          1.0.0 │  1 months ago │           - │ -                    │
│ demultiplex               │    37 │          1.4.1 │  5 months ago │           - │ -                    │
[..truncated..]
```

#### Filtering available nf-core pipelines

If you supply additional keywords after the `list` sub-command, the listed pipeline will be filtered.
**Note:** that this searches more than just the displayed output, including keywords and description text.

Here we filter on the keywords **rna** and **rna-seq** .

```bash
nf-core list rna rna-seq
```

```output

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.14.1 - https://nf-co.re


┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name         ┃ Stars ┃ Latest Release ┃     Released ┃ Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ rnaseq                │   821 │         3.14.0 │ 6 months ago │           - │ -                    │
│ denovotranscript      │     0 │            dev │   3 days ago │           - │ -                    │
│ scnanoseq             │     5 │            dev │   4 days ago │           - │ -                    │
│ circrna               │    43 │            dev │   6 days ago │           - │ -                    │
│ smrnaseq              │    70 │          2.3.1 │ 3 months ago │           - │ -                    │
│ scrnaseq              │   185 │          2.7.0 │  4 weeks ago │           - │ -                    │
│ differentialabundance │    47 │          1.5.0 │ 2 months ago │           - │ -                    │
│ rnafusion             │   132 │          3.0.2 │ 3 months ago │           - │ -                    │
│ spatialvi             │    45 │            dev │ 2 months ago │           - │ -                    │
│ rnasplice             │    35 │          1.0.4 │ 2 months ago │           - │ -                    │
│ dualrnaseq            │    16 │          1.0.0 │  3 years ago │           - │ -                    │
│ marsseq               │     4 │          1.0.3 │  1 years ago │           - │ -                    │
│ lncpipe               │    30 │            dev │  2 years ago │           - │ -                    │
└───────────────────────┴───────┴────────────────┴──────────────┴─────────────┴──────────────────────┘
```

#### Sorting available nf-core pipelines

You can sort the results by adding the option `--sort` followed by a keyword.
For example, latest release (`--sort release`), when you last pulled a local copy (`--sort pulled`), alphabetically (`--sort name`), or number of GitHub stars (`--sort stars`).

```bash
nf-core list rna rna-seq --sort stars
```

```output
                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--~\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

    nf-core/tools version 2.14.1 - https://nf-co.re


┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name         ┃ Stars ┃ Latest Release ┃     Released ┃ Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ rnaseq                │   821 │         3.14.0 │ 6 months ago │           - │ -                    │
│ scrnaseq              │   185 │          2.7.0 │  4 weeks ago │           - │ -                    │
│ rnafusion             │   132 │          3.0.2 │ 3 months ago │           - │ -                    │
│ smrnaseq              │    70 │          2.3.1 │ 3 months ago │           - │ -                    │
│ differentialabundance │    47 │          1.5.0 │ 2 months ago │           - │ -                    │
│ spatialvi             │    45 │            dev │ 2 months ago │           - │ -                    │
│ circrna               │    43 │            dev │   6 days ago │           - │ -                    │
│ rnasplice             │    35 │          1.0.4 │ 2 months ago │           - │ -                    │
│ lncpipe               │    30 │            dev │  2 years ago │           - │ -                    │
│ dualrnaseq            │    16 │          1.0.0 │  3 years ago │           - │ -                    │
│ scnanoseq             │     5 │            dev │   4 days ago │           - │ -                    │
│ marsseq               │     4 │          1.0.3 │  1 years ago │           - │ -                    │
│ denovotranscript      │     0 │            dev │   3 days ago │           - │ -                    │
└───────────────────────┴───────┴────────────────┴──────────────┴─────────────┴──────────────────────┘
```

:::::::::::::::::::::::::::::::::::::::::  callout

### Archived pipelines

Archived pipelines are not returned by default. To include them, use the `--show_archived` flag.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise: listing nf-core pipelines

1. Use the `--help` flag to print the list command usage.
2. Sort all pipelines by popularity (stars) and find out which is the most popular?.
3. Filter pipelines for those that work with RNA and find out how many there are?

:::::::::::::::  solution

### Solution

Use the `--help` flag to print the list command usage.

```bash
$ nf-core list --help
```

Sort all pipelines by popularity (stars).

```bash
$ nf-core list --sort stars
```

Filter pipelines for those that work with RNA.

```bash
$ nf-core list rna
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Running nf-core pipelines

#### Software requirements for nf-core pipelines

nf-core pipeline software dependencies are specified using either Docker, Singularity or Conda. It is Nextflow that handles the downloading of containers and creation of conda environments. In theory it is possible to run the pipelines with software installed by other methods (e.g. environment modules, or manual installation), but this is not recommended.

#### Fetching pipeline code

Unless you are actively developing pipeline code, you should use Nextflow's [built-in functionality](https://www.nextflow.io/docs/latest/sharing.html) to fetch nf-core pipelines. You can use the following command to pull the latest version of a remote workflow from the nf-core github site.;

```bash
$ nextflow pull nf-core/<PIPELINE>
```

**Note** Nextflow will also automatically fetch the pipeline code when use `nextflow run nf-core/<PIPELINE>` command.


For the best reproducibility, it is good to explicitly reference the pipeline version number that you wish to use with the `-revision`/`-r` flag.

In the example below we are pulling the rnaseq pipeline version 3.0

```bash
nextflow pull nf-core/rnaseq -revision 3.14.0
```

We can check the pipeline has been pulled using the `nf-core list` command.

```bash
nf-core list rnaseq -s pulled
```

We can see from the output we have the latest release.

```output
                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--~\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

    nf-core/tools version 2.14.1 - https://nf-co.re


┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Pipeline Name         ┃ Stars ┃ Latest Release ┃     Released ┃    Last Pulled ┃ Have latest release? ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│ rnaseq                │   821 │         3.14.0 │ 6 months ago │ 59 seconds ago │ Yes (v3.14.0)        │
[..truncated..]
```

:::::::::::::::::::::::::::::::::::::::::  callout

### Development Releases

If not specified, Nextflow will fetch the default git branch. For pipelines with a stable release this the default branch is `master` - this branch contains code from the latest release. For pipelines in early development that don't have any releases, the default branch is `dev`.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise: Fetch the latest RNA-Seq pipeline

1. Use the `nextflow pull` command to download the latest `nf-core/rnaseq` pipeline

2. Use the `nf-core list` command to see if you have the latest version of the pipeline

:::::::::::::::  solution

### Solution

Use the `nextflow pull` command to download the latest `nf-core/rnaseq` pipeline

```bash
$ nextflow pull nf-core/rnaseq
```

Use the `nf-core list` command to see if you have the latest version of the pipeline

```bash
$ nf-core list rnaseq --sort pulled
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Usage instructions and documentation

You can find general documentation and instructions for Nextflow and nf-core on the [nf-core website](https://nf-co.re/) .
Pipeline-specific documentation is bundled with each pipeline in the /docs folder. This can be read either locally, on GitHub, or on the nf-core website.

Each pipeline has its own webpage at [https://nf-co.re/](https://nf-co.re/)\<pipeline\_name> e.g. [nf-co.re/rnaseq](https://nf-co.re/rnaseq/usage)

In addition to this documentation, each pipeline comes with basic command line reference. This can be seen by running the pipeline with the parameter `--help` , for example:

```bash
nextflow run -r 3.14.0 nf-core/rnaseq --help
```

```output
N E X T F L O W  ~  version 20.10.0
Launching `nf-core/rnaseq` [silly_miescher] - revision: 964425e3fd [3.4]
------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/rnaseq v3.14.0-gb89fac3
------------------------------------------------------
Typical pipeline command:

  nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker

Input/output options
  --input                            [string]  Path to comma-separated file containing information about the samples in the experiment.
  --outdir                           [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                               infrastructure.
  --email                            [string]  Email address for completion summary.
  --multiqc_title                    [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.
..truncated..
```

### The nf-core launch command

As can be seen from the output of the help option nf-core pipelines have a number of flags that need to be passed on the command line: some mandatory, some optional.

To make it easier to launch pipelines, these parameters are described in a JSON file, `nextflow_schema.json` bundled with the pipeline.

The `nf-core launch` command uses this to build an interactive command-line wizard which walks through the different options with descriptions of each, showing the default value and prompting for values.

Once all prompts have been answered, non-default values are saved to a `params.json` file which can be supplied to Nextflow  using the `-params-file` option. Optionally, the Nextflow command can be launched there and then.

To use the launch feature, just specify the pipeline name:

```bash
nf-core launch -r 3.14.0 rnaseq
```

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise : Create nf-core params file

Use the `nf-core launch` command to create a params file named `nf-params.json`.

1. Use the `nf-core launch` command to launch the interactive command-line wizard.
2. Add an input file name `samples.csv`
3. Add a genome `GRCh38`
  \*\* Note \*\* : Do not run the command now.

:::::::::::::::  solution

### Solution

```bash
$ nf-core launch rnaseq
```

The contents of the `nf-params.json` file should be

```
{
  "input": "samples.csv",
  "genome": "GRCh38"
}
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Config files

nf-core pipelines make use of Nextflow's configuration files to specify how the pipelines runs, define custom parameters and what software management system to use e.g. docker, singularity or conda.

Nextflow can load pipeline configurations from multiple locations.  nf-core pipelines load configuration in the following order:

![A diagram explaining the structure and hierarchy of nextflow.config files. The diagram shows that the default 'base' config is always loaded. It also includes core profiles, such as docker, conda, and test, and server profiles from nf-core/config. Additionally, it highlights that your local config files, located in $HOME/.nextflow/config or specified with -c custom.config, are also considered.](fig/nfcore_config.png 'config')

1. Pipeline: Default 'base' config

- Always loaded. Contains pipeline-specific parameters and "sensible defaults" for things like computational requirements
- Does not specify any method for software packaging. If nothing else is specified, Nextflow will expect all software to be available on the command line.

2. Core config profiles

- All nf-core pipelines come with some generic config profiles. The most commonly used ones are for software packaging: docker, singularity and conda
- Other core profiles are debug and two test profiles. There two test profile, a small test profile (nf-core/test-datasets) for quick test and a full test profile which provides the path to full sized data from public repositories.

3. Server profiles

- At run time, nf-core pipelines fetch configuration profiles from the [configs remote repository](https://github.com/nf-core/configs). The profiles here are specific to clusters at different institutions.
- Because this is loaded at run time, anyone can add a profile here for their system and it will be immediately available for all nf-core pipelines.

4. Local config files given to Nextflow with the `-c` flag

```bash
nextflow run nf-core/rnaseq -r 3.14.0 -c mylocal.config
```

5\. Command line configuration: pipeline parameters can be passed on the command line using the `--<parameter>` syntax.

```bash
nextflow run nf-core/rnaseq -r 3.14.0 --email "my@email.com"`
```

#### Config Profiles

nf-core makes use of Nextflow configuration `profiles` to make it easy to apply a group of options on the command line.

Configuration files can contain the definition of one or more profiles. A profile is a set of configuration attributes that can be activated/chosen when launching a pipeline execution by using the `-profile` command line option. Common profiles are `conda`, `singularity` and `docker` that specify which software manager to use.

Multiple profiles are comma-separated. When there are differing configuration settings provided by different profiles, the right-most profile takes priority.

```bash
nextflow run nf-core/rnaseq -r 3.14.0 -profile test,conda
nextflow run nf-core/rnaseq -r 3.14.0 -profile <institutional_config_profile>, test, conda
```

**Note** The order in which config profiles are specified matters. Their priority increases from left to right.

:::::::::::::::::::::::::::::::::::::::::  callout

### Multiple Nextflow configuration locations

Be clever with multiple Nextflow configuration locations. For example, use `-profile` for your cluster  configuration, the file `$HOME/.nextflow/config` for your personal config such as `params.email` and a working directory >`nextflow.config` file for reproducible run-specific configuration.


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise  create a custom config

Add the `params.email`  to a file called `nfcore-custom.config`

:::::::::::::::  solution

### Solution

A line similar to one below in the file custom.config

```
params.email = "myemail@address.com"
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Running pipelines with test data

The nf-core config profile `test` is special profile, which defines a minimal data set and configuration, that runs quickly and tests the workflow from beginning to end. Since the data is minimal, the output is often nonsense. Real world  example output are instead linked on the nf-core pipeline web page, where the workflow has been run with a full size data set:

```bash
nextflow run nf-core/<pipeline_name> -profile test
```

:::::::::::::::::::::::::::::::::::::::::  callout

### Software configuration profile

Note that you will typically still need to combine this with a software configuration profile for your system - e.g.
`-profile test,conda`.
Running with the test profile is a great way to confirm that you have Nextflow configured properly for your system before attempting to run with real data


::::::::::::::::::::::::::::::::::::::::::::::::::

### Using nf-core pipelines offline

Many of the techniques and resources described above require an active internet connection at run time - pipeline files, configuration profiles and software containers are all dynamically fetched when the pipeline is launched. This can be a problem for people using secure computing resources that do not have connections to the internet.

To help with this, the `nf-core download` command automates the fetching of required files for running nf-core pipelines offline.
The command can download a specific release of a pipeline with `-r`/`--release` .  
By default, the pipeline will download the pipeline code and the institutional nf-core/configs files.

If you specify the flag `--singularity`, it will also download any singularity image files that are required (this needs Singularity to be installed). All files are saved to a single directory, ready to be transferred to the cluster where the pipeline will be executed.

```bash
nf-core download nf-core/rnaseq -r 3.14.0
```

```output

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.14.1 - https://nf-co.re


WARNING  Could not find GitHub authentication token. Some API requests may fail.
? Include the nf-core's default institutional configuration files into the download? No

In addition to the pipeline code, this tool can download software containers.
? Download software container images: none

If transferring the downloaded files to another system, it can be convenient to have everything compressed in a single file.
? Choose compression type: none
INFO     Saving 'nf-core/rnaseq'
          Pipeline revision: '3.14.0'
          Use containers: 'none'
          Container library: 'quay.io'
          Output directory: 'nf-core-rnaseq_3.14.0'
          Include default institutional configuration: 'False'
INFO     Downloading workflow files from GitHub
```

:::::::::::::::::::::::::::::::::::::::  challenge

### Exercise  Run a demo nf-core pipeline

Run the `nf-core/demo` pipeline release 1.0.0  with the provided test data using the profile `test`.
Add the parameters `--max_memory 3G`, `--skip_trim`, and --outdir "nfcore-demo-out" on the command line.
Include the config file, `nfcore-custom.config`, from the previous exercise using the option `-c`, to send an email when your pipeline finishes.

```bash
$ nextflow run nf-core/demo -r 1.0.0  -profile test --skip_trim --max_memory 3.GB -c nfcore-custom.config
```

The `nf-core/demo` pipleine is a simple nf-core style bioinformatics pipeline for workshops and demonstrations that runs FASTQC and multiqc.

:::::::::::::::  solution

### Solution

```output
 N E X T F L O W   ~  version 24.04.3

Launching `https://github.com/nf-core/demo` [silly_bohr] DSL2 - revision: 705f18e4b1 [master]

WARN: Access to undefined parameter `monochromeLogs` -- Initialise it to a default value eg. `params.monochromeLogs = some_value`


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/demo v1.0.0-g705f18e
------------------------------------------------------
Core Nextflow options
  revision                  : master
  runName                   : silly_bohr
  launchDir                 : /home/nf-training
  workDir                   : /home/nf-training/work
  projectDir                : /home/nf-training/.nextflow/assets/nf-core/demo
  userName                  : nf-training
  profile                   : test
  configFiles               :

Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
  outdir                    : out
  email                     : myemail@address.com

Process skipping options
  skip_trim                 : true

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Max job request options
  max_cpus                  : 2
  max_memory                : 3.GB
  max_time                  : 6.h

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
If you use nf-core/demo for your analysis please cite:

* The pipeline

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/nf-core/demo/blob/master/CITATIONS.md
------------------------------------------------------
executor >  local (4)
[08/94286a] process > NFCORE_DEMO:DEMO:FASTQC (SAMPLE2_PE) [100%] 3 of 3 ✔
[df/45c5c7] process > NFCORE_DEMO:DEMO:MULTIQC             [100%] 1 of 1 ✔
-[nf-core/demo] Sent summary e-mail to graeme.grimes@ed.ac.uk (sendmail)-
-[nf-core/demo] Pipeline completed successfully-
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Troubleshooting

If you run into issues running your pipeline you can you the nf-core  website  to troubleshoot common mistakes and issues [https://nf-co.re/usage/troubleshooting](https://nf-co.re/usage/troubleshooting) .

#### Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact the nf-core community via the Slack channel .
The nf-core Slack organisation has channels dedicated for each pipeline, as well as specific topics (eg. `#help`, `#pipelines`, `#tools`, `#configs` and much more).
The nf-core Slack can be found at [https://nfcore.slack.com](https://nfcore.slack.com) (NB: no hyphen in nfcore!). To join you will need an invite, which you can get at [https://nf-co.re/join/slack](https://nf-co.re/join/slack).

You can also get help by opening an issue in the respective pipeline repository on GitHub asking for help.

If you have problems that are directly related to Nextflow and not our pipelines or the nf-core framework tools then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).

### Referencing a Pipeline

#### Publications

If you use an nf-core pipeline in your work you should cite the main publication for the main nf-core paper, describing the community and framework,
as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso \& Sven Nahnsen.
> Nat Biotechnol. 2020 Feb 13. doi: [10\.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x). ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

Many of the pipelines have a publication listed on the nf-core website that can be found [here](https://nf-co.re/publications).

#### DOIs

In addition, each release of an nf-core pipeline has a digital object identifiers (DOIs) for easy referencing in literature
The DOIs are generated by Zenodo from the pipeline's github repository.



:::::::::::::::::::::::::::::::::::::::: keypoints

- nf-core is a community-led project to develop a set of best-practice pipelines built using the Nextflow workflow management system.
- The nf-core tool (`nf-core`) is a suite of helper tools that aims to help people run and develop nf-core pipelines.
- nf-core pipelines can be found using `nf-core list`, or by checking the nf-core website.
- `nf-core launch nf-core/<pipeline>` can be used to write a parameter file for an nf-core pipeline. This can be supplied to the pipeline using the `-params-file` option.
- An nf-core workflow is run using `nextflow run nf-core/<pipeline>` syntax.
- nf-core pipelines can be reconfigured by using custom config files and/or adding command line parameters.

::::::::::::::::::::::::::::::::::::::::::::::::::


