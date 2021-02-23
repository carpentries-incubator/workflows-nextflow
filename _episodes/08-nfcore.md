---
title: "nf-core"
teaching: 20
exercises: 10
questions:
- "What is nf-core?"
- "How do you find nf-core pipelines?"
- "How do you run nf-core pipelines?"
- "How do you use nf-core pipelines offline?"
- "How do oyu configure nf-core pipleines for your envirnoment?"
objectives:
- "Understand what nf-core core is how it relates to Nextflow."
- "List and search nf-core pipelines using the nf-core helper tool ."
- "Run test nf-core RNA-Seq pipeline."
keypoints:
- "nf-core is a community-led project to develop a set of best-practice pipelines built using Nextflow. "
- "nf-core pipelines can found using the nf-core helper tool or from he nf-core website."
---


## What is nf-core

nf-core is a community-led project to develop a set of best-practice pipelines built using Nextflow. Pipelines are governed by a set of guidelines, enforced by community code reviews and automatic linting (code testing). A suite of helper tools aim to help people run and develop pipelines.

## Where to get help

The beauty of nf-core is that there is lots of help on offer! The main place for this is Slack - an instant messaging service.

The nf-core Slack can be found at https://nfcore.slack.com (NB: no hyphen in nfcore!). To join you will need an invite, which you can get at https://nf-co.re/join/slack.

The nf-core Slack organisation has channels dedicated for each pipeline, as well as specific topics (eg. #help, #pipelines, #tools, #configs and much more).

> ## TLDR
> One additional tool which we like a lot is [TLDR](https://tldr.sh/) - it gives concise command line reference through example commands for most linux tools, including nextflow, docker, singularity, conda, git and more. There are many clients, but [raylee/tldr](raylee/tldr ) is arguably the simplest - just a single bash script.
{: .callout}

## What is nf-core tools

The nf-core tools package is written in Python and can be imported and used within other packages. 

### Automatic version check

nf-core/tools automatically checks the web to see if there is a new version of nf-core/tools available. If you would prefer to skip this check, set the environment variable NFCORE_NO_VERSION_CHECK. For example:

~~~
export NFCORE_NO_VERSION_CHECK=1
~~~
{: .language-bash}

### nf-core tools sub-commands

You can use the `--help` option to see the range  of sub-commands.

~~~
nf-core --help
~~~


## Listing available nf-core pipelines

As you saw from the `--help` output, the tool has a range of sub-commands. The simplest is `nf-core list`, which lists all available nf-core pipelines. The output shows the latest version number, when that was released. If the pipeline has been pulled locally using Nextflow, it tells you when that was and whether you have the latest version.

~~~
nf-core list
~~~
{: .language-bash}



If you supply additional keywords after the command, the listed pipeline will be filtered. Note that this searches more than just the displayed output, including keywords and description text. 

~~~
nf-core list rna rna-seq
~~~
{: .language-bash}

The `--sort` flag allows you to sort the list (default is by most recently released) and `--json` gives JSON output for programmatic use.

> ## Exercise 2 listing pipelines
>
>  Use the `--help` flag to print the list command usage
>  Sort pipelines alphabetically, then by popularity (stars)
>  Fetch one of the pipelines using `nextflow pull`
>  Use `nf-core list` to see if the pipeline you pulled is up to date
>  Filter pipelines for those that work with RNA
>  ave these pipeline details to a JSON file
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

## Running nf-core pipelines

### Software requirements for nf-core pipelines

In order to run nf-core pipelines, you will need to have [Nextflow](https://www.nextflow.io) installed . The only other requirement is a software packaging tool: Conda, Docker or Singularity. In theory it is possible to run the pipelines with software installed by other methods (e.g. environment modules, or manual installation), but this is not recommended. Most people find either Docker or Singularity the best options.

## Fetching pipeline code

Unless you are actively developing pipeline code, we recommend using the Nextflow [built-in functionality](https://www.nextflow.io/docs/latest/sharing.html) to fetch `nf-core pipelines`. Nextflow will automatically fetch the pipeline code when you run nextflow `run nf-core/PIPELINE`. For the best reproducibility, it is good to explicitly reference the pipeline version number that you wish to use with the `-revision`/`-r` flag. For example:

~~~
nextflow run nf-core/rnaseq -revision 3.0
~~~
{: .language-bash}

If not specified, Nextflow will fetch the default branch. For pipelines with a stable release this the default branch is `master` - this branch contains code from the latest release. For pipelines in early development that don't have any releases, the default branch is `dev`.

> ## latest flag
> Note that once pulled, Nextflow will use the local cached version for subsequent runs. Use the -latest flag when running the pipeline to always fetch the latest version. Alternatively, you can force Nextflow to pull a pipeline again using the nextflow pull command:
> ~~~
>  nextflow pull nf-core/rnaseq
> ~~~
{: .callout}

# Usage instructions and documentation

You can find general documentation and instructions for Nextflow and nf-core on the nf-core website: https://nf-co.re/. Pipeline-specific documentation is bundled with each pipeline in the /docs folder. This can be read either locally, on GitHub, or on the nf-core website. Each pipeline has its own webpage at https://nf-co.re/<pipeline_name> (e.g. nf-co.re/rnaseq)

In addition to this documentation, each pipeline comes with basic command line reference. This can be seen by running the pipeline with the `--help` flag, for example:

~~~
nextflow run nf-core/rnaseq --help
~~~
{: .language-bash}

## Config profiles

Nextflow can load pipeline configurations from multiple locations. To make it easy to apply a group of options on the command line, Nextflow uses the concept of config profiles. nf-core pipelines load configuration in the following order:

1. Pipeline: Default 'base' config
* Always loaded. Contains pipeline-specific parameters and "sensible defaults" for things like computational requirements
* Does not specify any method for software packaging. If nothing else is specified, Nextflow will expect all software to be available on the command line.
1. Pipeline: Core config profiles
* All nf-core pipelines come with some generic config profiles. The most commonly used ones are for software packaging: docker, singularity and conda
* Other core profiles are debug and test
1. nf-core/configs: Server profiles
* At run time, nf-core pipelines fetch configuration profiles from the configs remote repository. The profiles here are specific to clusters at different institutions.
* Because this is loaded at run time, anyone can add a profile here for their system and it will be immediately available for all nf-core pipelines.
1. Local config files given to Nextflow with the -c flag
1. Command line configuration

Multiple comma-separate config profiles can be specified in one go, so the following commands are perfectly valid:

~~~
nextflow run nf-core/rnaseq -profile test,docker
nextflow run nf-core/rnaseq -profile singularity,debug
~~~
{: .language-bash}

Note that the order in which config profiles are specified matters. Their priority increases from left to right.

Our tip: Be clever with multiple Nextflow configuration locations. For example, use `-profile` for your cluster configuration, `~/.nextflow/config` for your personal config such as `params.email` and a working directory `nextflow.config` file for reproducible run-specific configuration.

# Running pipelines with test data

The test config profile is a bit of a special case. Whereas all other config profiles tell Nextflow how to run on different computational systems, the test profile configures each nf-core pipeline to run without any other command line flags. It specifies URLs for test data and all required parameters. Because of this, you can test any nf-core pipeline with the following command:

~~~
nextflow run nf-core/<pipeline_name> -profile test
~~~
{: .language-bash}

Note that you will typically still need to combine this with a configuration profile for your system - e.g. `-profile test,docker`. Running with the test profile is a great way to confirm that you have Nextflow configured properly for your system before attempting to run with real data

## The nf-core launch command

Most nf-core pipelines have a number of flags that need to be passed on the command line: some mandatory, some optional. To make it easier to launch pipelines, these parameters are described in a JSON file bundled with the pipeline. The nf-core launch command uses this to build an interactive command-line wizard which walks through the different options with descriptions of each, showing the default value and prompting for values.

Once all prompts have been answered, non-default values are saved to a `params.json` file which can be supplied to Nextflow to run the pipeline. Optionally, the Nextflow command can be launched there and then.

To use the launch feature, just specify the pipeline name:

~~~
nf-core launch <pipeline_name>
~~~
{: .language-bash}

## Using nf-core pipelines offline

Many of the techniques and resources described above require an active internet connection at run time - pipeline files, configuration profiles and software containers are all dynamically fetched when the pipeline is launched. This can be a problem for people using secure computing resources that do not have connections to the internet.

To help with this, the `nf-core download` command automates the fetching of required files for running nf-core pipelines offline. The command can download a specific release of a pipeline with `-r`/`--release` and fetch the singularity container if `--singularity` is passed (this needs Singularity to be installed). All files are saved to a single directory, ready to be transferred to the cluster where the pipeline will be executed.

> ## Exercise 3 (using pipelines)
>
>  Print the command-line usage instructions for the nf-core/rnaseq pipeline
>  In a new directory, run the nf-core/rnaseq pipeline with the provided test data
>  Try launching the RNA pipeline using the nf-core launch command
>  Download the nf-core/rnaseq pipeline for offline use using the nf-core download command
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}




{% include links.md %}
