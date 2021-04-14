---
layout: page
title: "Setup"
permalink: /setup/
root: ..
---

# Requirements

Nextflow can be used on any [POSIX](https://en.wikipedia.org/wiki/POSIX) compatible system (Linux, OS X, etc). It requires Bash and Java 8 (or later, up to 12) to be installed.

Windows systems may be supported using a POSIX compatibility layer like Cygwin (unverified) or, alternatively, installing it into a Linux VM using virtualization software like VirtualBox or VMware.



## Nextflow installation

Install the latest version of Nextflow copy & pasting the following snippet in a terminal window:

~~~
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
export NXF_VER=20.10.0
curl get.nextflow.io | bash
~~~



## Add Nextflow binary to your user's PATH:
~~~
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
~~~
{: .language-bash}

Check the correct installation running the following command:

~~~
nextflow info
~~~
{: .language-bash}


You can also install Nextflow using Bioconda:

~~~
conda install -c bioconda nextflow
~~~

Check the correct installation running the following command:

~~~
nextflow info
~~~
{: .language-bash}

## Training material

Download the training material copy & pasting the following command in the terminal:

~~~
wget -q -O- https://s3-eu-west-1.amazonaws.com/seqeralabs.com/public/nf-training.tar.gz | tar xvz
~~~
{: .language-bash}


Or in alternative if you have wget instead of curl:
~~~
curl https://s3-eu-west-1.amazonaws.com/seqeralabs.com/public/nf-training.tar.gz | tar xvz
~~~
{: .language-bash}

## Atom text editor setup

Any text editor can be used to write Nextflow scripts. A recommended text editor is [Atom](https://atom.io/).

Go to https://atom.io and you should see a download button. The button or buttons should be specific to your platform and the download package should be  installable.

### MacOS

Atom follows the standard Mac zip installation process. You can either press the download button from the https://atom.io site or you can go to the Atom releases page to download the atom-mac.zip file explicitly. Once you have that file, you can click on it to extract the application and then drag the new Atom application into your "Applications" folder.

### Nextflow language support in Atom

You can add Nextflow language support in Atom by clicking the [install](atom://settings-view/show-package?package=language-nextflow) button on the  atom package site https://atom.io/packages/language-nextflow.

### Atom terminal package

https://atom.io/packages/atom-ide-terminal

You can enable a terminal window within Atom by installing the atom-ide-terminal package from https://atom.io/packages/atom-ide-terminal. Click the [install](atom://settings-view/show-package?package=atom-ide-terminal) button and follow the instructions.

Once installed enable the terminal window by selecting the Packages menu -> Atom IDE terminal -> Toggle menu item.

## nf-core/tools installation

### Bioconda

You can install nf-core/tools using conda using the bioconda channel.


First, install conda and configure the channels to use bioconda (see the bioconda documentation).
To install conda see [here](https://carpentries-incubator.github.io/introduction-to-conda-for-data-scientists/setup/).
Then, just run the conda installation command:

~~~
conda install bioconda::nf-core=1.13
~~~
{: .language-bash}

Alternatively, you can create a new environment with both nf-core/tools and nextflow:

~~~
pip install nf-core
~~~
{: .language-bash}


# Optional Requirements

## Pipeline software

An analysis pipeline chains the execution of multiple tools together. Historically, all tools would have to be manually installed, often a source of great frustration and a key step where reproducibility between analyses is lost. nf-core pipelines utilise the built-in support for software packaging that Nextflow offers: all can work with Docker and Singularity, and most pipelines also have support for Conda.

* [Docker](https://docs.docker.com/install/), 1.10.x (or later):  
  * Typically used locally / on single-user servers and the cloud.
Analysis runs in a container, which behaves like an isolated operating system
Previously required system root access, though a "rootless mode" is available since late 2019
* [Singularity](https://www.sylabs.io/), 2.5.x (or later, optional):
  * Often used as an alternative to Docker on multi-user systems such as HPC systems.
Also runs containers and can create these from Docker images
Does not need root access or any daemon processes - images built from files
* [Conda](https://conda.io/),Conda 4.5 (or later, optional):
  * Packaging system that manages environments instead of running analysis in containers.
Poorer reproducibility than Docker / Singularity
There can be changes in low-level package dependencies over time
The software still runs in your native operating system environment and so core system functions can differ

## Visualisation software

* [Graphviz](http://www.graphviz.org/) (optional): Graphviz is open source graph visualization software. Graph visualization is a way of representing structural information as diagrams of abstract graphs and network


{% include links.md %}
