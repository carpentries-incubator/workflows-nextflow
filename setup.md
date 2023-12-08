---
title: Setup
permalink: /setup/
---

# Running the lessons on your local machine

There are three items that you need to download:

1. The training material.
2. The training dataset.
3. The workshop scripts.

## Training material

Download the training material copy \& pasting the following command in the terminal:

```bash
$ git clone https://github.com/ggrimes/nf-training
$ cd nf-training
```

## Training software

The simplest way to install the software for this course is using conda.

To install conda see [here](https://carpentries-incubator.github.io/introduction-to-conda-for-data-scientists/setup/).

To create the training environment run:

```bash
$ conda env create -f environment.yml
```

Then activate the environment by running

```bash
$ conda activate nf-training
```

## Training scripts

To aid in the delivery of the lesson, the scripts mentioned in each episode, can be found in the respective episode folders in the github repository.

[https://github.com/carpentries-incubator/workflows-nextflow/tree/gh-pages/files/scripts](https://github.com/carpentries-incubator/workflows-nextflow/tree/gh-pages/files/scripts)

### Data

Inside the `nf-training` folder download the workshop dataset from Figshare, [https://figshare.com/articles/dataset/RNA-seq\_training\_dataset/14822481](https://figshare.com/articles/dataset/RNA-seq_training_dataset/14822481)

```bash
$ wget --content-disposition https://ndownloader.figshare.com/files/28531743
```

Unpack gzipped tar file:

```bash
$ tar -xvf  data.tar.gz
$ rm data.tar.gz
```

## Atom text editor setup

Any text editor can be used to write Nextflow scripts. A recommended text editor is [Atom](https://atom.io/).

Go to [https://atom.io](https://atom.io) and you should see a download button. The button or buttons should be specific to your platform and the download package should be  installable.

### MacOS

Atom follows the standard Mac zip installation process. You can either press the download button from the [https://atom.io](https://atom.io) site or you can go to the Atom releases page to download the atom-mac.zip file explicitly. Once you have that file, you can click on it to extract the application and then drag the new Atom application into your "Applications" folder.

### Nextflow language support in Atom

You can add Nextflow language support in Atom by clicking the [install](atom://settings-view/show-package?package=language-nextflow) button on the  atom package site [https://atom.io/packages/language-nextflow](https://atom.io/packages/language-nextflow).

### Atom terminal package

[https://atom.io/packages/atom-ide-terminal](https://atom.io/packages/atom-ide-terminal)

You can enable a terminal window within Atom by installing the atom-ide-terminal package from [https://atom.io/packages/atom-ide-terminal](https://atom.io/packages/atom-ide-terminal). Click the [install](atom://settings-view/show-package?package=atom-ide-terminal) button and follow the instructions.

Once installed enable the terminal window by selecting the Packages menu -> Atom IDE terminal -> Toggle menu item.

## Nextflow install without conda

Nextflow can be used on any [POSIX](https://en.wikipedia.org/wiki/POSIX) compatible system (Linux, OS X, etc). It requires Bash and Java 8 (or later, up to 12) to be installed.

Windows systems may be supported using a POSIX compatibility layer like Cygwin (unverified) or, alternatively, installing it into a Linux VM using virtualization software like VirtualBox or VMware.

## Nextflow installation

Install the latest version of Nextflow copy \& pasting the following snippet in a terminal window:

```
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
export NXF_VER=20.10.0
curl get.nextflow.io | bash
```

## Add Nextflow binary to your user's PATH:

```bash
$ mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

Check the correct installation running the following command:

```bash
$ nextflow info
```

## nf-core/tools installation without conda

### Pip

```bash
pip install nf-core
```




