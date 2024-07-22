---
title: Setup
permalink: /setup/
---




# Running the lessons on your local machine

## Training directory

Each learner should setup a training folder e.g. `nf-training`

```bash
mkdir nf-training
cd nf-training
```


There are three items that you need to download:


1. The training software.
2. The training dataset.
3. The workshop scripts.



## Training software

A list of software with version required for this training is listed below:

|Software|Version|
|--------|-------|
|Nextflow|20.10.0|
|nf-core/tools|1.12.1|
|salmon|1.5|
|fastqc|0.11|
|multiqc|1.10|
|python|3.8|

### conda

The simplest way to install the software for this course is using conda.


To install conda see [here](https://carpentries-incubator.github.io/introduction-to-conda-for-data-scientists/setup/). 

An environment file is provided here [environment.yml](https://raw.githubusercontent.com/carpentries-incubator/workflows-nextflow/main/episodes/data/environment.yml)

```bash
wget
wget https://raw.githubusercontent.com/carpentries-incubator/workflows-nextflow/main/episodes/data/environment.yml

# or curl 
curl -L -o environment.yml https://raw.githubusercontent.com/carpentries-incubator/workflows-nextflow/main/episodes/data/environment.yml
```

To create the training environment run:

```bash
conda env create -n nf-training -f environment.yml
```

Then activate the environment by running

```bash
conda activate nf-training
```

## Training scripts


To aid in the delivery of the lesson, the scripts mentioned in each episode, can be found in the respective episode folders in the github repository.
[https://github.com/carpentries-incubator/workflows-nextflow/tree/main/episodes/files/scripts](https://github.com/carpentries-incubator/workflows-nextflow/tree/gh-pages/files/scripts)

To get the scripts associated with each episode you will need to download the scripts folder from the github repository.

Below is a series of commands to download and unpack scripts folder.

```bash
# get the gitrepo as a zip file
wget https://github.com//carpentries-incubator/workflows-nextflow/archive/main.zip

#or
curl -L -o main.zip https://github.com//carpentries-incubator/workflows-nextflow/archive/main.zip

# unzip the script file
unzip main.zip 'workflows-nextflow-main/episodes/files/scripts*' -d  .

# mv the scripts folder to the nf-training folder 
mv workflows-nextflow-main/episodes/files/scripts .

# remove the zip file and the git repo
rm -r workflows-nextflow-main main.zip
```

The nextflow scripts for each episode, can be found in the respective episode folders inside this the scripts folder.


### Data

Inside the `nf-training` folder download the workshop dataset from Figshare, [https://figshare.com/articles/dataset/RNA-seq\_training\_dataset/14822481](https://figshare.com/articles/dataset/RNA-seq_training_dataset/14822481)

```bash
wget --content-disposition https://ndownloader.figshare.com/files/28531743

# or curl
curl -L -o  data.tar.gz https://ndownloader.figshare.com/files/28531743
```

Unpack gzipped tar file:

```bash
tar -xvf  data.tar.gz
rm data.tar.gz
```

## Visual Studio Code editor setup

Any text editor can be used to write Nextflow scripts. A recommended  code editor is [Visual Studio Code](https://code.visualstudio.com/).

Go to [Visual Studio Code](https://code.visualstudio.com/) and you should see a download button. The button or buttons should be specific to your platform and the download package should be  installable.


### Nextflow language support in Visual Studio Code

You can add Nextflow language support in Visual Studio Code by clicking the [install](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) button on the Nextflow language extension.


## Nextflow install without conda

Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through WSL. It requires Bash 3.2 (or later) and Java 11 (or later, up to 22) to be installed

## Nextflow installation

Install the latest version of Nextflow copy \& pasting the following snippet in a terminal window:

```bash
# Make sure that Java v11 or later is installed:
java -version

# Install Nextflow
curl -s https://get.nextflow.io | bash
```

## Add Nextflow binary to your user's PATH:

```bash
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

Check the correct installation running the following command:

```bash
nextflow info
```

## nf-core/tools installation without conda

### Pip

```bash
pip install nf-core
```




