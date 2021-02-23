---
title: Setup
---
# Requirements

Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash and Java 8 (or later, up to 12) to be installed.

Optional requirements for this workshop:

* Docker engine 1.10.x (or later)

* Singularity 2.5.x (or later, optional)

* Conda 4.5 (or later, optional)Nextflow installation

* Graphviz (optional)

* AWS CLI (optional)

* AWS Batch computing environment properly configured (optional)

## Nextflow installation

Install the latest version of Nextflow copy & pasting the following snippet in a terminal window:

~~~
curl get.nextflow.io | bash
mv nextflow ~/bin
~~~
{: .language-bash}

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
wget -q -O- https://s3-eu-west-1.amazonaws.com/seqeralabs.com/public/nf-training.tar.gz | tar xvz
~~~
{: .language-bash}


{% include links.md %}
