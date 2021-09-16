---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---


[Nextflow](https://www.nextflow.io/) is workflow management software which
enables the writing of scalable and reproducible scientific workflows. It
can integrate various software package and environment management systems
such as Docker, Singularity, and Conda. It allows for existing pipelines
written in common scripting languages, such as R and Python, to
be seamlessly coupled together. It implements a Domain Specific Language
(DSL) that simplifies the implementation and running of workflows on
cloud or high-performance computing (HPC) infrastructures.

This lesson also introduces [nf-core](https://nf-co.re/): a
community-driven platform, which provide peer reviewed
best practice analysis pipelines written in Nextflow.

This lesson motivates the use of Nextflow and nf-core as development tools
for building and sharing reproducible data science workflows.

## lesson objectives

1. The learner will understand the fundamental components of a Nextflow
script, including channels, processes and operators.
1. The learner will write a multi-step workflow script to align, quantify, and perform QC on an RNA-Seq data in Nextflow DSL.
1. The learner will be able to write a Nextflow configuration file to alter the computational resources allocated to a process.
1. The learner will use nf-core to run a community curated pipeline, on an RNA-Seq dataset.

{% comment %} This is a comment in Liquid {% endcomment %}


> ## Prerequisites
>
> This is an intermediate lesson and assumes familiarity with the core materials covered in the
> [Software Carpentry Lessons] [swc-lessons]. In particular learners need to be familiar with
> material covered in [The Unix Shell](http://swcarpentry.github.io/shell-novice),
> and either
> [Plotting and Programming in Python](http://swcarpentry.github.io/python-novice-gapminder) or
> [R for Reproducible Scientific Analysis](http://swcarpentry.github.io/r-novice-gapminder).
{: .prereq}


{: .prereq}

{% include links.md %}
