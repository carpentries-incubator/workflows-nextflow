---
title: Background and Software
teaching: 10
exercises: 0
---

::::::::::::::::::::::::::::::::::::::: objectives

- Why study *Saccharomyces cerevisiae?*?
- Understand the data set.
- What is RNA-Seq?

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What data are we using?
- Why is this experiment important?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Background

We are going to use a  stress response  dataset from a population of *Saccharomyces cerevisiae*.

- **What is *Saccharomyces cerevisiae?*?**
  - *Saccharomyces cerevisiae* is a species of yeast used extensively in baking, brewing, and scientific research. It is a single-celled fungus that has become one of the most comprehensively studied eukaryotic model organisms in molecular and cell biology.
 
![](fig/172px-EscherichiaColi_NIAID.jpg){alt='Wikimedia'}

<!-- https://species.wikimedia.org/wiki/Escherichia_coli#/media/File:EscherichiaColi_NIAID.jpg -->

## The data

- The data we are going to use is part of a stress reponse experiment led by [Petri-Jaan Lahtvee](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4966989/).

- The experiment was designed to measure the effect of the three most commonly encountered industrial stresses for yeast (ethanol, salt, and temperature) to identify the mechanisms of general and stress-specific responses under chemostat conditions in which specific growth rate–dependent changes are eliminated.

- The experiment found varying phenotypic responses, with all stresses leading to increased cellular maintenance costs. This prompted further RNA sequencing analysis for deeper insights.

- RNA-Seq data was generated using  Illumina TruSeq sample preparation kit, with poly-A selection,  sequenced on  Illumina HiSeq 2500, generating paired end reads (2× 100 base pairs).

- Read data was downsampled to 1% for the lesson.

### Software 

We will be using four bioinformatic software tools


| Software           | Description                                     | 
| ---------------- | ----------------------------------------------- |
| FastQC           | A quality control tool for high throughput sequence data.                                   | 
| salmon           | A tool for quantifying the expression of transcripts using RNA-seq data               | 
| samtools         | A suite of programs for interacting with high-throughput sequencing data                  | 
| multiqc          | A suite of programs for interacting with high-throughput sequencing data | 


:::::::::::::::::::::::::::::::::::::::: keypoints

- It is important to record and understand your experiment's metadata.

::::::::::::::::::::::::::::::::::::::::::::::::::
