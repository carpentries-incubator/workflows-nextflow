---
title: "Instructor Notes"
---


## When to use a workflow system?

You may be asked by the learners when you would need to use a workflow system. The paper [Workflow systems turn raw data into scientific knowledge](https://pubmed.ncbi.nlm.nih.gov/31477884/) has a view on this.

"So, do you need a workflow system? Not every task requires one, and there is a learning curve. Scripting usually suffices for one-off tasks and when working out the pipeline itself. The tipping point, most agree, comes when you need to run the same workflow over and over again, or if the data are likely to be published."


## Nextflow echo option

It can be useful to use the `-process.echo` option to echo output from a process when you run  Nextflow . 

~~~
nextflow run main.nf -process.echo
~~~


## Nextflow video

The nextflow lesson material has been adapted from the sequera-labs training course [link here](https://seqera.io/training) there is video of this material being presented [here]( https://youtu.be/8_i8Tn335X0)

## Editors

The [Atom editor](https://atom.io/) has [syntax highlighting](https://atom.io/packages/language-nextflow) for Nextflow.

## Help

Both Nextflow and nf-core have community forums for asking questions and getting help.

* Nextflow [Gitter](https://gitter.im/nextflow-io/nextflow)


## References

* [Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319.](https://www.nature.com/articles/nbt.3820) doi: 10.1038/nbt.3820). doi: 10.1038/nbt.3820.
* [Ewels PA, Peltzer A, Fillinger S, et al. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020;38(3):276-278.] (https://www.nature.com/articles/s41587-020-0439-x)  DOI: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
* [Perkel JM. Workflow systems turn raw data into scientific knowledge. Nature. 2019;573(7772):149-150.](https://www.nature.com/articles/d41586-019-02619-z) DOI: 10.1038/d41586-019-02619-z


{% include links.md %}
