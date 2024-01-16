---
title: Instructor Notes
---

## When to use a workflow system?

You may be asked by the learners when you would need to use a workflow system. The paper [Workflow systems turn raw data into scientific knowledge](https://pubmed.ncbi.nlm.nih.gov/31477884/) has a view on this:

> So, do you need a workflow system? Not every task requires one, and there is a learning curve. Scripting usually suffices for one-off tasks and when working out the pipeline itself. The tipping point, most agree, comes when you need to run the same workflow over and over again, or if the data are likely to be published.

### Which workflow system to use?

The paper [Using prototyping to choose a bioinformatics workflow management system](https://doi.org/10.1371/journal.pcbi.1008622) discusses why they chose nextflow compared to other workflow systems. It also has a view about when to use workflow systems:

> We conclude that many [...] multistep data analysis workflows can be rewritten in a workflow management system, and we advocate prototyping as a low-cost (both time and effort) way of making an informed selection of software for use within a research project.

## Nextflow echo option

It can be useful to use the `-process.echo` option to echo output from a process when you run  Nextflow .

```
nextflow run main.nf -process.echo
```

## Nextflow video

The nextflow lesson material has been adapted from the sequera-labs training course [link here](https://seqera.io/training) there is video of this material being presented [here](https://youtu.be/8_i8Tn335X0)

## Code Editor

The preferred code editor is [Visual Studio Code](https://code.visualstudio.com).  Visual Studio Code has a [syntax highlighting](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) extension that adds the Nextflow language support.

## Getting Help

For questions about the main Nextflow tool, use the Nextflow Gitter chat community: [https://gitter.im/nextflow-io/nextflow](https://gitter.im/nextflow-io/nextflow)

Many nf-core resources are to be found on this website. If in doubt, get in touch with the nf-core community via [Slack](https://nf-co.re/join).

The nf-core Slack channel has a #carpentries-course channel for asking questions about the course.

## Nextflow REPL console

Nextflow has a console graphical interface. The Nextflow console is a REPL (read-eval-print loop) environment that allows one to quickly test part of a script or pieces of Nextflow code in an interactive manner.

It is a handy tool that allows one to evaluate fragments of Nextflow/Groovy code or fast prototype a complete pipeline script.

See the blog post [here](https://www.nextflow.io/blog/2015/introducing-nextflow-console.html) for more information.

```bash
nextflow console
```

## Next Steps

Nextflow produces a [blog post](https://www.nextflow.io/blog/2023/learn-nextflow-in-2023.html) each year with links to training materials and workshops.
This is a good website to show learners who want to learn more.



## References

- [Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319.](https://www.nature.com/articles/nbt.3820) doi: 10.1038/nbt.3820). doi: 10.1038/nbt.3820.
- [Ewels PA, Peltzer A, Fillinger S, et al. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020;38(3):276-278.](https://www.nature.com/articles/s41587-020-0439-x)  DOI: [10\.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
- [Perkel JM. Workflow systems turn raw data into scientific knowledge. Nature. 2019;573(7772):149-150.](https://www.nature.com/articles/d41586-019-02619-z) DOI: 10.1038/d41586-019-02619-z
- [Jackson M, Wallace EWJ, Kavoussanakis K, PLoS Computational Biology. Using prototyping to choose a bioinformatics workflow management system. 2021; 17 :e1008622.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008622) [doi: 10.1371/journal.pcbi.1008622](https://doi.org/10.1371/journal.pcbi.1008622)




