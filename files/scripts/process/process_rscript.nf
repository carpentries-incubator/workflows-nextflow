nextflow.enable.dsl = 2

process RSCRIPT {

    script:
    """
    #!/usr/bin/env Rscript
    library("ShortRead")
    countFastq(dirPath="data/yeast/reads/ref1_1.fq.gz")
    """
}

workflow {
    RSCRIPT()
}
