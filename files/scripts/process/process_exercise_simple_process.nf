nextflow.enable.dsl = 2

process SALMON_VERSION {

    script:
    """
    salmon --version
    """
}

workflow {
    SALMON_VERSION()
}
