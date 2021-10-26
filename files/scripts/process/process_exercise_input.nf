nextflow.enable.dsl = 2

process FASTQC {

    // add input channel

    script:
    """
    mkdir fastqc_out
    fastqc -o fastqc_out ${reads}
    ls -1 fastqc_out
    """
}


workflow {
    reads_ch = Channel.fromPath( 'data/yeast/reads/ref1*_{1,2}.fq.gz' )
    FASTQC( reads_ch )
}
