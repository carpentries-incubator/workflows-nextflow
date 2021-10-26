nextflow.enable.dsl = 2

process FASTQC {

    input:
    tuple val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "fastqc_out" )

    script:
    """
    mkdir fastqc_out
    fastqc $reads -o fastqc_out -t 1
    """
}

workflow{
    reads_ch = Channel.fromFilePairs( 'data/yeast/reads/ref*_{1,2}.fq.gz' )
    FASTQC( reads_ch )
    FASTQC.out.view()
}
