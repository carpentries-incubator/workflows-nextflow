nextflow.enable.dsl = 2

process FASTQC {

    //add tag directive
    //add cpu directive

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

workflow {
    read_pairs_ch = Channel.fromFilePairs( 'data/yeast/reads/ref*_{1,2}.fq.gz' )
    FASTQC( read_pairs_ch )
    FASTQC.out.view()
}
