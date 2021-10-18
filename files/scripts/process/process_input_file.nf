nextflow.enable.dsl = 2

process NUMLINES {

    input:
    path read

    script:
    """
    printf '${read} '
    gunzip -c ${read} | wc -l
    """
}


workflow {
    reads_ch = Channel.fromPath( 'data/yeast/reads/ref*.fq.gz' )
    NUMLINES( reads_ch )
}
