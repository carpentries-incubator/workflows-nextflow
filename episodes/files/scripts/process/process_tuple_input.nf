nextflow.enable.dsl = 2

process TUPLEINPUT{

    input:
    tuple val(sample_id), path(reads)

    script:
    """
    echo $sample_id
    echo $reads
    """
}

workflow {
    reads_ch = Channel.fromFilePairs( 'data/yeast/reads/ref1_{1,2}.fq.gz' )
    TUPLEINPUT( reads_ch )
}
