nextflow.enable.dsl = 2

process INDEX {

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index -t $transcriptome -i index
    """
}

process QUANT {

    input:
    each path(index)
    tuple( val(pair_id), path(reads) )

    output:
    path pair_id

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

workflow {

    transcriptome_ch = Channel.fromPath( 'data/yeast/transcriptome/*.fa.gz', checkIfExists: true )
    read_pairs_ch = Channel.fromFilePairs( 'data/yeast/reads/*_{1,2}.fq.gz', checkIfExists: true )

    //index process takes 1 input channel as a parameter
    INDEX( transcriptome_ch )
    index_obj = INDEX.out

    //quant channel takes 2 input channels as parameters
    QUANT( index_obj, read_pairs_ch ).view()
}
