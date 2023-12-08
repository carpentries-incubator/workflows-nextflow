nextflow.enable.dsl = 2

process INDEX {

    input:
    path transcriptome

    output:
    path 'index', emit: salmon_index

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
    transcriptome_ch = Channel.fromPath( 'data/yeast/transcriptome/*.fa.gz' )
    read_pairs_ch = Channel.fromFilePairs( 'data/yeast/reads/*_{1,2}.fq.gz' )
    INDEX( transcriptome_ch )
    QUANT( INDEX.out.salmon_index, read_pairs_ch ).view()
}
