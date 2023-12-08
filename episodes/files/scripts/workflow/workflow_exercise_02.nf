nextflow.enable.dsl = 2

params.bam = 'data/yeast/bams/*.bam'

process FLAGSTAT {

    input:
    path bam

    output:
    path "${bam}.flagstats.txt"

    script:
    //flagstat simple stats on bam file
    """
    samtools flagstat ${bam} > ${bam}.flagstats.txt
    """
}

process MERGEFLAGSTAT {

    publishDir "results/flagstats", mode:"copy"

    input:
    path flagstats

    output:
    path 'flagstats.txt'

    script:
    """
    cat ${flagstats} | grep 'mapped (' > flagstats.txt
    """
}

workflow {
    bam_ch = Channel.fromPath( params.bams )
}
