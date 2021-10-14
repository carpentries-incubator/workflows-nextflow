bam_ch = Channel.fromPath("data/yeast/bams/*.bam")

process index {

    input:
    path bam from bam_ch

    output:
    path "${bam}*" into index_out_ch

    script:
    """
    samtools index $bam
    """
}
/*
*The flatMap operator applies a function to every item emitted by a channel, and returns the items so obtained as a new channel
*/
index_out_ch
    .flatMap()
    .view({ "File: ${it.name}" })
