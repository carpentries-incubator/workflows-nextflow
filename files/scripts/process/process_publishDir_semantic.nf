nextflow.enable.dsl = 2

process QUANT {

    publishDir "results/bams", pattern: "*.bam", mode: "copy"
    publishDir "results/quant", pattern: "${sample_id}_salmon_output", mode: "copy"

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}.bam")
    path "${sample_id}_salmon_output"

    script:
    """
    salmon quant -i $index \\
        -l A \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${sample_id}_salmon_output \\
        --writeMappings \\
        | samtools sort \\
        | samtools view -bS -o ${sample_id}.bam
    """
}

workflow {
    reads_ch = Channel.fromFilePairs( 'data/yeast/reads/ref1_{1,2}.fq.gz' )
    index_ch = Channel.fromPath( 'data/yeast/salmon_index' )
    QUANT( reads_ch, index_ch )
}
