nextflow.enable.dsl = 2

params.reads = 'data/yeast/reads/*_{1,2}.fq.gz'

process FASTQC {
    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs/*.zip"

    script:
    //flagstat simple stats on bam file
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
    """
}

process PARSEZIP {

    publishDir "results/fqpass", mode:"copy"

    input:
    path flagstats

    output:
    path 'pass_basic.txt'

    script:
    """
    for zip in *.zip; do
        zipgrep 'Basic Statistics' \$zip \\
        | grep 'summary.txt'
    done > pass_basic.txt
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
    PARSEZIP( FASTQC( read_pairs_ch ).collect() )
}
