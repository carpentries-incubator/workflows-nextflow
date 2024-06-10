//workflow_01.nf



 process FASTQC {
    input:
      tuple(val(sample_id), path(reads))
    output:
      path "fastqc_${sample_id}_logs"
    script:
      """
      mkdir fastqc_${sample_id}_logs
      fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
      """
}

process MULTIQC {
    publishDir "results/mqc"
    input:
      path transcriptome
    output:
      path "*"
    script:
      """
      multiqc .
      """
}

workflow {
    read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz',checkIfExists: true)

    //index process takes 1 input channel as a argument
    //assign process output to Nextflow variable fastqc_obj
    fastqc_obj = FASTQC(read_pairs_ch)

    //quant channel takes 1 input channel as an argument
    //We use the collect operator to gather multiple channel items into a single item
    MULTIQC(fastqc_obj.collect()).view()
}
