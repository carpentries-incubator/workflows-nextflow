nextflow.enable.dsl = 2

include { INDEX; QUANT } from './modules/rnaseq-tasks.nf'

workflow {

    reads = Channel.fromFilePairs( 'data/yeast/reads/*_{1,2}.fq.gz' )
    transcriptome_ch = Channel.fromPath( 'data/yeast/transcriptome/*.fa.gz' )
    INDEX( transcriptome_ch )
    QUANT( index.out,reads )
}
