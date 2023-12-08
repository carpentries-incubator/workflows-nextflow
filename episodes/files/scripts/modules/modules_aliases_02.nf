nextflow.enable.dsl = 2

include { INDEX; INDEX as SALMON_INDEX } from './modules/rnaseq-tasks.nf'

workflow {

    transcriptome_ch = Channel.fromPath( '/data/yeast/transcriptome/*.fa.gz' )
    INDEX( transcriptome )
    SALMON_INDEX( transcriptome )

}
