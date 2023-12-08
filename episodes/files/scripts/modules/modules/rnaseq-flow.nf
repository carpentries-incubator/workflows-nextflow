/*
 * include requires tasks
 */
include { INDEX; QUANT; FASTQC; MULTIQC  } from './rnaseq-tasks.nf'

/*
 * define the data analysis workflow
 */
workflow RNASEQFLOW {
    // required inputs
    take:
    transcriptome
    read_files
    // workflow implementation
    main:
    INDEX( transcriptome )
    QUANT( INDEX.out, read_files )
    FASTQC( read_files )
    MULTIQC( QUANT.out.mix( FASTQC.out ).collect() )
}
