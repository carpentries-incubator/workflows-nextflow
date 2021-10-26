nextflow.enable.dsl = 2

include { FASTQC } from './modules/rnaseq-tasks.nf'

params.reads = "$projectDir/data/yeast/reads/ref1_{1,2}.fq.gz"

workflow {

    read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )
    FASTQC( read_pairs_ch )
}
