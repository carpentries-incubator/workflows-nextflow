nextflow.enable.dsl = 2

params.reads = "$projectDir/data/yeast/reads/ref1_{1,2}.fq.gz"

workflow {

    read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )
    FASTQC( read_pairs_ch )
    MULTIQC( fastqc.out.collect() )

}
