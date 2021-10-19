nextflow.enable.dsl = 2

process FASTQC {

    input:
    path read

    output:
    path "fqc_res/*"

    script:
    """
    mkdir fqc_res
    fastqc $read -o fqc_res
    """
}


workflow {
    read_ch = Channel.fromPath( "data/yeast/reads/ref1*.fq.gz" )
    FASTQC( read_ch )
    FASTQC.out.view()
}
