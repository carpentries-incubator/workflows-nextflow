nextflow.enable.dsl = 2

process PRINTCHR {

    input:
    val chr

    script:
    """
    echo processing chromosome $chr
    """
}


workflow {
    chr_ch = Channel.of( 1..22, 'X', 'Y' )
    PRINTCHR( chr_ch )
}
