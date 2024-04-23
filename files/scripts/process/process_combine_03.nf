nextflow.enable.dsl = 2

process COMBINE {

    echo true

    input:
    val x
    val y

    script:
    """
    echo $x and $y
    """
}


workflow {
    ch_num = Channel.value( 1 )
    ch_letters = Channel.of( 'a', 'b', 'c' )
    COMBINE( ch_num, ch_letters )
}
