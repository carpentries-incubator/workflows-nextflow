nextflow.enable.dsl = 2

process COMBINE {

    input:
    val x
    each y

    script:
    """
    echo $x and $y
    """
}


workflow {
    ch_num = Channel.of( 1, 2 )
    ch_letters = Channel.of( 'a', 'b', 'c', 'd' )
    COMBINE( ch_num, ch_letters )
}
