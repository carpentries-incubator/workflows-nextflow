nextflow.enable.dsl = 2

process METHOD {

    input:
    val x

    output:
    path 'method.txt'

    """
    echo $x > method.txt
    """
}

workflow {
    methods_ch = Channel.of( 'salmon', 'kallisto' )
    METHOD( methods_ch )
    // use the view operator to display contents of the channel
    METHOD.out.view( { "Received: $it" } )
}
