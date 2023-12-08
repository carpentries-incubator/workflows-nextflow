nextflow.enable.dsl = 2

params.message = 'hello'

workflow {
    PRINT_MESSAGE(params.message)
}

process PRINT_MESSAGE {
    echo true

    input:
    val my_message

    script:
    """
    echo $my_message
    """
}
