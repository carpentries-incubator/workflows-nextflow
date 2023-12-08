#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*  Comments are uninterpreted text included with the script.
    They are useful for describing complex parts of the workflow
    or providing useful information such as workflow usage.

    Usage:
       nextflow run wc.nf --input <input_file>

    Multi-line comments start with a slash asterisk /* and finish with an asterisk slash. */
//  Single line comments start with a double slash // and finish on the same line

/*  Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator. */
params.input = "data/yeast/reads/ref1_1.fq.gz"

//  The default workflow
workflow {

    //  Input data is received through channels
    input_ch = Channel.fromPath(params.input)

    /*  The script to execute is called by it's process name,
        and input is provided between brackets. */
    NUM_LINES(input_ch)

    /*  Process output is accessed using the `out` channel.
        The channel operator view() is used to print
        process output to the terminal. */
    NUM_LINES.out.view()
}

/*  A Nextflow process block
    Process names are written, by convention, in uppercase.
    This convention is used to enhance workflow readability. */
process NUM_LINES {

    input:
    path read

    output:
    stdout

    script:
    /* Triple quote syntax """ (triple double quotes),
    or ''' (triple single quotes) are known as multi-line strings.
    The content of the string can be written on multiple lines
    preserving newlines, and removes the need to concatenate
    multiple short strings. Nextflow variables (variables beginning
    with `$`) are interpolated inside triple double-quoted strings.
    This then means bash variables (which also begin with $) should
    be escaped using `\` e.g., "\$HOME". */
    """
    printf '${read} '
    gunzip -c ${read} | wc -l
    """
}
