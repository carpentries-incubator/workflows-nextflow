#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "data/yeast/reads/ref1_1.fq.gz"

process NUM_LINES {

    input:
    path read

    output:
    stdout

    script:
    """
    printf '${read} '
    gunzip -c ${read} | wc -l
    fastp -i ${read}  -o out.fq 2>&1
    """
}

workflow {

    input_ch = Channel.fromPath(params.input)
    NUM_LINES(input_ch).out.view()

}

