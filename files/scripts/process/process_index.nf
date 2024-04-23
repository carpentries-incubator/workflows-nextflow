nextflow.enable.dsl = 2

process INDEX {

    script:
    "salmon index -t ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i data/yeast/salmon_index --kmerLen 31"
}

workflow {
    //process is called like a function in the workflow block
    INDEX()
}
