nextflow.enable.dsl = 2

params.kmer = 31

process INDEX {

    shell:
    '''
    #set bash variable KMERSIZE
    KMERSIZE=!{params.kmer}
    salmon index -t !{projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer ${KMERSIZE}
    echo "kmer size is  !{params.kmer}"
    '''
}

workflow {
    INDEX()
}
