nextflow.enable.dsl = 2

kmer = 31

process INDEX {

    script:
    """
    salmon index \\
        -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \\
        -i index \\
        --kmer $kmer
    echo "kmer size is $kmer"
    """
}

workflow {
    INDEX()
}
