nextflow.enable.dsl = 2

params.aligner = 'kallisto'
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.kmer = 31

process INDEX {

    script:
    if( params.aligner == 'kallisto' ) {
        """
        echo indexed using kallisto
        kallisto index -i index  -k $params.kmer $params.transcriptome
        """
    } else if( params.aligner == 'salmon' ) {
        """
        echo indexed using salmon
        salmon index -t $params.transcriptome -i index --kmer $params.kmer
        """
    } else {
        """
        echo Unknown aligner $params.aligner"
        """
    }
}

workflow {
    INDEX()
}
