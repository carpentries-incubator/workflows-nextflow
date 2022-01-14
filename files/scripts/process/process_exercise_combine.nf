nextflow.enable.dsl = 2

process COMBINE {

    input:
    <qualifier> <name>
    <qualifier> <name>

    script:
    """
    salmon index -t $transcriptome -i index -k $kmer
    """
}

workflow {

    transcriptome_ch = Channel.fromPath( 'data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz', checkIfExists: true )
    kmer_ch = Channel.of( 21 )
    COMBINE( transcriptome_ch, kmer_ch )
}
