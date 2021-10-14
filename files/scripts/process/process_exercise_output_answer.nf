nextflow.enable.dsl = 2

process INDEX {

    input:
    path transcriptome
    each kmer

    output:
    path "index_${kmer}"

    script:
    """
    salmon index -t $transcriptome -i index_$kmer -k $kmer
    """
}

workflow {
    transcriptome_ch = Channel.fromPath( 'data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz' )
    kmer_ch = Channel.of( 21, 27, 31 )
    INDEX( transcriptome_ch, kmer_ch )
    INDEX.out.view()
}
