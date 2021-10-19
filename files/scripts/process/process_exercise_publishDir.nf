nextflow.enable.dsl = 2

params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

process INDEX {

    //add publishDir directive here

    input:
    path transcriptome

    output:
    path "index"

    script:
    """
    salmon index -t $transcriptome -i index
    """
}


workflow{
    transcriptome_ch = Channel.fromPath( params.transcriptome, checkIfExists: true )
    INDEX( transcriptome_ch )
}
