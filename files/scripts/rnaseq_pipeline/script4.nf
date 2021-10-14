nextflow.enable.dsl = 2

/*
* pipeline input parameters
*/
params.reads = "$projectDir/data/yeast/reads/ref1_{1,2}.fq.gz"
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
 * define the `index` process that create a binary index
 * given the transcriptome file
 */
process INDEX {

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process QUANT {

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path(pair_id)

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

workflow {
    Channel.fromFilePairs( params.reads, checkIfExists:true )
        .set { read_pairs_ch }

    index_ch = INDEX( params.transcriptome )
    quant_ch = QUANT( index_ch, read_pairs_ch )
}
