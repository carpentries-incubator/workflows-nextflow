nextflow.enable.dsl = 2

/*
* pipeline input parameters
*/
params.reads = "data/yeast/reads/ref1_{1,2}.fq.gz"
params.transcriptome = "data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
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
 * define the `INDEX` process that create a binary index
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
    tag "quantification on $pair_id"
    publishDir params.outdir, mode:'copy'

    input:
    each index
    tuple val(pair_id), path(reads)

    output:
    path(pair_id)

    script:
    """
    salmon quant \\
        --threads $task.cpus \\
        --libType=U -i $index \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o $pair_id
    """
}

/*
 * Run fastQC to check quality of reads files
 */
process FASTQC {

    tag "FASTQC on $sample_id"
    cpus 1

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
    """
}

workflow {
   read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )
   transcriptome_ch = Channel.fromPath( params.transcriptome, checkIfExists:true )

   index_ch = INDEX( transcriptome_ch )
   quant_ch = QUANT( index_ch, read_pairs_ch )
}
