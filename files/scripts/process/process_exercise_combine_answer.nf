// process_exercise_combine_answer.nf
nextflow.enable.dsl=2
process COMBINE {
 input:
 path transcriptome
 val chr

 script:
 """
 zgrep -c ">Y{chr}" ${transcriptome}
 """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz', checkIfExists: true)
kmer_ch = channel.of(21)

workflow {
  COMBINE(transcriptome_ch, chr_ch)
}
