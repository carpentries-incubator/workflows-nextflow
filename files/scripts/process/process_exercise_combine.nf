// process_exercise_combine_answer.nf


process COMBINE {
 input:


 script:
 """
 zgrep -c ">Y${chr}" ${transcriptome}
 """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz', checkIfExists: true)
chr_ch = channel.of("A")

workflow {

}
