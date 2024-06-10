//process_exercise_repeat_answer.nf


process COMBINE {
  input:
  path transcriptome
  each chr
 
  script:
  """
  printf "Number of sequences for chromosome $chr: "
  zgrep -c "^>Y${chr}" ${transcriptome}
  """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz', checkIfExists: true)
chr_ch = channel.of('A'..'P')

workflow {
  COMBINE(transcriptome_ch, chr_ch)
}
