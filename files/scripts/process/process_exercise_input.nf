//process_exercise_input.nf


params.chr = "A"
params.transcriptome = "${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
process CHR_COUNT {

 script:
 """
 printf  'Number of sequences for chromosome '${params.chr}':'
 zgrep  -c '^>Y'${params.chr} ${params.transcriptome}
 """
}

workflow {
 CHR_COUNT()
}
