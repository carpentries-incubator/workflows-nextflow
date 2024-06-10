

params.chr = "A"
params.transcriptome = "${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"> >

process CHR_COUNT {
 input:
 path transcriptome

 script:
 """
 printf  'Number of sequences for chromosome '${params.chr}':'
 zgrep  -c '^>Y'${params.chr} ${transcriptome}
 """
}

transcriptome_ch = channel.fromPath(params.transcriptome)

workflow {
 CHR_COUNT(transcriptome_ch)
}
