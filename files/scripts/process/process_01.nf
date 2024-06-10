

process NUMSEQ {
  script:
  "zgrep -c '^>' ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
}

workflow {
  //process is called like a function in the workflow block
  NUMSEQ()
}
