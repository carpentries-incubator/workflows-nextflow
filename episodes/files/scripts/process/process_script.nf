//process_script.nf


chr = "A"

process CHR_COUNT {

  script:
  """
  printf "Number of sequences for chromosome ${chr} :"
  zgrep -c '>Y'${chr} ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
  """
}

workflow {
  CHR_COUNT()
}
