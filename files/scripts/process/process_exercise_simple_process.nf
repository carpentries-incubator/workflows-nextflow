

process COUNT_BASES {
   
  script:
  """
  zgrep -v '^>' ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz|tr -d '\n'|wc -m
  """
}

workflow {
  COUNT_BASES()
}
