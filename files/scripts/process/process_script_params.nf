//process_script_params.nf


params.chr = "A"

process CHR_COUNT {

  script:
  """
  printf  'Number of sequences for chromosome '${params.chr}':'
  zgrep  -c '^>Y'${params.chr} ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
  """
}

workflow {
  CHR_COUNT()
}
