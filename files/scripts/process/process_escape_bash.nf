//process_escape_bash.nf


process NUM_IDS {

  script:
  """
  #set bash variable NUMIDS
  NUMIDS=`zgrep -c '^>' $params.transcriptome`

  echo 'Number of sequences'
  printf "%'d\n" \$NUMIDS
  """
}

params.transcriptome = "${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

workflow {
  NUM_IDS()
}
