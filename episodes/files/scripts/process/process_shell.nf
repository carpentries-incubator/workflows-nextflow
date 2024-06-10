//process_shell.nf


process NUM_IDS {

  shell:
  //Shell script definition requires the use of single-quote ' delimited strings
  '''
  #set bash variable NUMIDS
  NUMIDS=`zgrep -c '^>' !{params.transcriptome}`

  echo 'Number of sequences'
  printf  $NUMIDS
  '''
}

params.transcriptome = "${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

workflow {
  NUM_IDS()
}
