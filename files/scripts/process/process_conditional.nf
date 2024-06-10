//process_conditional.nf


params.method = 'ids'
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"


process COUNT {
  script:
  if( params.method == 'ids' ) {
    """
    echo Number of sequences in transciptome
    zgrep -c "^>" $params.transcriptome
    """
  }  
  else if( params.method == 'bases' ) {
    """
    echo Number of bases in transciptome
    zgrep -v "^>" $params.transcriptome|grep -o "."|wc -l
    """
  }  
  else {
    """
    echo Unknown method $params.method
    """
  }  
}

workflow {
  COUNT()
}
