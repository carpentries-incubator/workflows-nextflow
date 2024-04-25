//process_script_params.nf
nextflow.enable.dsl=2

params.base='A'

process COUNT_BASES {
 
script:
 """
 zgrep -v  '^>'   ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz|grep -o ${params.base}|wc -l   
 """
}

workflow {
  COUNT_BASES()
}
