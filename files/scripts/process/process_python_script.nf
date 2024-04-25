//process_python_script.nf
nextflow.enable.dsl=2

process PROCESS_READS {

  script:
  """
  process_reads.py ${projectDir}/data/yeast/reads/ref1_1.fq.gz
  """
}

workflow {
  PROCESS_READS()
}