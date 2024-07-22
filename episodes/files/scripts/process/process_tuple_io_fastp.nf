//process_tuple_io_fastp.nf


process FASTP {
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("*FP*.fq.gz")
  
  script:
  """
  fastp \
   -i ${reads[0]} \
   -I ${reads[1]} \
   -o ${sample_id}_FP_R1.fq.gz \
   -O ${sample_id}_FP_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  FASTP(reads_ch)
  FASTP.out.view()
}
