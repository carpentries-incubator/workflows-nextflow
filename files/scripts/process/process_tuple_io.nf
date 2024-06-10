//process_tuple_io.nf


process COMBINE_FQ {
  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.fq.gz")

  script:
  """
  cat $reads > ${sample_id}.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  COMBINE_FQ(reads_ch)
  COMBINE_FQ.out.view()
}
