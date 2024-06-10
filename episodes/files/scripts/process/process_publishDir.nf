//process_publishDir.nf


process COMBINE_READS {
  publishDir "results/merged_reads"

  input:
  tuple val(sample_id), path(reads)

  output:
  path("${sample_id}.merged.fq.gz")

  script:
  """
  cat ${reads} > ${sample_id}.merged.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')


workflow {
  COMBINE_READS(reads_ch)
}
