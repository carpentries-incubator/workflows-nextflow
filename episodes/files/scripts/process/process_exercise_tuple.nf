//process_exercise_tuple.nf


process COMBINE_REPS {
  input:
  tuple ___(sample_id), ___(reads)

  output:
  tuple ___(sample_id), ___("*.fq.gz")

  script:
  """
  cat *_1.fq.gz > ${sample_id}_R1.fq.gz
  cat *_2.fq.gz > ${sample_id}_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref{1,2,3}*.fq.gz',size:-1)

workflow{
  COMBINE_REPS(reads_ch)
  COMBINE_REPS.out.view()
}
