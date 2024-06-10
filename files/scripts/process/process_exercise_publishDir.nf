//process_exercise_publishDir.nf


params.reads= "data/yeast/reads/ref{1,2,3}*{1,2}.fq.gz"

process MERGE_REPS {
  
  input:
  tuple val(sample_id), path(reads)
  output:
  path("*fq.gz")

  script:
  """
  cat *1.fq.gz > ${sample_id}.merged.R1.fq.gz
  cat *2.fq.gz > ${sample_id}.merged.R2.fq.gz
  """
}
reads_ch = Channel.fromFilePairs(params.reads,checkIfExists:true,size:6)

workflow {
  MERGE_REPS(reads_ch)
}
