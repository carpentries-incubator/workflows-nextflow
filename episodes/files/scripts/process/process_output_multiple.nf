//process_output_multiple.nf


params.transcriptome="${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

process SPLIT_FASTA {
  input:
  path transcriptome

  output:
  path "*"

  script:
  """
  zgrep  '^>' $transcriptome > sequence_ids.txt
  zgrep -v '^>' $transcriptome > sequence.txt
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
transcriptome_ch = channel.fromPath(params.transcriptome)

workflow {
  SPLIT_FASTA(transcriptome_ch)
  // use the view operator to display contents of the channel
  SPLIT_FASTA.out.view()
}
