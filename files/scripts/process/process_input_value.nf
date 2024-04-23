//process_input_value.nf
nextflow.enable.dsl=2

process PRINTCHR {

  input:
  val chr

  script:
  """
  echo processing chromosome $chr
  """
}

chr_ch = Channel.of( 'A' .. 'P' )

workflow {

  PRINTCHR(chr_ch)
}
