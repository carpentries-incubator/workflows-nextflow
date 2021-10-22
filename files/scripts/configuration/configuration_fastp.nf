// configure_fastp.nf
nextflow.enable.dsl = 2

params.input = "data/yeast/reads/ref1_1.fq.gz"

workflow {
    FASTP( Channel.fromPath( params.input ) ).out.view()
}

process FASTP {

   input:
   path read

   output:
   stdout

   script:
   """
   fastp -A -i ${read} -o out.fq 2>&1
   """
}

