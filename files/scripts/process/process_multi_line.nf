//process_multi_line.nf
nextflow.enable.dsl=2

process NUMSEQ_CHR {
    script:
    """
    zgrep  '^>' ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz > ids.txt
    zgrep -c '>YA' ids.txt
    """
}

workflow {
  NUMSEQ_CHR()
}

