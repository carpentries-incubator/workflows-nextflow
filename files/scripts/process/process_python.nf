nextflow.enable.dsl = 2

process PYSCRIPT {

    script:
    """
    #!/usr/bin/env python
    import gzip

    reads = 0
    bases = 0

    with gzip.open('${projectDir}/data/yeast/reads/ref1_1.fq.gz', 'rb') as read:
      for id in read:
          seq = next(read)
          reads += 1
          bases += len(seq.strip())
          next(read)
          next(read)

    print("reads", reads)
    print("bases", bases)
    """
}

workflow {
    PYSCRIPT()
}
