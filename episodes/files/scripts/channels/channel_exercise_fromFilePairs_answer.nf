pairs_ch = Channel.fromFilePairs('data/yeast/reads/temp33*_{1,2}.fq.gz')
pairs_ch.view()
