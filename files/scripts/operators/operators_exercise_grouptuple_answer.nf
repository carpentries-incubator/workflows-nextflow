ch = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )
    .map { file -> [ file.getName().split('_')[0], file ] }
    .groupTuple()
    .view()
