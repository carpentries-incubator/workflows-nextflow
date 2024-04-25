ch = Channel.fromPath( 'data/yeast/reads/*.fq.gz' )
    .map ( { file -> [ file, file.getName() ] } )
    .view( { file, name -> "file's name: $name" } )
