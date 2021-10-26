csv_ch = Channel.fromPath( 'data/yeast/samples.csv' )
        .splitCsv( header:true )

csv_ch.view( { it.sample_id } )
