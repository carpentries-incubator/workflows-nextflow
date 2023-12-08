ids = ['ERR908507', 'ERR908506', 'ERR908505']

queue_ch = Channel.fromList(ids)
value_ch = Channel.value(ids)
queue_ch.view()
value_ch.view()
