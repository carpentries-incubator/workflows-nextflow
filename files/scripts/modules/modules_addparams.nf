nextflow.enable.dsl = 2

params.message = 'parameter from workflow script'

include { sayMessage } from './modules/module.nf' addParams( message: 'using addParams' )

workflow {
    sayMessage()
}
