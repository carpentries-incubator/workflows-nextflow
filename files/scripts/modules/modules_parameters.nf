nextflow.enable.dsl = 2

params.message = 'parameter from workflow script'

include { sayMessage } from './modules/module.nf'

workflow {
    sayMessage()
}
