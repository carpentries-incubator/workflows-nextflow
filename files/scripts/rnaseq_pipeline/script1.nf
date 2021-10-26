/*
 * pipeline input parameters
 */

params.reads = "$projectDir/data/yeast/reads/*_{1,2}.fq.gz"
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.multiqc = "$projectDir/multiqc"

println "reads: $params.reads"
