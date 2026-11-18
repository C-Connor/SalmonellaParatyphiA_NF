process RUN_PARATYPE {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	errorStrategy 'retry'
	maxRetries 3
	memory { 4.GB * task.attempt }

	publishDir "results/BasicGenomics/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

	input:
    path(paratypescript)
    path(paratyperef)
	tuple val(id), path(fwread), path(rvread)
	
	output:
	tuple val(id), path("${id}_paratype.txt"), emit: paratypetxt
	
	script:
	"""
    python ${paratypescript} \
        --id ${id} \
        --mode fastq \
        --fastq ${fwread} ${rvread} \
        --threads ${task.cpus} \
        --output ${id}_paratype.txt \
		--ref ${paratyperef}
	"""
}