process RUN_ABRITAMR {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/BasicGenomics/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_abritamr"), emit: abritamrdir
	
	script:
	"""
	abritamr run \
		-c ${contigs} \
		-px ${id}_abritamr \
		-j ${task.cpus} \
		--species Salmonella
	"""
}

process RUN_ABRITAMR_LONG {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/LongRead/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_abritamr"), emit: abritamrdir
	
	script:
	"""
	abritamr run \
		-c ${contigs} \
		-px ${id}_abritamr \
		-j ${task.cpus} \
		--species Salmonella
	"""
}