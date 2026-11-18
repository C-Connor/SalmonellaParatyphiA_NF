process RUN_QUAST {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	//errorStrategy 'ignore'
	
	publishDir "results/BasicGenomics/${id}", \
		mode: params.publish_dir_mode
		//saveAs: {filename -> "${id}_${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_quast/"), emit: quastdir
	
	script:
	"""
	quast \
		-o ${id}_quast/ \
		--fast \
		-t ${task.cpus} \
		${contigs}
	"""
}