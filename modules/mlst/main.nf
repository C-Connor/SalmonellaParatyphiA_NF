process RUN_MLST {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/BasicGenomics/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_mlst.tsv"), emit: mlsttsv
	
	script:
	"""
	mlst --threads ${task.cpus} \
		--nopath \
		${contigs} \
		> ${id}_mlst.tsv
	"""
}

process RUN_MLST_LONG {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/LongRead/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_mlst.tsv"), emit: mlsttsv
	
	script:
	"""
	mlst --threads ${task.cpus} \
		--nopath \
		${contigs} \
		> ${id}_mlst.tsv
	"""
}