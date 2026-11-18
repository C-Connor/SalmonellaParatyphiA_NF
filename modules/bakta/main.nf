// have to run Bakta in base conda env
// running through wave container causes RAM to steadily climb until pipeline is killed
// unclear why
// have tried setting memory and cpu constraints - no difference
process RUN_BAKTA {
	//conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	//cpus 4
	tag "${id}"
	//errorStrategy 'ignore'
	//memory 8.GB

	publishDir "results/BasicGenomics/${id}/bakta", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${filename}" }

	input:
	tuple val(id), path(contigs)
    path(bakta_db)
	
	output:
	tuple val(id), path("${id}_bakta.gff3"), emit: baktagff3
	path("${id}_bakta.*")
	
	script:
	"""
    bakta \
        ${contigs} \
        --db ${bakta_db} \
        --prefix ${id}_bakta \
        --threads ${task.cpus}
	"""
}

process RUN_BAKTA_LONG {
	//conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	//cpus 4
	tag "${id}"
	//errorStrategy 'ignore'
	//memory 8.GB

	publishDir "results/LongRead/${id}/bakta", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${filename}" }

	input:
	tuple val(id), path(contigs)
    path(bakta_db)
	
	output:
	tuple val(id), path("${id}_bakta.gff3"), emit: baktagff3
	path("${id}_bakta.*")
	
	script:
	"""
    bakta \
        ${contigs} \
        --db ${bakta_db} \
        --prefix ${id}_bakta \
        --threads ${task.cpus}
	"""
}