process RUN_MOBRECON {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_1'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/BasicGenomics/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_mobrecon"), emit: mobrecondir
	
	script:
	"""
	mv ${contigs} ${id}_${contigs}
    mob_recon \
        --infile ${id}_${contigs} \
        --outdir ${id}_mobrecon
	"""
}

process RUN_MOBRECON_LONG {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_1'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/LongRead/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

	input:
	tuple val(id), path(contigs)
	
	output:
	tuple val(id), path("${id}_mobrecon"), emit: mobrecondir
	
	script:
	"""
	mv ${contigs} ${id}_${contigs}
    mob_recon \
        --infile ${id}_${contigs} \
        --outdir ${id}_mobrecon
	"""
}