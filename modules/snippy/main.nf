process RUN_SNIPPY {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/tree/snippy/${id}", \
		mode: params.publish_dir_mode, \
		saveAs: { filename -> "${filename}" }

	input:
	tuple val(id), path(fwread), path(rvread)
    path(refgenome)
	
	output:
	tuple val(id), path("${id}_snippyOut"), emit: snippyrun
	
	script:
	"""
    snippy \
        --cpus ${task.cpus} \
        --outdir ${id}_snippyOut \
        --ref ${refgenome} \
        --R1 ${fwread} \
        --R2 ${rvread} \
        --mincov 10 \
        --minfrac 0.9
	"""
}

process RUN_SNIPPYCORE {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8'
	//tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/tree/snippycore", \
		mode: params.publish_dir_mode, \
		saveAs: { filename -> "${filename}" }

	input:
	path(snippy)
	path(phagemask)
    path(refgenome)
	
	output:
	path("core.full.aln"), emit: snippycore
	path("core.*") //for publishing
	
	script:
	"""
    snippy-core \
		--mask ${phagemask} \
		--mask-char n \
		--ref ${refgenome} \
		${snippy.join(' ')}
	"""
}

process RUN_SAMTOOLSDEPTH {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_1'
	tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/read_coverage/", \
		mode: params.publish_dir_mode, \
		saveAs: { filename -> "${filename}" }

	input:
	tuple val(id), path(snippy)
	
	output:
	tuple val(id), path("${id}_depth.tab")
	
	script:
	"""
    samtools depth \
		-a \
		-H \
		-o ${id}_depth.tab \
		${snippy}/snps.bam
	"""
}