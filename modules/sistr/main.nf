process RUN_SISTR {
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
	tuple val(id), path("${id}_sistr_output.csv"), emit: sistrcsv
    path("${id}_sistr_fullOutput")
	
	script:
	"""
    sistr \
        --output-format csv \
        --output-prediction ${id}_sistr_output.csv \
        --cgmlst-profiles ${id}_sistr_cgMLSTprofiles.csv \
        --novel-alleles ${id}_sistr_novelAlleles.fasta \
        --alleles-output ${id}_sistr_outputAlleles.json \
        --run-mash \
        --qc \
        --threads ${task.cpus} \
		${contigs}
    mkdir ${id}_sistr_fullOutput
    cp ${id}_sistr_cgMLSTprofiles.csv ${id}_sistr_fullOutput/
    cp ${id}_sistr_novelAlleles.fasta ${id}_sistr_fullOutput/
    cp ${id}_sistr_outputAlleles.json ${id}_sistr_fullOutput/
	"""
}

process RUN_SISTR_LONG {
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
	tuple val(id), path("${id}_sistr_output.csv"), emit: sistrcsv
    path("${id}_sistr_fullOutput")
	
	script:
	"""
    sistr \
        --output-format csv \
        --output-prediction ${id}_sistr_output.csv \
        --cgmlst-profiles ${id}_sistr_cgMLSTprofiles.csv \
        --novel-alleles ${id}_sistr_novelAlleles.fasta \
        --alleles-output ${id}_sistr_outputAlleles.json \
        --run-mash \
        --qc \
        --threads ${task.cpus} \
		${contigs}
    mkdir ${id}_sistr_fullOutput
    cp ${id}_sistr_cgMLSTprofiles.csv ${id}_sistr_fullOutput/
    cp ${id}_sistr_novelAlleles.fasta ${id}_sistr_fullOutput/
    cp ${id}_sistr_outputAlleles.json ${id}_sistr_fullOutput/
	"""
}