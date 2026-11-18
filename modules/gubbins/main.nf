//get an error "cannot cache function"
//can fix with `export NXF_APPTAINER_HOME_MOUNT=true`
process RUN_GUBBINS {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8' //or more?
	//tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/tree/Gubbins/", \
		mode: params.publish_dir_mode, \
		saveAs: { filename -> "${filename}" }

	input:
	path(aln)
	
	output:
	path("CoreFull_Gubbins*"), emit: gubbinsdir
	path("CoreFull_Gubbins.recombination_predictions.gff"), emit: gubbinsgff

	script:
	//export home to fix "cannot cache function"
	"""
	export HOME=\$(pwd)

	run_gubbins.py \
		--prefix CoreFull_Gubbins \
    	--threads ${task.cpus} \
    	--tree-builder raxml \
    	--converge-method recombination \
    	--model GTRGAMMA \
    	${aln}
	"""
}

process MASK_GUBBINS_RECOM {
	conda "${moduleDir}/environment.yml"
	cache 'lenient'
	label 'process_8' //or more?
	//tag "${id}"
	//errorStrategy 'ignore'

	publishDir "results/tree/Gubbins/", \
		mode: params.publish_dir_mode, \
		saveAs: { filename -> "${filename}" }

	input:
	path(aln)
	path(gubbinsgff)
	
	output:
	path("GubbinsRecombMasked.aln"), emit: gubbinsmask

	script:
	//export home to fix "cannot cache function"
	"""
	mask_gubbins_aln.py \
		--aln ${aln} \
		--gff ${gubbinsgff} \
		--out GubbinsRecombMasked.aln
	"""
}