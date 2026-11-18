process RUN_KRAKEN2 {
    conda "${moduleDir}/environment.yml"
    cache 'lenient'
    label 'process_8'
    tag "${id}"
    //errorStrategy 'ignore'

	publishDir "results/BasicGenomics/${id}", \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

    input:
    tuple val(id), path(fwread), path(rvread)
    val(gtdbversion)
    path(gtdb_path)
    
    output:
    tuple val(id), path("${id}_kraken2_gtdbr${params.gtdbversion}.tab"), emit: krakenout
    
    script:
    """
	kraken2 \
		--db ${gtdb_path} \
		--threads ${task.cpus} \
		--output - \
		--memory-mapping \
		--report ${id}_kraken2_gtdbr${gtdbversion}.tab \
		--paired ${fwread} ${rvread} \
		2> /dev/null
    """
}