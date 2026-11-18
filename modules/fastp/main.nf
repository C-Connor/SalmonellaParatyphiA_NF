process RUN_FASTP {
    conda "${moduleDir}/environment.yml"
    cache 'lenient'
    label 'process_8'
    tag "${id}"
    errorStrategy 'ignore'

    publishDir "results/BasicGenomics/${id}",\
		pattern: '*.fastq.gz', \
		mode: params.publish_dir_mode_fastpreads
		//saveAs: { filename -> "${id}_${filename}" }
	
	publishDir "results/BasicGenomics/${id}", \
		pattern: '*fastp.*', \
		mode: params.publish_dir_mode
		//saveAs: { filename -> "${id}_${filename}" }

    input:
    tuple val(id), path(fwread), path(rvread)
    
    output:
    tuple val(id), \
		path("${id}_R1.trimmed.fastq.gz"), \
		path("${id}_R2.trimmed.fastq.gz"), \
		emit: fastpreads
	tuple val(id), path("${id}_fastp.json"), emit: fastpjson
	tuple val(id), path("${id}_fastp.html"), emit: fastphtml
    
    script:
    """
    fastp \
		-w 3 \
		-r \
		-W 4 \
		-M 20 \
		-l 36 \
		-j ${id}_fastp.json \
		-h ${id}_fastp.html \
		-R '${id}' \
		-i ${fwread} \
		-I ${rvread} \
		-o ${id}_R1.trimmed.fastq.gz \
		-O ${id}_R2.trimmed.fastq.gz
    """
}