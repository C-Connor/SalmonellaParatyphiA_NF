process RUN_SHOVILL {
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
    
    output:
    path("${id}_shovill/"), emit: shovilldir
	tuple val(id), path("${id}_contigs.fa"), emit: contigs 
    
    script:
    """
	shovill \
		--cpus ${task.cpus} \
		--outdir ${id}_shovill \
		--R1 ${fwread} --R2 ${rvread} \
		--minlen 200 \
		&& \
		cp ${id}_shovill/contigs.fa ${id}_contigs.fa
    """
}