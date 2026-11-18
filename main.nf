include {BasicGenomics_shortRead;BasicGenomics_longAssembly} from './workflows/BasicGenomics'
include {BuildTree} from './workflows/BuildTree'
include {Pangenome} from './workflows/Pangenome'

sample_reads = Channel
		.fromPath(params.input_file)
		.splitCsv ( sep:'\t', header:['id','fastq1','fastq2'] )
		.map { row ->
			[row.id, 
			file(row.fastq1, checkIfExists: true), 
            file(row.fastq2, checkIfExists: true)
			]
		}

workflow {
    main:
    BasicGenomics_shortRead(
        sample_reads,
        params.bakta_db,
        file("${projectDir}/assets/Paratype/paratype.py",checkIfExists:true),
        file("${projectDir}/assets/Paratype/SParatyphiAKU12601.fasta",checkIfExists:true), // paratype ref genome
        file(params.snippy_ref,checkIfExists:true)
    )
    
    if (params.tree_samples_list) {
        log.info "making tree alignment"
        tree_samples_file = file(params.tree_samples_list,checkIfExists:true)
        tree_samples_id = Channel
                            .fromPath(tree_samples_file)
                            .splitCsv()
        //filter to just samples in tree list file
        tree_samples = BasicGenomics_shortRead.out.genomic_result_files
                            .join(tree_samples_id)

        //tree_samples_id.view() // [AUSMDU00122741]
        //BasicGenomics_shortRead.out.genomic_result_files.view()

        BuildTree(
            tree_samples,
            file(params.snippy_ref,checkIfExists:true),
            file(params.snippy_mask,checkIfExists:true)
        )
        
    } else {
        log.info "No list of IDs for tree alignment given."
    }

    if (params.long_read_assemblies) {
        long_assemblies = Channel
		    .fromPath(params.long_read_assemblies)
		    .splitCsv ( sep:'\t', header:['id','assembly'] )
		    .map { row ->
			    [row.id, 
			    file(row.assembly, checkIfExists: true)
			    ]
		    }
        
        BasicGenomics_longAssembly(
            long_assemblies,
            params.bakta_db
        )
    }

}