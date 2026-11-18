include {RUN_FASTP} from '../modules/fastp'
include {RUN_KRAKEN2} from '../modules/kraken2'
include {RUN_SHOVILL} from '../modules/shovill'
include {RUN_QUAST} from '../modules/quast'
include {RUN_BAKTA;RUN_BAKTA_LONG} from '../modules/bakta'
include {RUN_ABRITAMR;RUN_ABRITAMR_LONG} from '../modules/abritamr'
include {RUN_MLST;RUN_MLST_LONG} from '../modules/mlst'
include {RUN_MOBRECON;RUN_MOBRECON_LONG} from '../modules/mobsuite'
include {RUN_SISTR;RUN_SISTR_LONG} from '../modules/sistr'
include {RUN_SNIPPY;RUN_SAMTOOLSDEPTH} from '../modules/snippy'
include {RUN_PARATYPE} from '../modules/paratype'


workflow BasicGenomics_shortRead {
    take:
    sample_reads
    bakta_db
    paratype_script
    paratype_refgenome
    snippy_ref

    main:
    RUN_FASTP(
        sample_reads
    )
    RUN_KRAKEN2(
        RUN_FASTP.out.fastpreads,
        params.gtdbversion,
        params.gtdb_path
    )
    RUN_SHOVILL(
        RUN_FASTP.out.fastpreads
    )
    RUN_QUAST(
        RUN_SHOVILL.out.contigs
    )
    RUN_BAKTA(
        RUN_SHOVILL.out.contigs,
        bakta_db
    )
    RUN_ABRITAMR(
        RUN_SHOVILL.out.contigs
    )
    RUN_MLST(
        RUN_SHOVILL.out.contigs
    )
    RUN_MOBRECON(
        RUN_SHOVILL.out.contigs
    )
    RUN_SISTR(
        RUN_SHOVILL.out.contigs
    )
    // paratype puts the reference index files in the assets/Paratype dir for some reason
    RUN_PARATYPE(
        paratype_script,//"${projectDir}/assets/Paratype/paratype.py", 
        paratype_refgenome, //"${projectDir}/assets/Paratype/SParatyphiAKU12601.fasta",
        RUN_FASTP.out.fastpreads
    )
    RUN_SNIPPY(
        sample_reads,
        snippy_ref
    ) 
    RUN_SAMTOOLSDEPTH(
        RUN_SNIPPY.out.snippyrun
    )

    // emit all genomic results together at once?
    // qc check each sample as it's completed? - need to generate csv ouput of combined data
    // could write individual lines out to results/${id} and merge all lines at end?
    // collectFile -> tab file with sample ID and filepaths of results
    // python QC script to collect and process results
    genomic_result_files = RUN_FASTP.out.fastpreads
                                .join(RUN_KRAKEN2.out.krakenout)
                                .join(RUN_SHOVILL.out.contigs)
                                .join(RUN_QUAST.out.quastdir)
                                .join(RUN_BAKTA.out.baktagff3)
                                .join(RUN_ABRITAMR.out.abritamrdir)
                                .join(RUN_MLST.out.mlsttsv)
                                .join(RUN_MOBRECON.out.mobrecondir)
                                .join(RUN_SISTR.out.sistrcsv)
                                .join(RUN_PARATYPE.out.paratypetxt)
                                .join(RUN_SNIPPY.out.snippyrun)
                                .join(RUN_SAMTOOLSDEPTH.out)

    emit:
    //emit QC passed samples
    genomic_result_files
}

// repeat basic genomics for hybrid assembled samples
workflow BasicGenomics_longAssembly {
    take:
    long_assemblies
    bakta_db

    main:
    RUN_ABRITAMR_LONG(
        long_assemblies
    )

    RUN_MLST_LONG(
        long_assemblies
    )

    RUN_MOBRECON_LONG(
        long_assemblies
    )

    RUN_SISTR_LONG(
        long_assemblies
    )

    RUN_BAKTA_LONG(
        long_assemblies,
        bakta_db
    )

    genomic_result_files = RUN_ABRITAMR_LONG.out.abritamrdir
                                .join(RUN_MLST_LONG.out.mlsttsv)
                                .join(RUN_MOBRECON_LONG.out.mobrecondir)
                                .join(RUN_SISTR_LONG.out.sistrcsv)
                                .join(RUN_BAKTA_LONG.out.baktagff3)
    
    emit:
    genomic_result_files
}