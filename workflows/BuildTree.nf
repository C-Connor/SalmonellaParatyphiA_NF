include {RUN_GUBBINS;MASK_GUBBINS_RECOM} from '../modules/gubbins'
include {RUN_SNIPPY;RUN_SNIPPYCORE;RUN_SAMTOOLSDEPTH} from '../modules/snippy'

workflow BuildTree {
    take:
    filtered_genomics //id, fastpfw, fastprv, kraken2, contigs, quastDir, abritamrDir, mlst, mobrecon, sistr, paratype, consensus
    snippy_ref
    phage_mask
    
    main:
    snippy_runs = filtered_genomics
                    .map{it[12]} //just take snippy 
                    .collect()

    RUN_SNIPPYCORE(
        snippy_runs,
        phage_mask,
        snippy_ref
    )
    RUN_GUBBINS(
        RUN_SNIPPYCORE.out.snippycore
    )
    RUN_SAMTOOLSDEPTH(
        RUN_SNIPPY.out.snippyrun
    )
    
    MASK_GUBBINS_RECOM(
        RUN_SNIPPYCORE.out.snippycore,
        RUN_GUBBINS.out.gubbinsgff
    )

    //emit:
    //gubbinsoutput?
}