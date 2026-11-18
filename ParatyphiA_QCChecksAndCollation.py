import marimo

__generated_with = "0.15.2"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    import pandas as pd
    import numpy as np
    import os, sys
    from matplotlib import pyplot as plt
    import seaborn as sns
    import glob
    from pathlib import Path
    import json
    return Path, glob, json, mo, pd, sns


@app.cell
def _(sns):
    sns.set_theme(style='white')
    return


@app.cell
def _():
    cm = 1/2.54
    SenteriaPTyhiAGenomeSize = 4581797 #https://www.ncbi.nlm.nih.gov/nuccore/FM200053.1/
    return (SenteriaPTyhiAGenomeSize,)


@app.cell
def _(pd):
    MDUmetadata = pd.read_excel("../AUSMDUids_Paratyphi_ONT_repeat_travel_14112024.xlsx").add_prefix("MDULIMS_")

    MDUmetadata.rename(columns={'MDULIMS_AUSMDU':"SampleID"},inplace=True)
    return (MDUmetadata,)


@app.cell
def _(MDUmetadata):
    MDUmetadata
    return


@app.cell
def _(MDUmetadata):
    #Fix Missing and missing in MDULIMS_Year column
    MDUmetadata['MDULIMS_Year'] = MDUmetadata['MDULIMS_Year'].replace('missing','Missing')
    return


@app.cell
def _():
    #Collect genomics
    return


@app.cell
def _(MDUID, Path, json, pd):
    def CollectFastp(filepath):
        id = Path(filepath).name
        with open(f"{filepath}/{id}_fastp.json","r") as filehandle:
            indfastpjson = json.load(filehandle)
        indfastp = pd.DataFrame.from_dict({
            "BeforeFiltering":indfastpjson['summary']['before_filtering'],
            "AfterFiltering":indfastpjson['summary']['after_filtering']
        },orient='index')
        indfastpresults = pd.DataFrame(data={
            "SampleID":id,
            "PercentReadsFilterPass":round(indfastp.loc['AfterFiltering','total_reads'] / indfastp.loc['BeforeFiltering','total_reads'] * 100,4),
            "TotalFilteredReads":indfastp.loc['AfterFiltering','total_reads'],
            "TotalFilteredBases":indfastp.loc['AfterFiltering','total_bases'],
            "Q30Rate":indfastp.loc['AfterFiltering','q30_rate'],
            "GCContent":indfastp.loc['AfterFiltering','gc_content']
        },index=[1])
        return indfastpresults

    def CollectKraken(filepath,krakendbversion):
        id = Path(filepath).name
        indkraken = pd.read_table(f"{filepath}/{id}_kraken2_{krakendbversion}.tab",header=None)
        indkraken = indkraken.loc[indkraken[3]=="S",[5,0]].sort_values(by=0,ascending=False).head(3).copy()
        indkraken[5] = indkraken[5].str.strip().str.replace("s__","")
        indkraken = indkraken.reset_index(drop=True).rename(
            columns={0:"Percentage",5:"Species"},
            index={0:"FirstMatch",1:"SecondMatch",2:"ThirdMatch"}).unstack().to_frame().sort_index(level=1).T
        indkraken.columns = indkraken.columns.map('_'.join)
        indkraken['SampleID'] = id
        #return pd.concat([AllKraken,indkraken],ignore_index=True)
        return indkraken

    def CollectAssemblyStats(filepath):
        id = Path(filepath).name
        try:
            indassemblystats = pd.read_table(f"{filepath}/{id}_quast/transposed_report.tsv")
            indassemblystats = indassemblystats.loc[:,[#"Assembly",
                "# contigs","Total length","N50","auN","# N's per 100 kbp"]]
            indassemblystats = indassemblystats.rename(columns={"# contigs":"NumContigs","Total length":"TotalLength","# N's per 100 kbp":"NumNsPer100kbp"})
            indassemblystats['SampleID'] = id
            return indassemblystats
        except FileNotFoundError:
            indassemblystats = pd.DataFrame({
                'SampleID':[id],
                'NumContigs':[0],
                'TotalLength':[0],
                'N50':[0],
                'auN':[0],
                'NumNsPer100kbp':[0]
            })
            return indassemblystats

    def CollectMLST(filepath):
        id = Path(filepath).name
        try:
            indmlst = pd.read_table(f"{filepath}/{id}_mlst.tsv",header=None)
            indmlst["MLST_alleles"] = indmlst.loc[:,3:].agg(";".join,axis=1)
            indmlst['SampleID'] = id
            indmlst = indmlst.rename(columns={0:"MLSTContigs",1:"MLSTScheme",2:"ST"})
            return indmlst.loc[:,['SampleID','MLSTScheme','ST','MLST_alleles']]
        except FileNotFoundError:
            indmlst = pd.DataFrame({
                'SampleID':[id],
                'MLSTScheme':[""],
                'ST':[""],
                'MLST_alleles':[""]
            })
            return indmlst

    def CollectAMR(filepath):
        id = Path(filepath).name
        try:
            indamr = pd.read_table(f"{filepath}/{id}_abritamr/abritamr.txt") #file not always present?
            indamr['SampleID'] = id
            return indamr
        except FileNotFoundError:
            indamr = pd.DataFrame({
                'Isolate':[f"{MDUID}_abritamr"],
                'SampleID':[id]
            })
            return indamr

    def CollectMOBRecon(filepath):
        #Just collect num plasmid contigg, num of cluster_ids
        id = Path(filepath).name
        try:
            mobrecon = pd.read_table(f"{filepath}/{id}_mobrecon/contig_report.txt")
        except FileNotFoundError:
            mobresult = pd.DataFrame({
                'SampleID':[id],
                'NumPlasmidContigs':[],
                'UniqueMOBClusterIDs':[]
            })
        numplasmidcontigs = mobrecon.loc[mobrecon['molecule_type'].str.contains('plasmid')].shape[0]
        clusterids = ';'.join(set(mobrecon.loc[mobrecon['molecule_type'].str.contains('plasmid'),'primary_cluster_id']))
        mobresult = pd.DataFrame({
            'SampleID':[id],
            'NumPlasmidContigs':[numplasmidcontigs],
            'UniqueMOBClusterIDs':[clusterids]
        })
        return mobresult

    def CollectParatype(filepath):
        id = Path(filepath).name
        try:
            paratype = pd.read_table(f"{filepath}/{id}_paratype.txt")
        except FileNotFoundError:
            paratype = pd.DataFrame({
                'SampleID':[id],
                'Primary_clade':[],
                'Secondary_clade':[],
                'Subclade':[],
                'Genotype':[],
                'Support':[],
                'Mutations':[]
            })
        paratype = paratype.add_prefix("ParaType_")
        paratype.rename(columns={'ParaType_Strain':'SampleID'},inplace=True)
        return paratype

    def CollectSistr(filepath):
        id = Path(filepath).name
        try:
            sistr = pd.read_csv(f"{filepath}/{id}_sistr_output.csv")
        except FileNotFoundError:
            sistr = pd.DataFrame({'SampleID':[id]})
        sistr = sistr.add_prefix("Sistr_")
        sistr['SampleID'] = id
        return sistr

    def CollectContigReport(filepath):
        id = Path(filepath).name
        try:
            contigreport = pd.read_table(f"{filepath}/{id}_mobrecon/contig_report.txt")
        except FileNotFoundError:
            contigreport = pd.DataFrame({'SampleID':[id]})
        contigreport['SampleID'] = id
        return contigreport

    def CollectVirulence(filepath):
        id = Path(filepath).name
        try:
            virulencereport = pd.read_table(f"{filepath}/{id}_abritamr/summary_virulence.txt")
        except FileNotFoundError:
            virulencereport = pd.DataFrame({"SampleID":[id]})
        virulencereport['SampleID'] = id
        return virulencereport
    return (
        CollectAMR,
        CollectAssemblyStats,
        CollectContigReport,
        CollectFastp,
        CollectKraken,
        CollectMLST,
        CollectMOBRecon,
        CollectParatype,
        CollectSistr,
        CollectVirulence,
    )


@app.cell
def _(
    CollectAMR,
    CollectAssemblyStats,
    CollectContigReport,
    CollectFastp,
    CollectKraken,
    CollectMLST,
    CollectMOBRecon,
    CollectParatype,
    CollectSistr,
    CollectVirulence,
    Path,
    glob,
    mo,
    pd,
):
    resultsdir = '../Nextflow/results/BasicGenomics'
    SampleIDlist = [Path(x).name for x in glob.glob(f"{resultsdir}/AUSMDU*")]

    FastpList = []
    Kraken2List = []
    AssemblyStatsList = []
    MLSTList = []
    AMRList = []
    MOBList = []
    ParatypeList = []
    SistrList = []
    mobRecon_contigreportlist = []
    VirulenceList = []

    krakendbversion = "gtdbr214"

    for fp in mo.status.progress_bar(glob.glob(f"{resultsdir}/AUSMDU*")):
        #FastP
        FastpList.append(CollectFastp(fp))
        #Kraken2
        Kraken2List.append(CollectKraken(fp,krakendbversion))
        #Quast
        AssemblyStatsList.append(CollectAssemblyStats(fp))
        #MLST
        MLSTList.append(CollectMLST(fp))
        #Abritamr
        AMRList.append(CollectAMR(fp))
        #MOBSuite
        MOBList.append(CollectMOBRecon(fp))
        #Paratype
        ParatypeList.append(CollectParatype(fp))
        #Sistr
        SistrList.append(CollectSistr(fp))
        #MOB Recon contig report
        mobRecon_contigreportlist.append(CollectContigReport(fp))
        #Virulence
        VirulenceList.append(CollectVirulence(fp))


    CombinedFastP = pd.concat([x for x in FastpList],ignore_index=True)
    CombinedKraken2 = pd.concat([x for x in Kraken2List],ignore_index=True)
    CombinedAssembyStats = pd.concat([x for x in AssemblyStatsList],ignore_index=True)
    CombinedMLST = pd.concat([x for x in MLSTList],ignore_index=True)
    CombinedAMR = pd.concat([x for x in AMRList],ignore_index=True)
    CombinedMOB = pd.concat([x for x in MOBList],ignore_index=True)
    CombinedParatype = pd.concat([x for x in ParatypeList],ignore_index=True)
    CombinedSistr = pd.concat([x for x in SistrList],ignore_index=True)
    CombinedContigReport = pd.concat([x for x in mobRecon_contigreportlist],ignore_index=True)
    CombinedVirulence = pd.concat([x for x in VirulenceList],ignore_index=True)
    return (
        CombinedAMR,
        CombinedAssembyStats,
        CombinedContigReport,
        CombinedFastP,
        CombinedKraken2,
        CombinedMLST,
        CombinedMOB,
        CombinedParatype,
        CombinedSistr,
        CombinedVirulence,
        SampleIDlist,
    )


@app.cell
def _(CombinedContigReport, SampleIDlist, pd):
    plasmid_rep = CombinedContigReport.loc[
        CombinedContigReport['rep_type_accession(s)'] != '-'
    ].groupby('SampleID').agg(
        RepGeneContigIDs = ('contig_id',lambda x: ','.join(x.str[0:11])),
        RepGeneAccessions = ('rep_type_accession(s)', lambda x: ','.join(x))
    )

    plasmid_relaxase = CombinedContigReport.loc[
        CombinedContigReport['relaxase_type_accession(s)'] != '-'
    ].groupby('SampleID').agg(
        RelaxaseGeneContigIDs = ('contig_id',lambda x: ','.join(x.str[0:11])),
        RelaxaseGeneAccessions = ('relaxase_type_accession(s)', lambda x: ','.join(x))
    )

    plasmid_mpf = CombinedContigReport.loc[
        CombinedContigReport['mpf_type_accession(s)'] != '-'
    ].groupby('SampleID').agg(
        MpfGeneContigIDs = ('contig_id',lambda x: ','.join(x.str[0:11])),
        MpfGeneAccessions = ('mpf_type_accession(s)', lambda x: ','.join(x))
    )

    plasmid_orit = CombinedContigReport.loc[
        CombinedContigReport['orit_accession(s)'] != '-'
    ].groupby('SampleID').agg(
        OriTGeneContigIDs = ('contig_id',lambda x: ','.join(x.str[0:11])),
        OriTGeneAccessions = ('orit_accession(s)', lambda x: ','.join(x))
    )

    #CombinedPlasmidGenes = pd.concat([plasmid_rep,plasmid_relaxase,plasmid_mpf,plasmid_orit]).reset_index()
    CombinedPlasmidGenes = pd.DataFrame({'SampleID':SampleIDlist})
    for z in [plasmid_rep,plasmid_relaxase,plasmid_mpf,plasmid_orit]:
        CombinedPlasmidGenes = CombinedPlasmidGenes.merge(z,on='SampleID',how='outer')
    return (CombinedPlasmidGenes,)


@app.cell
def _(
    CombinedAMR,
    CombinedAssembyStats,
    CombinedFastP,
    CombinedKraken2,
    CombinedMLST,
    CombinedMOB,
    CombinedParatype,
    CombinedPlasmidGenes,
    CombinedSistr,
    CombinedVirulence,
    SampleIDlist,
    SenteriaPTyhiAGenomeSize,
    pd,
):
    GenomicData = pd.DataFrame({'SampleID':SampleIDlist})
    for x in [CombinedFastP, CombinedKraken2, CombinedAssembyStats, CombinedMLST, CombinedAMR, CombinedMOB, CombinedParatype, CombinedSistr, CombinedPlasmidGenes, CombinedVirulence]:
        GenomicData = GenomicData.merge(x, on="SampleID",how='outer') 

    GenomicData['SeqCoverage'] = GenomicData['TotalFilteredBases'] / SenteriaPTyhiAGenomeSize
    return (GenomicData,)


@app.cell
def _(GenomicData):
    GenomicData
    return


@app.cell
def _(GenomicData, SenteriaPTyhiAGenomeSize):
    #QC Checks
    QCFail = {}

    #Lots of reads trimmed / removed
    QCFail['High read loss'] = [x for x in GenomicData.loc[GenomicData['PercentReadsFilterPass'] < 70,'SampleID'].values]
    #Low coverage
    QCFail['Low sequencing coverage'] = [x for x in GenomicData.loc[GenomicData['SeqCoverage'] < 30,'SampleID'].values]
    #Not Salmonella enterica
    QCFail['Not S. enterica'] = [x for x in GenomicData.loc[GenomicData['Species_FirstMatch'] != 'Salmonella enterica','SampleID'].values]
    #Possible contamination
    QCFail['High second Kraken hit'] = [x for x in GenomicData.loc[GenomicData['Percentage_SecondMatch'] > 15,'SampleID'].values]
    #High number of contigs
    QCFail['High number of contigs'] = [x for x in GenomicData.loc[GenomicData['NumContigs'] > 500,'SampleID'].values]
    #Not paratyphi A
    QCFail['Not Paratyphi A'] = [x for x in GenomicData.loc[GenomicData['Sistr_serovar'] != 'Paratyphi A','SampleID'].values]
    #Genome size deviates
    QCFail['Genome size deviation over 20percent'] = [x for x in GenomicData.loc[~GenomicData['TotalLength'].between(SenteriaPTyhiAGenomeSize*0.8,SenteriaPTyhiAGenomeSize*1.2),'SampleID'].values]

    return (QCFail,)


@app.cell
def _(GenomicData):
    GenomicData
    return


@app.cell
def _(QCFail):
    QCFail
    return


@app.cell
def _(GenomicData, QCFail):
    QCFailInv = {}
    for Reason, ids in QCFail.items():
        for id in ids:
            if id in QCFailInv.keys():
                QCFailInv[id].append(Reason)
            else:
                QCFailInv[id] = [Reason]

    ExclusionMappingDict = {k:"; ".join(v) for k,v in QCFailInv.items()}

    GenomicData['InclusionExclusion'] = 'included'
    GenomicData.loc[GenomicData['SampleID'].isin(QCFailInv.keys()),'InclusionExclusion'] = 'EXCLUDED'
    GenomicData['ExclusionReason'] = GenomicData['SampleID'].map(ExclusionMappingDict)
    return


@app.cell
def _(GenomicData):
    GenomicData
    return


@app.cell
def _(GenomicData):
    #Write out sample that have passed QC
    GenomicData.loc[GenomicData['InclusionExclusion']=='included','SampleID'].to_csv('../QCPassedAUSMDUIDs.list',header=False,index=False)
    return


@app.cell
def _(GenomicData):
    #Write out full short read genomic data
    GenomicData.to_csv('../ShortReadCollatedNextflowResults.csv',header=True,index=False)
    return


@app.cell
def _():
    #Collect Long Read
    return


@app.cell
def _(
    CollectAMR,
    CollectMLST,
    CollectMOBRecon,
    CollectSistr,
    Path,
    glob,
    mo,
    pd,
):
    Longresultsdir = '../Nextflow/results/LongRead'
    LongSampleIDlist = [Path(x).name for x in glob.glob(f"{Longresultsdir}/DMG*")]

    #FastpList = []
    #Kraken2List = []
    #AssemblyStatsList = []
    LongMLSTList = []
    LongAMRList = []
    LongMOBList = []
    #ParatypeList = []
    LongSistrList = []

    #krakendbversion = "gtdbr214"

    for lfp in mo.status.progress_bar(glob.glob(f"{Longresultsdir}/DMG*")):
        #FastP
        #FastpList.append(CollectFastp(fp))
        #Kraken2
        #Kraken2List.append(CollectKraken(fp,krakendbversion))
        #Quast
        #AssemblyStatsList.append(CollectAssemblyStats(fp))
        #MLST
        LongMLSTList.append(CollectMLST(lfp))
        #Abritamr
        LongAMRList.append(CollectAMR(lfp))
        #MOBSuite
        LongMOBList.append(CollectMOBRecon(lfp))
        #Paratype
        #ParatypeList.append(CollectParatype(fp))
        #Sistr
        LongSistrList.append(CollectSistr(lfp))



    #CombinedFastP = pd.concat([x for x in FastpList],ignore_index=True)
    #CombinedKraken2 = pd.concat([x for x in Kraken2List],ignore_index=True)
    #CombinedAssembyStats = pd.concat([x for x in AssemblyStatsList],ignore_index=True)
    LongCombinedMLST = pd.concat([x for x in LongMLSTList],ignore_index=True)
    LongCombinedAMR = pd.concat([x for x in LongAMRList],ignore_index=True)
    LongCombinedMOB = pd.concat([x for x in LongMOBList],ignore_index=True)
    #CombinedParatype = pd.concat([x for x in ParatypeList],ignore_index=True)
    LongCombinedSistr = pd.concat([x for x in LongSistrList],ignore_index=True)
    return (
        LongCombinedAMR,
        LongCombinedMLST,
        LongCombinedMOB,
        LongCombinedSistr,
        LongSampleIDlist,
    )


@app.cell
def _(LongCombinedMLST):
    LongCombinedMLST
    return


@app.cell
def _(
    LongCombinedAMR,
    LongCombinedMLST,
    LongCombinedMOB,
    LongCombinedSistr,
    LongSampleIDlist,
    pd,
):
    LongGenomicData = pd.DataFrame({'SampleID':LongSampleIDlist})
    for y in [LongCombinedMLST, LongCombinedAMR, LongCombinedMOB, LongCombinedSistr]:
        LongGenomicData = LongGenomicData.merge(y, on="SampleID",how='outer') 
    return (LongGenomicData,)


@app.cell
def _(LongGenomicData):
    LongGenomicData
    return


@app.cell
def _(LongGenomicData):
    #Write out full short read genomic data
    LongGenomicData.to_csv('../LongAssemblyCollatedNextflowResults.csv',header=True,index=False)
    return


@app.cell
def _():
    #Process Sequencing depth 
    return


@app.cell
def _(glob, mo, pd):
    dflist = []
    rawlist = []
    for f in mo.status.progress_bar(glob.glob("../Nextflow/results/read_coverage/*_depth.tab")):
        depth = pd.read_table(f)
        #depth.drop(columns="#CHROM",inplace=True)
        sampleid = depth.columns[2].strip("_snippyOut/snps.bam")
        depth.columns = ['chrom','pos','depth']
        mean_depth = depth['depth'].mean()
        #genome_mean_depth = (depth['depth'] / mean_depth).mean()
        rfbXdepth = depth.loc[
            (depth['pos'] >= 869620) & (depth['pos'] <= 870919),'depth'
        ].mean()
        rfbVdepth = depth.loc[
            (depth['pos'] >= 870920) & (depth['pos'] <= 871922),'depth'
        ].mean()
        rfbUdepth = depth.loc[
            (depth['pos'] >= 872496) & (depth['pos'] <= 873558),'depth'
        ].mean()
        rfbEdepth = depth.loc[
            (depth['pos'] >= 868518) & (depth['pos'] <= 869550),'depth'
        ].mean()
        rfbNdepth = depth.loc[
            (depth['pos'] >= 873558) & (depth['pos'] <= 874503),'depth'
        ].mean()
        cyaAdepth = depth.loc[
            (depth['pos'] >= 3920255) & (depth['pos'] <= 3922801),'depth'
        ].mean()
        depth.drop(columns="chrom",inplace=True)
        depth['SampleID'] = sampleid
        depth['NormDepth'] = depth['depth'] / mean_depth
        rawlist.append(
            depth.loc[
                (depth['pos'] >= 860063) & (depth['pos'] <= 884690)
            ]
        )
        dflist.append(
            pd.DataFrame({
                'SampleID':[sampleid],
                'MeanDepth':[mean_depth],
                'rfbXdepth':[rfbXdepth],
                'rfbVdepth':[rfbVdepth],
                'rfbUdepth':[rfbUdepth],
                'rfbEdepth':[rfbEdepth],
                'rfbNdepth':[rfbNdepth],
                'cyaAdepth':[rfbNdepth]
            })
        )
    return dflist, rawlist


@app.cell
def _(dflist, pd, rawlist):
    alldepth = pd.concat([x for x in dflist],ignore_index=True)
    rawdepth = pd.concat([x for x in rawlist],ignore_index=True)
    return (rawdepth,)


@app.cell
def _(rawdepth):
    rawdepth.to_csv("../rfbLocusSeqDepth.csv",header=True,index=False)
    return


@app.cell
def _():
    #check how many have 0 read depth across genomes

    return


if __name__ == "__main__":
    app.run()
