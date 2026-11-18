# Nextflow pipeline for S. Paratyphi A analysis

Pipeline used to generate genomic data for publication.


To run pipeline, clone this repository `git clone`. 

Create conda environment using yaml file `conda env create -f environment.yml`. And activate `conda activate parayphiA`

Make input files for the pipeline:

```bash
cd fastq/
for f in *.fastq.gz ; do echo "${f:0:14}" ; done | sort | uniq > ../id.list 

realpath *.fastq.gz | grep "_R1_001.fastq.gz\|_R1.fastq.gz\|_1.trim.fastq.gz\|_R1.trim.fastq.gz" > ../R1.list
realpath *.fastq.gz | grep "_R2_001.fastq.gz\|_R2.fastq.gz\|_2.trim.fastq.gz\|_R2.trim.fastq.gz" > ../R2.list

paste id.list R1.list R2.list > input.list
```
Assumes that the fastq files have some common ending, and the first 14 characteris are the unique IDs.


Initial nextflow run:

```bash
cd Nextflow
conda activate paratyphiA
nextflow run main.nf --input_file ../input.list -c nextflow.config -resume --long_read_assemblies ../LongInput.list
```

This will generate the basic genomic characterisation output files in `results/`. Collect the results and apply QC filters. This was done in Python / Marimo notebook: `ParatyphiA_QCChecksAndCollation.py`.

Once the data has been QC checked, can run the pipeline a second time providing the QC pass list to build a core genome alignment for phylogenetic tree constructions (`--tree_samples_list`). Can also provide long read genome assemblies for characterisation as well (`--long_read_assemblies`).

Provide reference genome with `--snippy_ref` and phage regions for masking with `--snippy_mask`. 

```bash
cd Nextflow
nextflow run main.nf --input_file ../input.list -c nextflow.config -resume --long_read_assemblies ../LongInput.list --tree_samples_list ../SampleIDsForTreeAndPangenome.list --snippy_ref assets/GCF_000026565.1_ASM2656v1_genomic.gbff --snippy_mask assets/PHASTEST_AKU12601_PhageRegions.bed
```

Outputs will be in `results/`.