# HISAT2-pipeline [![DOI](https://zenodo.org/badge/475739892.svg)](https://zenodo.org/badge/latestdoi/475739892)
This pipeline runs a RNA-seq analysis automatically, using HISAT2/Stringtie with default settings. The pipeline should be particularly useful for people who lack the necessary skills to script their own pipeline or people who run RNA-seqs frequently. Please note, that it will not trim your reads, as modern mappers such as HISAT2 can handle untrimmed reads very well in almost all cases. If you wish to trim your reads anyways, please do so before running this pipeline.

## 1. Installation:
The preferred way to install the pipeline is via (bio)conda, as this will also install all required dependencies. Simply run the following command to install the pipeline:

```
conda install -c bioconda hisat2-pipeline
```
*If you encounter issues during the installation, try running these commands before trying the command above again:*

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## 2. Running an RNA-seq analysis

If not otherwise specified, the pipeline will expect the following file structure:
```
your_folder
 └── genome
 └── reads
``` 
The *reads* folder should contain your read data (*.fastq |.fq |.fastq.gz|.fq.gz*)<br>
The *genome* will have to contain the **genome fasta file** (can be gzipped), as well as the **corresponding gff file**.
<br><br>
Then simply run the following command to run the analysis:
```
hisat2-pipeline
```

The pipeline will ask you for feedback before proceeding to the actual analysis, to confirm whether it was able to infer all settings correctly.
<br>
Upon completion, you will find a folder called *./mapping* containing the mapped \*.bam files, as well as a a **mapping_summary** containing the merged mapping summary for all files as tab-separated and Excel file.

Gene expression values are found in the folder *./abundance*, as well as the merged_FPKM and merged_TPM files (both as tab-separated and Excel files)

## 3. Advanced parameters
You can manually specify the reads, genome and output folders via:
```
hisat2-pipeline --reads_folder [YOUR_READS_FOLDER] --genome_folder [YOUR_GENOME_FOLDER] --outfolder [YOUR_OUTPUT_FOLDER]
```
This pipeline uses the default parameters for HISAT2 and Stringtie, but additional options can be specified. The example below shows how:
```
hisat2-pipeline --hisat_options "--very-sensitive --no-spliced-alignment" --stringtie_options "-m 150 -t"
```
To see all run parameters, use:
```
hisat2-pipeline -h
```

## 4. Citing

If you use this software for your scientific research, please cite it.
Refer to the **about** section for informations on how to cite this software or see the **CITATION.cff** file.

