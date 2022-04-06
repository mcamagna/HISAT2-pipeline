# HISAT2-pipeline [![DOI](https://zenodo.org/badge/475739892.svg)](https://zenodo.org/badge/latestdoi/475739892)
The following document will give you an explanation of how to use the pipeline. If you encounter problems or would like me to make changes to the pipeline, feel free to contact me.

## 1. Prerequisites:

The script will work on Linux and MacOS. In order to run the program, the following programs/libraries need to be installed on your system:

*Python (any version ≥ 3)\
HISAT2\
stringtie\
samtools\
pandas (optional, but strongly recommended to get result reports)\
openpyxl (optional, allows excel output, otherwise, output will be text files only)*

Using conda / anaconda you should be able to install these programs using the command:

```
conda install hisat2 stringtie samtools pandas openpyxl
```

*If conda cannot install the packages, try to run these commands before trying the command above again:*

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## 2. Running the program

In order to run the script, you should prepare your data in the following way:

any_folder/

* ../Pipeline.py
* ../reads
* ../genome

Copy the Pipeline.py into a folder of your choice. Inside this folder, make two new folders “reads” and “genome”. Copy your read data (*.fastq |*.fq | *.fastq.gz |* .fq.gz) into the reads folder. Copy the genome fasta file into the genome folder, as well as the gff/gtf file, which contains the gene locations for that genome.

Open a terminal or a console and go to the folder containing the script. Run it using:

```
python Pipeline.py
```

The program will first check whether all required software is installed and notify you if it cannot find a software in question. It will then analyse your reads and determine whether they are split into different files, and whether they are paired or unpaired. You’ll be asked to confirm the detection of read types.

If it cannot find a genome index in the genome folder, it will create one, otherwise it will skip building an index.

Once mapping is complete, you will find a folder called “mapping” containing the mapped \*.sam files, as well as a “**mapping_summary.tsv**” containing the merged mapping summary for all files. If pandas and openpyxl were installed, you’ll also find an Excel file.

Once stringtie has finished, you’ll find the gene expression values in the folder “abundance”, as well as a file called “**merged_FPKM.tsv**” (and potentially an Excel file)

## 3. Additional parameters

Run python Pipeline.py -h to see additional commands


## 4. Citing

If you use this software for your scientific research, please cite it.
Refer to the **about** section for informations on how to cite this software.
