# Single cell RNA-Seq barcode error correction

## User Guide

#### Get help information
#### `bash correct_barcode.sh -h`
#### ./correct_barcode.sh [options]
#### -h --help
#### Please specify the following options:
#### -fq1 --fastq_cell=SRR1853178_1.fastq   [gzipped or plain fastq R1 cell index reads]
#### -fq2 --fastq_biology=SRR1853178_2.fastq   [gzipped or plain fastq R2 biology reads]
#### -bcs --bc_start=1   [start position of cell barcode in cell index reads]
#### -bcl --bc_length=12   [length of cell barcodes in cell index reads, 12 for dropseq data, 16 for 10XGv2,v3 data]
#### -t --threads=8 [number of threads to use]
#### -o --output=BCFIX   [output file name]
#### -w --whitelist=NO [whitelist of barcodes, optional for 10XG data]

#### An example to run the program on a 10XGv2 dataset
#### `bash barcode_correct4.sh -fq1=hgmm_6k_R1.fastq.gz -fq2=hgmm_6k_R1.fastq.gz -bcs=1 -bcl=16 -t=8 -o=BCFIX -w=10Xv2_whitelist.txt`
