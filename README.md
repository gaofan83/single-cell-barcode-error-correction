# Single cell RNA-Seq barcode error correction

##User Guide

####Get help information
`bash barcode_correct4.sh -h`
./barcode_correct.sh [options]
-h --help
Please specify the following options:
-fq1 --fastq_cell=SRR1853178_1.fastq   [gzipped fastq R1 or R2 for cell index reads]
-fq2 --fastq_biology=SRR1853178_2.fastq   [gzipped fastq R1 or R2 for biology reads]
-bcs --bc_start=1   [start position of cell barcode in cell index reads]
-bcl --bc_length=12   [length of cell barcodes in cell index reads]
-t --threads=8 [number of threads to use]
-o --output=BCFIX   [output file name]
-w --whitelist=NO [whitelist of barcodes, optional]


####An example to run the program on a 10XGv2 dataset
`bash barcode_correct4.sh -fq1=hgmm_6k_R1.fastq.gz -fq2=hgmm_6k_R1.fastq.gz -bcs=1 -bcl=16 -t=8 -o=BCFIX -w=10Xv2_whitelist.txt`
