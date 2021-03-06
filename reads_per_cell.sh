#Extract cell barcodes from raw R1 reads
zcat raw_R1.fastq.gz | awk '{if(NR%4==2) print substr($1,1,16)}' > bc_raw.txt

#Creat a filtered cell barcode file barcodes.tsv, and search the filtered barcodes in raw data
awk 'NR==FNR{Arr[$1]++;next} ($1 in Arr){print $1}' barcodes.tsv bc_raw.txt > bc_match.txt

#Count occurrence of cell barcodes
sort bc_match.txt | uniq -c > bc_occurrence.txt
