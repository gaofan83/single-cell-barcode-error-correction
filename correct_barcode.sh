#!/bin/bash
export LANG=C
export LC_ALL=C

FILE_SIZE=100000000
THREAD=8
CELL=50000
SAMPLING="YES"
READ_CELL_INDEX="SRR1853178_1.fastq"
READ_BIOLOGY="SRR1853178_2.fastq"
BC_START=1
BC_LENGTH=12
OUTPUT="BCFIX"
WHITELIST="NO"

function usage()
{
    echo "./correct_barcode.sh [options]"
    echo "-h --help"
    echo "Please specify the following options:"
    echo "-fq1 --fastq_cell=$READ_CELL_INDEX   [gzipped or plain fastq R1 cell index reads]"
    echo "-fq2 --fastq_biology=$READ_BIOLOGY   [gzipped or plain fastq R2 biology reads]"
    echo "-bcs --bc_start=$BC_START   [start position of cell barcode in cell index reads]"
    echo "-bcl --bc_length=$BC_LENGTH   [length of cell barcodes in cell index reads, 12 for dropseq data, 16 for 10XGv2,v3 data]"
    echo "-t --threads=$THREAD [number of threads to use]"
    echo "-o --output=$OUTPUT   [output file name]"
    echo "-w --whitelist=$WHITELIST [whitelist of barcodes, optional for 10XG data]"
    echo ""
    echo "Example below:"
    echo "bash correct_barcode.sh -fq1=SRR1853178_1.fastq -fq2=SRR1853178_2.fastq -bcs=1 -bcl=12 -t=8 -o=BCFIX"
}

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        -fq1 |--fastq_cell)
            READ_CELL_INDEX=$VALUE
            ;;
        -fq2 |--fastq_biology)
            READ_BIOLOGY=$VALUE
            ;;
        -bcs | --bc_start)
            BC_START=$VALUE
            ;;
        -bcl | --bc_length)
            BC_LENGTH=$VALUE
            ;;
        -t | --threads)
            THREAD=$VALUE
            ;;
        -o |--output)
            OUTPUT=$VALUE
            ;;
        -w |--whitelist)
            WHITELIST=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done

echo "READ_CELL_INDEX is $READ_CELL_INDEX";
echo "READ_BIOLOGY is $READ_BIOLOGY";
echo "BC_START is $BC_START";
echo "BC_LENGTH is $BC_LENGTH";
echo "THREAD is $THREAD";
echo "OUTPUT is $OUTPUT";
echo "WHITELIST is $WHITELIST";

N=$THREAD

#split the raw fastq file
rm -rf temp
rm -rf temp2
mkdir temp
mkdir temp2
READ_CELL_LIST=${READ_CELL_INDEX//[,]/ }
READ_BIOLOGY_LIST=${READ_BIOLOGY//[,]/ }

if [[ $READ_CELL_LIST == *.gz ]]; then
  zcat $READ_CELL_LIST | split -dl $FILE_SIZE - temp/ &
else
  cat $READ_CELL_LIST | split -dl $FILE_SIZE - temp/ &
fi

if [[ $READ_BIOLOGY_LIST == *.gz ]]; then
  zcat $READ_BIOLOGY_LIST | split -dl $FILE_SIZE - temp2/ &
else
  cat $READ_BIOLOGY_LIST | split -dl $FILE_SIZE - temp2/ &
fi

wait
ls temp > temp.ID
find temp/ -type f ! -name "*.*" -exec mv "{}" "{}".fq1 \;
find temp2/ -type f ! -name "*.*" -exec mv "{}" "{}".fq2 \;
mv temp2/* temp/
rm -r temp2

#extract barcodes 
  run_bc(){
        ID=$1
        awk -v start="$BC_START" -v len="$BC_LENGTH" -F\\t '{if(NR%4==2) printf "%s\t",substr($0,start,len); if(NR%4==0) printf "%s\n",substr($0,start,len);}' temp/$ID.fq1 > temp/$ID.bc
        sort -k1,1 temp/$ID.bc | awk '{print $1}' > temp/$ID.bc_sort
          }

  while read line;
  do
    ((i=i%N)); ((i++==0)) && wait
    run_bc "$line" &
  done < temp.ID
  wait
  echo "Cell barcodes were extracted..."

if [ $WHITELIST == "NO" ]; then
#get corrected whitelist 
  if [ $SAMPLING = "YES" ]; then
    uniq -c temp/00.bc_sort | sort -nr -k1,1 | awk -v cell="$CELL" '{if(FNR <= cell) print FNR" "$2}' > temp/bc.whitelist_start
  else
    sort -m temp/*bc_sort | uniq -c | sort -nr -k1,1 | awk -v cell="$CELL" '{if(FNR <= cell) print FNR" "$2}' > temp/bc.whitelist_start
  fi
  rm temp/*bc_sort
  cp temp/bc.whitelist_start temp/bc.whitelist_process

  for (( i=$BC_START; i<=$((BC_START+BC_LENGTH-1)); i++ ))
  do
    awk -v pos="$i" '{print $1" "substr($2,1,pos-1)substr($2,pos+1)}' temp/bc.whitelist_process > temp/bc.whitelist_cut
    sort -k2,2 -k1,1n temp/bc.whitelist_cut | uniq -f1 | sort -k1,1n | awk '{print $1}' > temp/bc.list_keep
    awk 'NR==FNR{i++;Arr[$1]=$2;next} ($0 in Arr){print $0" "Arr[$0]}' temp/bc.whitelist_process temp/bc.list_keep > temp/bc.whitelist_keep
    mv temp/bc.whitelist_keep temp/bc.whitelist_process
  done
  awk '{print $2}' temp/bc.whitelist_process > temp/bc.whitelist_final
  rm temp/bc.whitelist_cut
  echo "A whitelist of enriched cell barcodes were generated..."

else 
  rm temp/*bc_sort
  cp $WHITELIST temp/bc.whitelist_final
fi

#generate 1bp mismatch whitelist
for (( i=$BC_START; i<=$((BC_START+BC_LENGTH-1)); i++ ))
do
  rm -f temp/bc.whitelist_del.$i
  awk -v pos=$i -v base="$base" '{print substr($0,pos,1)"\t"substr($0,1,pos-1)substr($0,pos+1)}' temp/bc.whitelist_final >> temp/bc.whitelist_del.$i
done


#search for barcode errors
run_error(){
     ID=$1
     BC_START=$2
     BC_LENGTH=$3
     rm -f temp/$ID.index_tofix
     awk 'NR==FNR{Arr[$0]++;next} !($1 in Arr){print $0"\t"FNR}' temp/bc.whitelist_final temp/$ID.bc > temp/$ID.bc_tofix
     for (( i=$BC_START; i<=$((BC_START+BC_LENGTH-1)); i++ ))
     do
       awk -v pos="$i" 'NR==FNR{j++;BASE[$2]=$1;next} (substr($1,1,pos-1)substr($1,pos+1) in BASE){bc=substr($1,1,pos-1)substr($1,pos+1); print pos"\t"BASE[bc]"\t"substr($2,pos,1)"\t"$3}' \
             temp/bc.whitelist_del.$i temp/$ID.bc_tofix >> temp/$ID.index_tofix
     done
     
     sort -k4,4n -k3,3 temp/$ID.index_tofix > temp/$ID.index_tofix_sort
#    uniq -f3 -u temp/$ID.index_tofix_sort > temp/$ID.index_tofix_uniq
     uniq -f2 -c temp/$ID.index_tofix_sort > temp/$ID.index_tofix_count
     uniq -f4 temp/$ID.index_tofix_count | awk '{if($1==1) print $2"\t"$3"\t"$4"\t"$5}' > temp/$ID.index_tofix_uniq

     awk -F\\t 'NR==FNR{j++;POS[$4]=$1;BASE[$4]=$2;next} (int((FNR+3)/4) in POS){if(FNR%4==2) {print substr($0,1,POS[(FNR+2)/4]-1)BASE[(FNR+2)/4]substr($0,POS[(FNR+2)/4]+1)} else {print $0}}' \
             temp/$ID.index_tofix_uniq temp/$ID.fq1 > temp/$ID.fq1_fixed
     awk -F\\t 'NR==FNR{Arr[$4]++;next} !(int((FNR+3)/4) in Arr){print $0}' temp/$ID.index_tofix_uniq temp/$ID.fq1 > temp/$ID.fq1_unfixed
     awk -F\\t 'NR==FNR{Arr[$4]++;next} (int((FNR+3)/4) in Arr){print $0}' temp/$ID.index_tofix_uniq temp/$ID.fq2 > temp/$ID.fq2_fixed
     awk -F\\t 'NR==FNR{Arr[$4]++;next} !(int((FNR+3)/4) in Arr){print $0}' temp/$ID.index_tofix_uniq temp/$ID.fq2 > temp/$ID.fq2_unfixed
           }

while read line;
do
  ((i=i%N)); ((i++==0)) && wait
  run_error "$line" "$BC_START" "$BC_LENGTH" &
done < temp.ID
wait

cat temp/*fq1_unfixed temp/*fq1_fixed > ${OUTPUT}_1.fastq &
#cat temp/*fq1_fixed > SRR1853178fixed_1.fastq &
cat temp/*fq2_unfixed temp/*fq2_fixed > ${OUTPUT}_2.fastq &
#cat temp/*fq2_fixed > SRR1853178fixed_2.fastq &
wait
echo "Cell barcode errors have been fixed. Please use "$OUTPUT"_1.fastq and "$OUTPUT"_2.fastq for downstream processing."
