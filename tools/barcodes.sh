#!/bin/bash

##### PARSE INPUT FROM C++ MAIN OR CMD LINE
OUT=$1
MAXGEN=$2
POPSIZE=$3
DT=$4
FORMAT=$5

LONG=0

##### CREATE DIRECTORIES FOR RESULTS
rm -rf out/$OUT/barcodes; mkdir out/$OUT/barcodes
rm -rf out/$OUT/graph; mkdir out/$OUT/graph
PREFIX=out/$OUT

echo Begin analysis.

if [ "$FORMAT" -eq "$LONG" ]; then
	echo Using long format.
else
	echo Using short format.
fi

echo Working in $PREFIX.

echo Extracting barcodes from $PREFIX/snapshots/...

### CONVERT BINARY SNAPSHOTS TO TEXT FILES
FILES=$PREFIX/snapshots/*.snap
for filename in $FILES
do
	y=${filename%.001}
	./sodasnap $filename $y.txt $FORMAT
done

#### EXTRACT AND SORT BARCODES
rm -f $PREFIX/avg_fitness.txt
FILES=$PREFIX/snapshots/*.snap.txt
for filename in $FILES
do
	y=${filename%%.txt}
	awk 'NR>1 {print $0}' $filename | awk 'NR%4==2' - | cut -f1 | sort > $PREFIX/barcodes/${y##*/}.barcodes
	#### SUM POPULATION FITNESS FOR EACH TIME POINT AND DIVIDE BY POP SIZE
	if [ "$FORMAT" -eq "$LONG" ]; then
		awk 'NR>1 {print $0}' $filename | awk 'NR%4==2' - | awk -v N=$3 '{sum += $3} END {print sum/N}' - >> $PREFIX/avg_fitness.txt
    else
		awk 'NR>1 {print $0}' $filename | awk 'NR%4==2' - | awk -v N=$3 '{sum += $4} END {print sum/N}' - >> $PREFIX/avg_fitness.txt
    fi
	
done

echo Parsing unique barcodes...

#### PARSE UNIQUE BARCODES
FILES=$PREFIX/barcodes/*.barcodes
for filename in $FILES
do
	y=${filename%%.barcodes}
	uniq -c $filename | awk -F' ' '{t = $1; $1 = $2; $2 = t; print; }' > $PREFIX/barcodes/${y##*/}.unique

	#### BREAK AT FIXATION POINT IF IT OCCURS
	if ! $(read -r && read -r)
	then
		#### OUTPUT FIXATION GENERATION TO FILE
	  	echo Fixation at ${y##*/} > $PREFIX/fixation.txt
	  	echo Fixation at ${y##*/}
	  	break
	fi < $PREFIX/barcodes/${y##*/}.unique
done

cat $PREFIX/barcodes/$OUT.gen0000000001.snap.unique > $PREFIX/barcodes/start.txt

echo Combining time series...

i=0
j=1

cat $PREFIX/barcodes/start.txt > $PREFIX/barcodes/series$i.txt

#### JOIN TIME FRAMES IN A SUITABLE FORMAT
for filename in `find $PREFIX/barcodes/ -maxdepth 1 -name "*.unique"`
do
	join -t' ' -e 0 -a 1 -1 1 -2 1 -o 2.2 $PREFIX/barcodes/series$i.txt $filename | paste -d' ' $PREFIX/barcodes/series$i.txt - > $PREFIX/barcodes/series$j.txt
	rm -f $PREFIX/barcodes/series$i.txt
	((i++))
	((j++))
done

cat $PREFIX/barcodes/series$i.txt | cut -d " " -f 1,3- > $PREFIX/ALL_generations.txt

rm -f $PREFIX/barcodes/series*.txt

rm -Rf $PREFIX/snapshots/*.snap

#### PLOT RESULTS IN R SCRIPTS
Rscript tools/polyclonal_structure.R /out/$OUT/ $DT

echo Done.