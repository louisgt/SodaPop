#!/bin/bash

##### PARSE INPUT FROM C++ MAIN OR CMD LINE
OUT=$1
MAXGEN=$2
POPSIZE=$3
DT=$4
LONG=$5

##### CREATE DIRECTORIES FOR RESULTS
rm -r out/$OUT/barcodes; mkdir out/$OUT/barcodes
rm -r out/$OUT/graph; mkdir out/$OUT/graph
PREFIX=out/$OUT
echo $PREFIX

echo Extracting barcodes from $PREFIX/snapshots/...

#### CONVERT BINARY SNAPSHOTS TO TEXT FILES
FILES=$PREFIX/snapshots/*.snap
for filename in $FILES
do
	y=${filename%.001}
	./sodasnap $filename $y.txt $LONG
done

#### EXTRACT AND SORT BARCODES
rm $PREFIX/avg_fitness.txt
FILES=$PREFIX/snapshots/*.snap.txt
for filename in $FILES
do
	y=${filename%%.txt}
	grep -w "C" $filename | cut -f1 | sort > $PREFIX/barcodes/${y##*/}.barcodes
	#### SUM POPULATION FITNESS FOR EACH TIME POINT AND DIVIDE BY POP SIZE
	grep -w "C" $filename | awk -v N=$3 '{sum += $5} END {print sum/N}' >> $PREFIX/avg_fitness.txt
done

echo Parsing unique barcodes...

#### PARSE UNIQUE BARCODES
FILES=$PREFIX/barcodes/*.barcodes
for filename in $FILES
do
	y=${filename%%.barcodes}
	uniq -c $filename | sed "s/^[ \t]*//" | awk -F' ' '{t = $1; $1 = $2; $2 = t; print; }' > $PREFIX/barcodes/${y##*/}.unique

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
####for filename in `ls -v ./$PREFIX/barcodes/*.unique`
for filename in `find $PREFIX/barcodes/ -maxdepth 1 -name "*.unique"`
do
	join -t' ' -e 0 -a 1 -1 1 -2 1 -o 2.2 $PREFIX/barcodes/series$i.txt $filename | paste -d' ' $PREFIX/barcodes/series$i.txt - > $PREFIX/barcodes/series$j.txt
	rm $PREFIX/barcodes/series$i.txt
	((i++))
	((j++))
done

cat $PREFIX/barcodes/series$i.txt | cut -d " " -f 1,3- > $PREFIX/ALL_generations.txt

rm $PREFIX/barcodes/series*.txt

#### PLOT RESULTS IN R SCRIPTS
Rscript tools/polyclonal_structure.R /out/$OUT/ $DT

echo Done.