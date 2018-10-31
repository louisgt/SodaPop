#!/bin/bash

##### PARSE INPUT FROM C++ MAIN OR CMD LINE
OUT=$1
MAXGEN=$2
POPSIZE=$3
DT=$4
ENCODING=$5
### Gene count *per cell*, not for the simulation as a whole
GENE_COUNT=$6
COL=3

##### CREATE DIRECTORIES FOR RESULTS
rm -rf out/$OUT/barcodes; mkdir out/$OUT/barcodes
rm -rf out/$OUT/graph; mkdir out/$OUT/graph
rm -rf out/$OUT/gene; mkdir out/$OUT/gene
PREFIX=out/$OUT
HOME=../..

echo Begin analysis.

export LC_NUMERIC=C

FACTOR=2
let "FACTOR += 2*$GENE_COUNT"

### implement: detect file encoding

if [ "$ENCODING" -eq "0" ]; then
	echo Using long output format.
	let "COL = 3"
else if [ "$ENCODING" -eq "1" ]; then
		echo Using short output format.
		let "COL = 4"
	 else 
	 	echo Using DNA output format.
	 	let "COL = 3"
	 fi
fi

echo Gene count = $6.
echo Working in $PREFIX.
cd $PREFIX

echo Gunzipping snapshots and converting to .txt...
## CONVERT BINARY SNAPSHOTS TO TEXT FILES
FILES=snapshots/*.snap.gz
for filename in $FILES
do
	gunzip -k $filename
	y="$(basename $filename .gz)"
	./$HOME/sodasnap snapshots/$y snapshots/$y.txt
	rm snapshots/$y
	gzip -f snapshots/$y.txt
done

echo Extracting barcodes from $PREFIX/snapshots/...
#### EXTRACT AND SORT BARCODES
rm -f avg_fitness.txt
FILES=snapshots/*.snap.txt.gz
for filename in $FILES
do
	y=${filename%%.txt.gz}
	if [ "$ENCODING" -eq "1" ]; then
		gunzip -c $filename | awk 'NR>3 {print $0}' - | cut -f1 | sort > barcodes/${y##*/}.barcodes
		#### SUM POPULATION FITNESS FOR EACH TIME POINT AND DIVIDE BY POP SIZE
		gunzip -c $filename | awk 'NR>3 {print $0}' - | awk -v N=$3 -v C=$COL '{sum += $C} END {printf "%.9f\n",sum/N}' - >> avg_fitness.txt

		#### FOR EACH GENE, CALCULATE MEAN STABILITY
		for i in $( seq 0 $((GENE_COUNT-1)) )
		do
			gunzip -c $filename | grep -w ^"$i" | awk -F ' ' -v N=$3 '{sum += $4} END {printf "%.9f\n",sum/N}' >> gene/stabil_$i.txt
		done
	else
		gunzip -c $filename | awk 'NR>2 {print $0}' - | awk -v N=$FACTOR 'NR%N==2' - | cut -f1 | sort > barcodes/${y##*/}.barcodes
		#### SUM POPULATION FITNESS FOR EACH TIME POINT AND DIVIDE BY POP SIZE
		gunzip -c $filename | awk 'NR>2 {print $0}' - | awk -v N=$FACTOR 'NR%N==2' - | awk -v N=$3 -v C=$COL '{sum += $C} END {printf "%.9f\n",sum/N}' - >> avg_fitness.txt

		#### FOR EACH GENE, CALCULATE MEAN STABILITY
		for i in $( seq 0 $((GENE_COUNT-1)) )
		do
			gunzip -c $filename | grep -w ^"$i" | awk -F ' ' -v N=$3 '{sum += $4} END {printf "%.9f\n",sum/N}' >> gene/stabil_$i.txt
		done
	fi
done

echo Parsing unique barcodes...

#### PARSE UNIQUE BARCODES
FILES=barcodes/*.barcodes
for filename in $FILES
do
	y=${filename%%.barcodes}
	uniq -c $filename | awk -F' ' '{t = $1; $1 = $2; $2 = t; print; }' > barcodes/${y##*/}.unique

	#### BREAK AT FIXATION POINT IF IT OCCURS
	if ! $(read -r && read -r)
	then
		#### OUTPUT FIXATION GENERATION TO FILE
	  	echo Fixation at ${y##*/} > fixation.txt
	  	echo Fixation at ${y##*/}
	  	break
	fi < barcodes/${y##*/}.unique
done

cat barcodes/$OUT.gen0000000001.snap.unique > barcodes/start.txt

echo Combining time series...

i=0
j=1

cat barcodes/start.txt > barcodes/series$i.txt

#### JOIN TIME FRAMES IN A SUITABLE FORMAT
for filename in `find barcodes/ -maxdepth 1 -name "*.unique"`
do
	join -t' ' -e 0 -a 1 -1 1 -2 1 -o 2.2 barcodes/series$i.txt $filename | paste -d' ' barcodes/series$i.txt - > barcodes/series$j.txt
	rm -f barcodes/series$i.txt
	((i++))
	((j++))
done

cat barcodes/series$i.txt | cut -d " " -f 1,3- > ALL_generations.txt

rm -f barcodes/series*.txt

#remove compressed binary snapshots
#rm -f snapshots/*.snap.gz

cd $HOME

#### PLOT RESULTS IN R SCRIPTS
Rscript tools/polyclonal_structure.R /out/$OUT/ $DT

echo Done.