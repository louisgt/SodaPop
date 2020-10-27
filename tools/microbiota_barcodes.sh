for file in barcodes/*.unique; do cat $file >> barcodes/all_unique.txt; done
sort barcodes/all_unique.txt > barcodes/sorted_unique.txt
cut -f1 -d" " barcodes/sorted_unique.txt | uniq > barcodes/start.txt

echo Combining time series of all barcodes... 

i=0
j=1

cat barcodes/start.txt > barcodes/series$i.txt

#### JOIN TIME FRAMES IN A SUITABLE FORMAT
for filename in `find barcodes/ -maxdepth 1 -name "*.unique" | sort`
do
        join -t' ' -e 0 -a 1 -1 1 -2 1 -o 2.2 barcodes/series$i.txt $filename | paste -d' ' barcodes/series$i.txt - > barcodes/series$j.txt
        rm -f barcodes/series$i.txt
        ((i++))
        ((j++))
done

cat barcodes/series$i.txt | cut -d " " -f 1,3- > ALL_generations.txt



### paste with seed table
join -t' ' -e 0 -a 1 -1 1 -2 1 -o 2.2 barcodes/series_seed.txt barcodes/start.txt | paste -d' ' barcodes/series_seed.txt - > barcodes/series1.txt
