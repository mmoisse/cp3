#!/bin/bash

files=$path/$input1/*

res_name=GBA_spans_def



cd $path/programs/GBA/

counter=0
for f in $files;
	do
	let counter=counter+1
	name="$(head -1 $f | cut -d ' ' -f 1 | tr -d '>')"
	name2="$(head -1 $f | cut -d ' ' -f 2)"
	file_name="$(basename $f)"
	echo ">"$name $name2 >> $res_name
	cp $f $path/programs/GBA/file_store/
	java applications.Gbm file_store/$file_name/ swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095 3 15 5  knowledge/blosum62Matrix 0 2>/dev/null | tail -n +2 | sed -e 's/\s\+/\n/g' >> $res_name
	echo '#' >> $res_name
	if ! (($counter % 1000)); then
		echo $counter $(date) 
	fi
done

mv $res_name $path/result_folder
rm file_store/*

echo "Process finished:" $(date)
