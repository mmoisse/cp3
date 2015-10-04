#!/bin/bash

files=$4/$1/*

cd $4/programs/VSL2/

for f in $files;
	do
		name="$(head -1 $f | cut -d ' ' -f 1 | tr -d '>')"
		name2="$(head -1 $f | cut -d ' ' -f 2)"
		echo $name $name2: > vsl2_$2.$name
		tail -n +2 $f > temp1_vsl2.fasta
		java -jar VSL2.jar -s:temp1_vsl2.fasta | python disconversion_vsl2.py >> vsl2_$2.$name
	done

rm temp*

mv vsl2_$2.* $4/$3/vsl2_files/

echo VSL2b completed, $(date)
