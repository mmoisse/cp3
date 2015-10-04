#!/bin/bash

files=$4/$1/*

cd $4/programs/papa/

for f in $files;
	do
		name="$(head -1 $f | cut -d ' ' -f 1 | tr -d '>')"
		python papa.py --ignore_fold_index -o papa_$2.$name $f
	done

mv papa_$2.* $4/$3/papa_files/

echo papa.py prion analysis completed, $(date)

