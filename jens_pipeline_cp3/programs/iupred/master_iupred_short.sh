#!/bin/bash

files=$4/$1/*

cd $4/programs/iupred/

for f in $files;
	do
		name="$(head -1 $f | cut -d ' ' -f 1 | tr -d '>')"
		name2="$(head -1 $f | cut -d ' ' -f 2)"
		echo $name $name2: > iup_$2.$name
		./iupred $f short | python disconversion_iupred.py >> iup_$2.$name
	done

mv iup_$2.* $4/$3/iupred_short_files/

echo IUPred short completed, $(date)
