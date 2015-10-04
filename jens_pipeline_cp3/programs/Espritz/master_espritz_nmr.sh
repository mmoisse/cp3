#!/bin/bash

cd $4/programs/Espritz/nmr_temp/

files=$4/programs/Espritz/nmr_run/*

nvar="_n"

for f in $files;
	do
		name="$(head -1 $f | cut -d ' ' -f 1 | tr -d '>')"
		name2="$(head -1 $f | cut -d ' ' -f 2)"
		echo $name $name2: > espr_$2$nvar.$name
	done

cd ..

cp disconversion_espr_nmr.py nmr_temp/

./espritz.pl nmr_run/ N 0 >/dev/null &

wait %1

cd nmr_run

rm *.fasta

cd ../nmr_temp/

giles=$4/programs/Espritz/nmr_run/*
for g in $giles;
	do
		name="$(basename $g | sed 's/.espritz//g')"
		cat $g | python disconversion_espr_nmr.py >> espr_$2$nvar.$name
	done

mv espr_$2$nvar.* $4/$3/espritz_n_files/
rm disconversion_espr_nmr.py
rm $4/programs/Espritz/nmr_run/*

echo Espritz nmr complete, $(date)
