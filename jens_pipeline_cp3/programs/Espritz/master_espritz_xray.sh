#!/bin/bash

cd $4/programs/Espritz/xray_temp/

files=$4/programs/Espritz/xray_run/*

xvar="_x"

for f in $files;
	do
		name="$(head -1 $f | cut -d ' ' -f 1 | tr -d '>')"
		name2="$(head -1 $f | cut -d ' ' -f 2)"
		echo $name $name2: > espr_$2$xvar.$name
	done

cd ..

cp disconversion_espr_xray.py xray_temp/

./espritz.pl xray_run/ N 0 >/dev/null &

wait %1

cd xray_run

rm *.fasta

cd ../xray_temp/

giles=$4/programs/Espritz/xray_run/*
for g in $giles;
	do
		name="$(basename $g | sed 's/.espritz//g')"
		cat $g | python disconversion_espr_xray.py >> espr_$2$xvar.$name
	done

mv espr_$2$xvar.* $4/$3/espritz_x_files/
rm disconversion_espr_xray.py
rm $4/programs/Espritz/xray_run/*

echo Espritz xray complete, $(date)
