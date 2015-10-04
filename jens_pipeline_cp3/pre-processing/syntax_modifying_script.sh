#!/bin/bash

input=$1
input2=$2
input3=$3
path="$(pwd)"

echo -------------------------------------------------------------------
echo Converting files in folder: $input, Field for file name: $input3
echo -------------------------------------------------------------------

files=$path/$input/*

for f in $files;
	do
		syst_name="$(grep -e'>' $f | cut -d ' ' -f $input3 | tr -d '>')"
		cat $f | cut -d ' ' -f 1,2 | sed 's/\*//g' > $syst_name.fasta
	done

mkdir $input2

mv *.fasta $input2
