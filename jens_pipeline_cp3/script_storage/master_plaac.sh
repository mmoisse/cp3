#!/bin/bash

cd $3/programs/plaac/cli/target/

files=$3/$1/*

for f in $files;
	do
		cat $f >> plaac_read.fasta
	done


java -jar plaac.jar -i plaac_read.fasta -a 0 > plaac_output

sed '1,9d' -i plaac_output

rm plaac_read.fasta
mv plaac_output $3/$2/
