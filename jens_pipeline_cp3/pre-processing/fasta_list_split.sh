#!/bin/bash

input1=$1

cat $input1 | cut -d " " -f 1,2 > cut_sequences.fasta

./bioperl_script.pl cut_sequences.fasta sub_fasta 1

rm sub_fasta.1

seq_number="$(grep -e'>' $input1 | wc -l)"

upper_seq=$(expr $seq_number + 1)

counter1=1
counter2=2
while [ $counter2 -lt $upper_seq ]; 
	do
		mv sub_fasta.$counter2 sub_fasta.$counter1
		let counter1=counter1+1
		let counter2=counter2+1
	done

mkdir fasta_split_results
mv sub_fasta* fasta_split_results
cd fasta_split_results

files=sub_fasta*

for f in $files;
	do
		sed -i 's/\*//g' $f
		sed -i '/^\s*$/d' $f
		file_name="$(head -1 $f| cut -d ' ' -f 1 | tr -d '>')"
		mv $f $file_name.fasta
	done

cd ..


echo "Completed processing!"
