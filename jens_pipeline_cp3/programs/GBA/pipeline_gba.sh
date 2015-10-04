#!/bin/bash

start="$(date)"

result_name=$1_lcr_spans

files=$1/*
for f in $files;
	do
		header_name="$(head -1 $f)"
		file_name="$(basename $f)"
		echo $header_name >> $result_name
		java applications.Gbm $1/$file_name/ swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095 3 15 5  knowledge/blosum62Matrix 0 2>/dev/null | tail -n +2 | sed -e 's/\s\+/\n/g' >> $result_name
		echo '#' >> $result_name

	done

rm -r $1
