#!/bin/bash

start=$(date)
files1=$3/$2/iupred_short_files/*
files2=$3/$2/espritz_x_files/*
files3=$3/$2/espritz_n_files/*
files4=$3/$2/vsl2_files/*
files5=$3/$2/papa_files/*

cd $3/$2/
echo Span result compilation starting, $start
touch iup_short_spans
(for f in $files1;
	do
		name="$(head -1 $f | tr -d :)"
		spans="$(grep -e'#' $f | tr -d '#')"
		echo ">"$name >> iup_short_spans
		while read -r line; do
			echo "$(echo $line | cut -d ':' -f 1)" >> iup_short_spans
		done <<< "$spans"
		echo "#" >> iup_short_spans
	done) &


touch espritz_x_spans
(for f in $files2;
	do
		name="$(head -1 $f | tr -d :)"
		spans="$(grep -e'#' $f | tr -d '#')"
		echo ">"$name >> espritz_x_spans
		while read -r line; do
			echo "$(echo $line | cut -d ':' -f 1)" >> espritz_x_spans
		done <<< "$spans"
		echo "#" >> espritz_x_spans
	done) &

touch espritz_n_spans
(for f in $files3;
	do
		name="$(head -1 $f | tr -d :)"
		spans="$(grep -e'#' $f | tr -d '#')"
		echo ">"$name >> espritz_n_spans
		while read -r line; do
			echo "$(echo $line | cut -d ':' -f 1)" >> espritz_n_spans
		done <<< "$spans"
		echo "#" >> espritz_n_spans
	done) &

touch vsl2_spans
(for f in $files4;
	do
		name="$(head -1 $f | tr -d :)"
		spans="$(grep -e'#' $f | tr -d '#')"
		echo ">"$name >> vsl2_spans
		while read -r line; do
			echo "$(echo $line | cut -d ':' -f 1)" >> vsl2_spans
		done <<< "$spans"
		echo "#" >> vsl2_spans
	done) &

touch papa_summary
(for f in $files5;
	do
		cat $f >> papa_summary
	done) &

		
		
wait %1 %2 %3 %4 %5

sed -i '/^$/d' iup_short_spans
sed -i '/^$/d' espritz_x_spans
sed -i '/^$/d' espritz_n_spans
sed -i '/^$/d' vsl2_spans
echo ----------------------------------
echo Span results compilation completed, $(date)
echo ----------------------------------
