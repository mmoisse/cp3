#!/bin/bash

if [ -z "$1" ];
	then echo "Input Error:"
	     echo "Please provide the sequence folder name"
	     exit 0
fi

if [ -z "$2" ];
	then echo "Input Error:"
	     echo "Please provide tag name for result files, e.g 'sc' for S.cerevisiae"
	     exit 0
fi

if [ -z "$3" ];
	then echo "Input Error:"
	     echo "Please provide name for result folder"
             exit 0
fi

if [ -d "$3" ];
	then echo "---------------------------------"
	     echo "Output directory already exists"
	     echo "Please try another directory name"
	exit 0
fi

echo Analysis Started: $(date)
echo "---------------------------------"

start="$(date)"

path="$(pwd)"

if [[ ! -e $path/$1 ]]; then
    echo "---------------------------------"
    echo "Input Error, folder not found..."
    echo "---------------------------------"
    exit 0
fi

export IUPred_PATH=$path/programs/iupred/

mkdir -p $3/{espritz_n_files,espritz_x_files,iupred_short_files,papa_files,vsl2_files}

cp script_storage/master_compiler.py $3
cp script_storage/rerun_analysis.py $3

files=$1/*

for f in $files;
	do
		cp $f programs/Espritz/xray_run/
		cp $f programs/Espritz/nmr_run/
	done

trim_input="$(echo $1 | tr -d '/')"

sh $path/programs/iupred/master_iupred_short.sh $trim_input $2 $3 $path &
sh $path/programs/Espritz/master_espritz_xray.sh $trim_input $2 $3 $path &
sh $path/programs/Espritz/master_espritz_nmr.sh $trim_input $2 $3 $path &
sh $path/programs/VSL2/master_vsl2.sh $trim_input $2 $3 $path &
sh $path/programs/papa/master_papa.sh $trim_input $2 $3 $path &
sh $path/programs/plaac/cli/target/master_plaac.sh $trim_input $3 $path &



wait %1 %2 %3 %4 %5 %6



echo -----------------------------------------------------
echo Intrinsic disorder analysis complete!
echo Analysis Initiated: $start
echo Analysis complete: $(date)
echo -----------------------------------------------------



./span_handler.sh $1 $3 $path &

wait %1

echo Low complexity analysis and final summaries compiling, $(date)
python $3/master_compiler.py $1 $2 $3 $path &

wait %1

cd $3/
rm master_compiler.py
#espritz_x_files iupred_short_files papa_files vsl2_files
echo Full analysis complete! $(date)
 
