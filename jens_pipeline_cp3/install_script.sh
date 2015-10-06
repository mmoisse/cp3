#!/bin/bash

#Installs the programs included in the pipeline

cd programs/iupred
cc iupred.c -o iupred
cd ../..


programs/plaac/cli/./build_plaac.sh
cp script_storage/master_plaac.sh programs/plaac/cli/target/
