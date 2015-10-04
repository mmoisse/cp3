#!/bin/bash

#Installs the programs included in the pipeline

cd programs/iupred
cc iupred.c -o iupred
cd ..


cd plaac/cli
./build_plaac.sh
cd ..



