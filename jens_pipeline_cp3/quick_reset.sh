#!/bin/bash

if [ -z "$1" ];
	then echo "Input Error:"
	     echo "Please provide the latest used tag name"
	     exit 0
fi

tag=$1

rm programs/iupred/iup_$tag.*
rm programs/Espritz/nmr_run/*
rm programs/Espritz/nmr_temp/*
rm programs/Espritz/xray_run/*
rm programs/Espritz/xray_temp/*
rm programs/papa/papa_$tag.*
rm programs/GBA/file_store/*
rm programs/VSL2/vsl2_$tag.*

