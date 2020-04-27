#!/bin/bash

genename=$1
basepath=$2

if [ -z "$1" ]
then
	echo "Provide gene name as first argument"
else
	suffix=".output.txt"

	cd "$2""$genename"/primers/input
	
	#parallel solution
	#ls | parallel primer3_core --output=../output/{.}.output.txt {}
	

	for i in ./*.txt
	do
		primer3_core < $i > ../output/$i$suffix	
	done

fi



