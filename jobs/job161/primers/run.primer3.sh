#!/bin/bash

genename=$1
basepath=$2

if [ -z "$1" ]
then
	echo "Provide gene name as first argument"
else
	suffix=".output.txt"

	cd "$2""$genename"/primers/input
	
	#ls | parallel echo primer3_core --output={.}.bazza.txt {}
	ls | parallel primer3_core --output=../output/{.}.output.txt {}

	#for i in ./*.txt;  do
	#primer3_core < $i > ../output/$i$suffix
	##echo $i$suffix
	#done

fi



