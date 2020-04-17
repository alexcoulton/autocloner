#!/bin/bash
if [ -z "$1" ]
then
	echo "Provide gene name as first argument"
else
	echo 
	echo "####################################"	
	echo "Creating folder structure..."	
	echo "####################################"	
	echo 
	mkdir -p jobs/"$1"	
	mkdir -p jobs/"$1"/blast.results
	mkdir -p jobs/"$1"/primers	
	mkdir -p jobs/"$1"/primers/best.primers
	mkdir -p jobs/"$1"/primers/input
	mkdir -p jobs/"$1"/primers/output
	mkdir -p jobs/"$1"/primers/sequencing.primers
	mkdir -p jobs/"$1"/seq
	mkdir -p jobs/"$1"/seq/extended
	mkdir -p jobs/"$1"/seq/extended/alignments
	mkdir -p jobs/"$1"/seq/extended/seqs
	mkdir -p jobs/"$1"/seq/extended/primers.set1
	mkdir -p jobs/"$1"/seq/extended/primers.set2


fi

