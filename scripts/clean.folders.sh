#!/bin/bash

genename=$1
basepath=$2

if [ -z "$1" ]
then
	echo "Provide gene name as first argument"
else
	echo "Cleaning directories"	
	rm -f "$2""$genename"/primers/input/*
	rm -f "$2""$genename"/primers/output/*
	rm -f -r "$2""$genename"/primers/best.primers/*
	rm -f -r "$2""$genename"/primers/penalties/*
	echo
fi



