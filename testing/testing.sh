#!/bin/bash

cd ..

#./main.R -n testseq1 -f testing/test_sequences/TraesCS1A01G007100.1.fa &> testing/test_output/TraesCS1A01G007100.1.fa.output

for i in testing/test_sequences/*
do
	#echo $i
	QUERY=$i
	querynopath=${QUERY##*/}
	query2=${querynopath%.*}
	#echo $query2
	./main.R -n $query2 -f $i &> testing/test_output/$query2.output.txt
done
