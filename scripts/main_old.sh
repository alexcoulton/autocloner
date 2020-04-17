#!/bin/bash

if [ -z "$1" ] # || [ -z "$2" ] #|| [ -z "$3" ] || [ -z "$4" ]
then
	echo "Arguments: gene.name, fasta.path" #path.to.multiple.alignment, product.size, full.gene"
else
	genename=$1
	fastapath=$2
	# multalignpath=$2
	# productsize=$3
	# fullgene=$4

	./create.folder.structure.sh $genename
	cp ./clean.folders.sh ./"$genename"/primers/
	cp ./run.primer3.sh ./"$genename"/primers/

	./make.fasta.indexes.rscript.R
	echo 
	echo "####################################"
	echo "Running perform.blast.rscript.R"
	echo "####################################"
	echo 
	./perform.blast.rscript.R $genename $fastapath
	echo 
	echo "####################################"
	echo "Running blast.scaffold.parser.rscript.R"
	echo "####################################"
	echo 
	./blast.scaffold.parser.rscript.R $genename $fastapath
	echo 
	echo "####################################"
	echo "Running muscle to generate alignment"
	echo "####################################"
	echo 
	muscle -in "$genename"/seq/extended/seqs/all.fa -out "$genename"/seq/extended/alignments/all.align.rev.fa

	./rearrange.muscle.rscript.R $genename

	echo 
	echo "####################################"
	echo "Running primer.selection.script.R"
	echo "####################################"
	echo 
	./primer.selection.rscript.R $genename
	echo 
	echo "####################################"
	echo "Running primer.evaluation.rscript.R"
	echo "####################################"
	echo 
	./primer.evaluation.rscript.R $genename
fi