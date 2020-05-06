# AutoCloner
### Gene cloning in polyploids made easy
Initial requirements

AutoCloner requires muscle, BLAST, primer3, and R 3.5.1 to be installed on the system. In addition, it requires R packages optparse, Biostrings, tibble and dplyr.

###Setting up the configuration file

Before using AutoCloner, it is first necessary to setup the configuration file, config.txt. This tells AutoCloner where your genome assembly, in fasta format, and BLAST database for the same assembly are. The file contains three variables, GENOME_N_NAME, which details the name of the genome assembly, GENOME_N_FA_PATH, the path to the fasta file of the assembly and GENOME_N_BLASTDB_PATH, the path to the BLAST database of the genome.

For example:
```
GENOME_1_NAME=IWGSC
GENOME_1_FA_PATH=/path/to/fasta_file
GENOME_1_BLASTDB_PATH=/path/to/blastdb
```
Additional genomes of other varieties of the same species can also be included, in this case the number after “GENOME_” is incrementally increased by one. For example:

```
GENOME_2_NAME=IWGSC
GENOME_2_FA_PATH=/path/to/fasta_file
GENOME_2_BLASTDB_PATH=/path/to/blastdb
GENOME_3_NAME=IWGSC
GENOME_3_FA_PATH=/path/to/fasta_file
GENOME_3_BLASTDB_PATH=/path/to/blastdb
```

###Running AutoCloner

To run AutoCloner, launch main.R with arguments -n for the job name and either -f for the path to the fasta sequence or -a for the path to the multiple sequence alignment if your wish to use your own alignment. AutoCloner returns two sets of primers covering the input sequence such that a redundant set is available in case one of the initial primers does not work.

For example:

./main.R -n myjob1 -f ./myfastafile.fa

Or

./main.R -n myjob1 -a ./mymultiplesequencealignment.fa

###Arguments

AutoCloner also has several other additional arguments:

Minimum product size (-m) – specifies the minimum allowed product size for any of the primers in base pairs. Default value is 800.
Maximum product size (-M) – specifies the maximum allowed product size for any of the primers in base pairs. Default value is 2000.
Start buffer (-s) – specifies the number of bases to include in the multiple sequence alignment upstream of the sequence to be cloned.
End buffer (-e) – specifies the number of bases to include in the multiple sequence alignment downstream of the sequence to be cloned.
Full gene product (-P) – Boolean flag. If present, only one set of primers will be returned whose product includes the entire input sequence.
Only perform primer selection (-O) – Boolean flag. If present, AutoCloner will skip the initial stages of the pipeline and perform primer selection again for the specified job. This is useful if you wish to try out some different settings for an existing job, for example a different maximum product size. NB. If this flag is present, the -f and -a flags are not required.

###Output

AutoCloner outputs several files into the jobs/yourjobname directory:

blast.results – A directory containing the results of the BLAST search in tabular output format.
primers/best.primers – Two sets of Primer3 output files. These two sets are intended to provide redundancy – if primers from the first set do not work then the second set can be used instead.
seq/extended/alignments – Some multiple sequence alignments in Fasta format containing the input sequence, various homologues, varietal sequences and primer sequences.
seq/extended/seqs – The sequences used to produce the multiple sequence alignment in fasta format.
seq/extended/primers.set1 – the sequences of the primers in fasta format.

