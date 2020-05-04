#generate sequences for testing
library(Biostrings)
seqs = readDNAStringSet("~/bioinf/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")
seqs2 = seqs[sample(length(seqs), 1000)]

for(x in 1:length(seqs2)){
	g = strsplit(names(seqs2[x]), " ")
	g = paste0(g[[1]][1], ".fa")
    writeXStringSet(seqs2[x], paste0("~/livesite/primer/pipeline/test_sequences/", g))
}

