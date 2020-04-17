
#perform BLAST search of sequence against genomes for primer picking
#provide path of file to be BLASTed against genomes as an argument

gene.name = opt$sequence.name

if(opt$own.alignment != ""){
	orig.seqs = readDNAStringSet(p("jobs/", gene.name, "/seq/extended/alignments/all.align.rev.nohyphen.fa"))
	muscle.orig = readDNAStringSet(p("jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa"))

	muscle.orig = muscle.orig[match(names(orig.seqs), names(muscle.orig))]
	writeXStringSet(muscle.orig, p("jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa"))
} else {
	orig.seqs = readDNAStringSet(p("jobs/", gene.name, "/seq/extended/seqs/all.fa"))
	muscle.orig = readDNAStringSet(p("jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa"))

	muscle.orig = muscle.orig[match(names(orig.seqs), names(muscle.orig))]
	writeXStringSet(muscle.orig, p("jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa"))
}




