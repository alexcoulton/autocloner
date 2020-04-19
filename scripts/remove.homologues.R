paste0("jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa")




system(p("muscle -in ", "jobs/", opt$sequence.name, "/seq/extended/seqs/all.fa -out ", "jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))	