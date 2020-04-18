if(opt$run.mode == "own.align"){ #if using own alignment
	system(p("muscle -in ", "jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.nohyphen.fa -out ", "jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))
} else {

	if(muscle.debug == T){
	file.copy("debug_files/all.align.rev.fa", paste0("jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))
	} else {
	system(p("muscle -in ", "jobs/", opt$sequence.name, "/seq/extended/seqs/all.fa -out ", "jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))	
	}

	}

