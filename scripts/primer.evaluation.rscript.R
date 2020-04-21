
gene.name = opt$sequence.name
project.path = base_directory

#the base path for auto primer picker
multiple.alignment.path = p(project.path, "jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa.trimmed.fa")


multi.align = readDNAMultipleAlignment(multiple.alignment.path)
multi.align1 = convert.to.character.data.frame(as.data.frame(as.matrix(multi.align)))



conv.matrix.dnastringset = function(mult.mat, row.names1){
  #converts a matrix to a DNAStringSet object
  #args:
  #	mult.mat - DNAMultipleAlignment object that has been converted to character df w/ convert.to.character.data.frame(as.matrix(dnamultiplealignment.object))
  #	row.names1 - Character vector of names for the sequences
  
  mult.mat3 = apply(mult.mat, 1, function(x) paste(x, collapse = ""))
  
  #convert each object in the list to a DNAStringSet object
  mult.mat4 = lapply(mult.mat3, DNAStringSet)
  
  #make new DNAStringSet object
  mult.mat5 = DNAStringSet()
  
  #add our sequences to the new DNAStringSet object one by one
  for(i in 1:length(mult.mat4)){
    mult.mat5 = c(mult.mat5, mult.mat4[[i]])
  }
  
  names(mult.mat5) = row.names1
  
  mult.mat5
}


#evaluate best primers

for(i in c("set1", "set2")){
	primer.base.path = p(project.path, "jobs/", gene.name, "/primers/best.primers/", i, "/") 
	primer.files = list.files(primer.base.path)	

	if(length(primer.files) == 0){
		print(p("No primers in ", i))
	} else {

	print("primer.files")
	print(primer.files)

	#parse number from start of primer.files filename and order according to this number
	g = strsplit(primer.files, "")
	num.coords = lapply(g, function(x){
	  q = grep("p", x)
	  q[1]
	})

	primer.files.order = unlist(Map(function(x, y){
	  as.numeric(paste(x[1:(y - 1)], collapse = ""))
	}, g, num.coords))

	primer.files = primer.files[match(1:length(primer.files.order), primer.files.order)]
	
	#start main loop
	for(x in primer.files){

		parse.coords = strsplit(x, "\\.")
		parse.coords = parse.coords[[1]][2]
		parse.coords = strsplit(parse.coords, "-")

		# print("parse.coords:")
		# print(parse.coords)



		primer.file1 = readLines(p(primer.base.path, x))
		primer.file1.seq = strsplit(primer.file1[grep("PRIMER_LEFT_0_SEQUENCE", primer.file1)], "=")
		primer.file1.seq = DNAString(primer.file1.seq[[1]][2])		

		primer.file1.multi.forward.pos = strsplit(primer.file1[grep("multiple.alignment.forward.coord", primer.file1)], "=")
		primer.file1.multi.forward.pos = as.numeric(primer.file1.multi.forward.pos[[1]][2])

		# print(primer.file1.seq)

		primer.file1.seq.rev = strsplit(primer.file1[grep("PRIMER_RIGHT_0_SEQUENCE", primer.file1)], "=")
		primer.file1.seq.rev = primer.file1.seq.rev[[1]][2]

		primer.file1.seq.rev = reverseComplement(DNAString(primer.file1.seq.rev))

		primer.file1.multi.reverse.pos = strsplit(primer.file1[grep("multiple.alignment.reverse.coord", primer.file1)], "=")
		primer.file1.multi.reverse.pos = as.numeric(primer.file1.multi.reverse.pos[[1]][2])

		both.seqs = c(DNAStringSet(primer.file1.seq), DNAStringSet(primer.file1.seq.rev))
		names(both.seqs) = c(p("forward.primer.", parse.coords[[1]][1]), p("reverse.primer.", parse.coords[[1]][2]))

		#add our forward primer to the multiple alignment via a matrix
		multi.align1 = add_row(multi.align1)
		multi.align1[nrow(multi.align1), ] = "-"
		f.seq = strsplit(as.character(primer.file1.seq), "")
		f.seq = f.seq[[1]]
		#browser()
		multi.align1[nrow(multi.align1), (primer.file1.multi.forward.pos - (length(primer.file1.seq) - 1)):primer.file1.multi.forward.pos] = f.seq	

		#add our reverse primer to the multiple alignment via a matrix
		multi.align1 = add_row(multi.align1)
		multi.align1[nrow(multi.align1), ] = "-"
		r.seq = strsplit(as.character(primer.file1.seq.rev), "")
		r.seq = r.seq[[1]]
		multi.align1[nrow(multi.align1), primer.file1.multi.reverse.pos:(primer.file1.multi.reverse.pos + (length(primer.file1.seq.rev) - 1))] = r.seq

		multi.align2 = conv.matrix.dnastringset(multi.align1, c(rownames(multi.align), p("forward.primer.", primer.file1.multi.forward.pos), p("reverse.primer.", primer.file1.multi.reverse.pos)))
		
		#return the reverse sequence to correct orientation after adding to 
		#multiple sequence alignment
		both.seqs[2] = reverseComplement(both.seqs[2])

		writeXStringSet(both.seqs, p(project.path, "jobs/", gene.name, "/seq/extended/primers.", i, "/", x, ".seqs.fa"))	
		print(p("Primers written to: ", project.path, gene.name, "/seq/extended/primers.", i, "/", x, ".seqs.fa"))
		writeXStringSet(multi.align2, p(project.path, "jobs/", gene.name, "/seq/extended/alignments/alignment.w.primers", i, ".fa"))	
		

		multi.align = readDNAMultipleAlignment(p(project.path, "jobs/", gene.name, "/seq/extended/alignments/alignment.w.primers", i, ".fa"))
		multi.align1 = convert.to.character.data.frame(as.data.frame(as.matrix(multi.align)))

		# write.csv(multi.align1, p(project.path, "jobs/", gene.name, "/seq/extended/primers.", i, "/", x, ".csv"))

		# writeXStringSet(both.seqs, p(project.path, "jobs/", gene.name, "/seq/extended/primers.", i, "/", x, ".seqs.fa"))	
	}

	#concatenate primer sequences and multiple alignment sequences into one file
	# system(p("cat ", multiple.alignment.path, " ", project.path, gene.name, "/seq/extended/primers.", i, "/* > ",
	# 	project.path, gene.name, "/seq/extended/alignments/alignment.w.primer.", i, ".fa"))

	#use muscle to perform a multiple sequence alignment that includes the primer sequences
	# system(p("muscle -in ", project.path, gene.name, "/seq/extended/alignments/alignment.w.primer.", i,
	#  ".fa -out ", project.path, gene.name, "/seq/extended/alignments/alignment.w.primer.", i, "done.fa"))
	print(p("Alignment written to: ", project.path, gene.name, "/seq/extended/alignments/alignment.w.primers.", i, ".fa"))
	}
}



