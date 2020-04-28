write('primer.evaluation.rscript.R', p("jobs/", opt$sequence.name, "/pipeline.checkpoint.txt"))
gene.name = opt$sequence.name
project.path = base_directory

#the base path for auto primer picker
multiple.alignment.path = p(project.path, "jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa")


multi.align = readDNAStringSet(multiple.alignment.path)
# multi.align1 = convert.to.character.data.frame(as.data.frame(as.matrix(multi.align)))

multi.align1 = as.matrix(multi.align)


interleave <- function(v1,v2)                                                                                    
{                                                                                                              
	ord1 <- 2*(1:length(v1))-1
	ord2 <- 2*(1:length(v2))
	c(v1,v2)[order(c(ord1,ord2))]
}   

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

	primer.files.order = sort(sapply(g, function(x) as.numeric(x[1])), index.return = T)$ix
	primer.files = primer.files[primer.files.order]

	multi.align1 = matrix(data = "-", nrow = length(multi.align) + (length(primer.files) * 2), ncol = length(multi.align[[1]]))
	multi.align1[1:length(multi.align), ] = as.matrix(multi.align)
	
	#start main loop	
	primer.index = length(multi.align) + 1
	primer.file1.multi.forward.pos.vector = as.numeric()
	primer.file1.multi.reverse.pos.vector = as.numeric()
	for(x in primer.files){

		parse.coords = strsplit(x, "\\.")
		parse.coords = parse.coords[[1]][2]
		parse.coords = strsplit(parse.coords, "-")

		primer.file1 = readLines(p(primer.base.path, x))
		primer.file1.seq = strsplit(primer.file1[grep("PRIMER_LEFT_0_SEQUENCE", primer.file1)], "=")
		
		primer.file1.seq = DNAString(primer.file1.seq[[1]][2])		

		primer.file1.multi.forward.pos = strsplit(primer.file1[grep("multiple.alignment.forward.coord", primer.file1)], "=")
		primer.file1.multi.forward.pos = as.numeric(primer.file1.multi.forward.pos[[1]][2])
		primer.file1.multi.forward.pos.vector = c(primer.file1.multi.forward.pos.vector, primer.file1.multi.forward.pos)

		primer.file1.seq.rev = strsplit(primer.file1[grep("PRIMER_RIGHT_0_SEQUENCE", primer.file1)], "=")
		primer.file1.seq.rev = primer.file1.seq.rev[[1]][2]

		primer.file1.seq.rev = reverseComplement(DNAString(primer.file1.seq.rev))

		primer.file1.multi.reverse.pos = strsplit(primer.file1[grep("multiple.alignment.reverse.coord", primer.file1)], "=")
		primer.file1.multi.reverse.pos = as.numeric(primer.file1.multi.reverse.pos[[1]][2])
		primer.file1.multi.reverse.pos.vector = c(primer.file1.multi.reverse.pos.vector, primer.file1.multi.reverse.pos)

		both.seqs = c(DNAStringSet(primer.file1.seq), DNAStringSet(primer.file1.seq.rev))
		names(both.seqs) = c(p("forward.primer.", parse.coords[[1]][1]), p("reverse.primer.", parse.coords[[1]][2]))

		#add our forward primer to the multiple alignment via a matrix		
		f.seq = strsplit(as.character(primer.file1.seq), "")
		f.seq = f.seq[[1]]
		
		multi.align1[primer.index, (primer.file1.multi.forward.pos - (length(primer.file1.seq) - 1)):primer.file1.multi.forward.pos] = f.seq	

		#add our reverse primer to the multiple alignment via a matrix		
		r.seq = strsplit(as.character(primer.file1.seq.rev), "")
		r.seq = r.seq[[1]]
		multi.align1[primer.index + 1, primer.file1.multi.reverse.pos:(primer.file1.multi.reverse.pos + (length(primer.file1.seq.rev) - 1))] = r.seq

		
		
		#return the reverse sequence to correct orientation after adding to 
		#multiple sequence alignment
		both.seqs[2] = reverseComplement(both.seqs[2])

		writeXStringSet(both.seqs, p(project.path, "jobs/", gene.name, "/seq/extended/primers.", i, "/", x, ".seqs.fa"))	
		print(p("Primers written to: ", project.path, gene.name, "/seq/extended/primers.", i, "/", x, ".seqs.fa"))
		primer.index = primer.index + 2
	}
	
	primer.file1.multi.forward.pos.vector = paste0("forward.", primer.file1.multi.forward.pos.vector)
	primer.file1.multi.reverse.pos.vector = paste0("reverse.", primer.file1.multi.reverse.pos.vector)
	primer.row.names = interleave(primer.file1.multi.forward.pos.vector, primer.file1.multi.reverse.pos.vector)
	multi.align2 = conv.matrix.dnastringset(multi.align1, c(names(multi.align), primer.row.names))
	writeXStringSet(multi.align2, p(project.path, "jobs/", gene.name, "/seq/extended/alignments/alignment.w.primers", i, ".fa"))	
	
	#concatenate primer sequences and multiple alignment sequences into one file
	# system(p("cat ", multiple.alignment.path, " ", project.path, gene.name, "/seq/extended/primers.", i, "/* > ",
	# 	project.path, gene.name, "/seq/extended/alignments/alignment.w.primer.", i, ".fa"))

	#use muscle to perform a multiple sequence alignment that includes the primer sequences
	# system(p("muscle -in ", project.path, gene.name, "/seq/extended/alignments/alignment.w.primer.", i,
	#  ".fa -out ", project.path, gene.name, "/seq/extended/alignments/alignment.w.primer.", i, "done.fa"))
	print(p("Alignment written to: ", project.path, gene.name, "/seq/extended/alignments/alignment.w.primers.", i, ".fa"))
	}
}





all.forward.primers = read.csv(p(project.path, "jobs/", gene.name, "/primers/penalties/forward.all.pen.csv"))
all.forward.primers.fasta = interleave(paste0(">", all.forward.primers$p.name), all.forward.primers$PRIMER_LEFT_0_SEQUENCE)
writeLines(all.forward.primers.fasta, p(project.path, "jobs/", gene.name, "/primers/penalties/all.forward.seq.fa"))

all.reverse.primers = read.csv(p(project.path, "jobs/", gene.name, "/primers/penalties/reverse.all.pen.csv"))
all.reverse.primers.fasta = interleave(paste0(">", all.reverse.primers$p.name), all.reverse.primers$PRIMER_RIGHT_0_SEQUENCE)
writeLines(all.reverse.primers.fasta, p(project.path, "jobs/", gene.name, "/primers/penalties/all.reverse.seq.fa"))