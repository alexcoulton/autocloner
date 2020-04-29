write('trim.alignment.rscript.R', p("jobs/", opt$sequence.name, "/pipeline.checkpoint.txt"))
gene.name = opt$sequence.name
full.gene.product = opt$product.full.gene
min.product.size = opt$min.product.size
max.product.size = opt$max.product.size
start.buffer = opt$start.buffer
end.buffer = opt$end.buffer

#the base path for auto primer picker
project.path = base_directory

config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

system(p(project.path, "jobs/", gene.name, "/primers/clean.folders.sh ", "jobs/", gene.name))

trim.multiple.alignment = function(input.path, main.seq.id, start.cut, end.cut){
  #takes a multiple alignment fasta file and trims it so that all sequences are the same length.
  #output is written to the same directory with .trimmed.fa appeneded to the original file name
  #
  #args:
  # input.path - character string of the path to the multiple alignment fasta file
  # main.seq.id - integer; number of the reference sequence to make cuts from in the fasta file
  # start.cut - integer; number of bases to cut upstream of the start of the reference 
  # end.cut - integer; number of bases to cut downstream of the end of the reference
  
  mult.align = readDNAMultipleAlignment(input.path)
  
  #convert the multiple alignment object to a matrix for easier
  #processing
  mult.mat = as.matrix(mult.align)
  
  #get the start coordinate of the chinese spring gene sequence
  #in the multiple alignment
  #start.main.seq = min(which(mult.mat[main.seq.id, ] == "A"))
  start.main.seq = min(which(mult.mat[main.seq.id, ] == "A" | mult.mat[main.seq.id, ] == "T" | mult.mat[main.seq.id, ] == "G" | mult.mat[main.seq.id, ] == "C"))
  #get the end coordinate
  end.main.seq = max(which(mult.mat[main.seq.id, ] == "A" | mult.mat[main.seq.id, ] == "T" | mult.mat[main.seq.id, ] == "G" | mult.mat[main.seq.id, ] == "C"))
  
  #subset the matrix
  start.cut1 = start.main.seq - start.cut
  end.cut1 = end.main.seq + end.cut
  if(start.cut1 < 1) start.cut1 = 1
  if(end.cut1 > ncol(mult.mat)) end.cut1 = ncol(mult.mat)
 
  mult.mat2 = mult.mat[, (start.cut1):(end.cut1)]
  
  #now begin conversion back to a multiple alignment object
  #first grab the sequences from the matrix and put them into a list
  mult.mat3 = apply(mult.mat2, 1, function(x) paste(x, collapse = ""))
  
  #convert each object in the list to a DNAStringSet object
  mult.mat4 = lapply(mult.mat3, DNAStringSet)
  
  #make new DNAStringSet object
  mult.mat5 = DNAStringSet()
  
  #add our sequences to the new DNAStringSet object one by one
  for(i in 1:length(mult.mat4)){
    mult.mat5 = c(mult.mat5, mult.mat4[[i]])
  }
  
  names(mult.mat5) = rownames(mult.align)
  
  #parse input path to make an output path - the original file with .trimmed.fa appended to the end
  output.path = strsplit(input.path, "\\/")
  output.path = output.path[[1]]
  file.name = output.path[length(output.path)]
  file.name = paste(file.name, ".trimmed.fa", sep = "")
  output.path = p(paste(output.path[1:(length(output.path)-1)], collapse = "/"), "/", file.name)
  
  
  writeXStringSet(mult.mat5, output.path)  
}




#   ____________________________________________________________________________
#   PRIMER3 MULTIPLE ALIGNMENT                                              ####

main.gene.row = 1
template.row = 2
print("test1")
trim.multiple.alignment(p(project.path, "jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa"), main.gene.row, start.buffer, end.buffer)
print("test2")
mult.align_orig = readDNAMultipleAlignment(p(project.path, "jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa"))
mult.align1 = readDNAMultipleAlignment(p(project.path, "jobs/", gene.name, "/seq/extended/alignments/all.align.rev.fa.trimmed.fa"))

homologue.rows = 3:(nrow(mult.align1) - (number.genomes - 1))
variety.rows = (max(homologue.rows) + 1):nrow(mult.align1)
