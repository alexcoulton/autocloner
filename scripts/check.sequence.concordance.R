alignment = readDNAStringSet(p("jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))

seq1 = strsplit(as.character(alignment[[1]]), "")[[1]]
seq2 = strsplit(as.character(alignment[[2]]), "")[[1]]

hyphencoord = which(seq1 == "-")

seq1 = seq1[-hyphencoord]
seq2 = seq2[-hyphencoord]

input_template_sim = (length(which(seq1 == seq2)) / length(seq1)) * 100

# blast.files = list.files(p("jobs/", opt$sequence.name, "/blast.results"), full.names = T)
# bf2 = blast.files[grep("IWGSC", blast.files)]
# bf2 = bf2[-grep("w_flanking", blast.files[grep("IWGSC", blast.files)])]

# bf3 = read.blast(bf2)
# top_blast_hit_sim = bf3[1, ]$percentage.identical

if(input_template_sim < 95){
    write('Input sequence less then 95% similar to extracted sequence', p("jobs/", opt$sequence.name, "/low.sim.input.extract.txt"))
    rm.error.txt()
    stop('Input sequence less then 95% similar to extracted sequence')    
}