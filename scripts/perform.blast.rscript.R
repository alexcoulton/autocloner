
#provide path of file to be BLASTed against genomes as an argument

#perform BLAST search of sequence against genomes for primer picking



gene.name = opt$sequence.name
query.fa.path = opt$fasta.path

file.copy(query.fa.path, p("jobs/", gene.name, "/seq/extended/seqs/input_seq.fa"))

query.fa.name = strsplit(query.fa.path, "/")
query.fa.name = query.fa.name[[1]]
query.fa.name = query.fa.name[length(query.fa.name)]

#read the configuration file
config.file = readLines(opt$alternate.config)
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))




#make arbitrary assignment so lapply doesn't print
g = lapply(1:number.genomes, function(x){
  #begin by parsing the config file for the genome name, fasta file path and blastdb path
  config.settings = config.file[grep(x, config.variables)]
  
  config.settings.temp = config.settings[[1]]
  genome.name = strsplit(config.settings.temp[1], "=")
  genome.name = genome.name[[1]][2]
  
  config.settings.temp = config.settings[[2]]
  fa.path = strsplit(config.settings.temp[1], "=")
  fa.path = fa.path[[1]][2]
  
  config.settings.temp = config.settings[[3]]
  blastdb.path = strsplit(config.settings.temp[1], "=")
  blastdb.path = blastdb.path[[1]][2]
  
  print(p("Running ", query.fa.name, " vs ", genome.name, " BLAST"))

if(blast.debug == T){
  file.copy("debug_files/debug.blast", paste0("./jobs/", gene.name, "/blast.results/", "1", ".", query.fa.name, ".vs.", genome.name, ".blast"))
} else {
  system(p("blastn -db ", blastdb.path, " -query ", query.fa.path,
           " -outfmt 6 -out ./jobs/", gene.name, "/blast.results/", x, ".", query.fa.name, ".vs.", genome.name, ".blast",
           " -num_threads 6 -culling_limit 10"))
}

})


