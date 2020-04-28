
#BLAST PARSER FOR SCAFFOLDS
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(Biostrings))

args = commandArgs(trailingOnly = T)
gene.name = opt$sequence.name
write('make.fasta.indexes.rscript.R', p("jobs/", opt$sequence.name, "/pipeline.checkpoint.txt"))

#read the configuration file
config.file = readLines(opt$alternate.config)
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

config.settings = config.file[grep(1, config.variables)]

config.settings.temp = config.settings[[1]]
genome.name = strsplit(config.settings.temp[1], "=")
genome.name = genome.name[[1]][2]

fa.index.files = list.files("./fasta.indexes/", pattern = genome.name)
if(length(fa.index.files) == 0){
	print("No indexes found in ./fasta.indexes/; making some...")
	lapply(1:number.genomes, function(x){
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
	  
	  print(p("Making index for ", genome.name))
	  fasta.index1 = fasta.index(fa.path)
	  write.csv(fasta.index1, p("./fasta.indexes/", genome.name, ".fa.idx"), row.names = F)	  
	})
}

