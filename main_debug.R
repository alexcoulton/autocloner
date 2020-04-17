#!/usr/bin/Rscript

#load libraries and required scripts


options(show.error.locations = T)
options(error = function() { traceback(2); if(!interactive()) quit('no', status = 1, runLast = FALSE) })

suppressMessages(library("optparse"))
source("./scripts/functions.R")
print("current working directory:")
print(getwd())

base_directory = paste(getwd(), "/", sep = "")
job_directory = paste(base_directory, "jobs", "/", sep = "")

setwd(base_directory)

args = commandArgs(trailingOnly = T)

config.file = readLines("./config.txt")
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

if(number.genomes < 1) {
	print("./config.txt is empty")
} else {

option_list = list(
	make_option(c("--sequence.name", "-n"), type = "character"),
	make_option(c("--fasta.path", "-f"), type = "character"),
	make_option(c("--only.primer.selection", "-O"), action = "store_true", default = F),
	make_option(c("--product.full.gene", "-P"), action = "store_true", default = F),
	make_option(c("--min.product.size", "-m"), type = "integer", default = 800),
	make_option(c("--max.product.size", "-M"), type = "integer", default = 1000),
	make_option(c("--start.buffer", "-s"), type = "integer", default = 2000),
	make_option(c("--end.buffer", "-e"), type = "integer", default = 2000),
	make_option(c("--own.alignment", "-a"), type = "character", default = "")
	)

#opt = parse_args(OptionParser(option_list = option_list))

opt = list()

opt$fasta.path = "debug_seq.fa"
opt$sequence.name = "debug_seq"
opt$only.primer.selection = F
opt$product.full.gene = F
opt$min.product.size = 400
opt$max.product.size = 2000
opt$start.buffer = 2000
opt$end.buffer = 2000
opt$own.alignment = ""

primer.selection.script = "./scripts/primer.selection.rscript.R"

# print(opt)

if(is.null(opt$sequence.name) == T | ((is.null(opt$fasta.path) == T & opt$own.alignment == "") & opt$only.primer.selection == F)) {
	print("Usage: -n sequence.name -f fasta.path (or -a multiple.alignment.path)")
} else {
	#load other packages after initial checks 
	suppressMessages(library(Biostrings))
	suppressMessages(library(dplyr))
	suppressMessages(library(tibble))

	if(opt$own.alignment != ""){
		system(p("./scripts/create.folder.structure.sh ", opt$sequence.name, " ", base_directory))
		system(p("cp ./scripts/clean.folders.sh ./", "jobs/", opt$sequence.name, "/primers/"))
		system(p("cp ./scripts/run.primer3.sh ./", "jobs/", opt$sequence.name, "/primers/"))

		source("./scripts/make.fasta.indexes.rscript.R")		

		file.copy(opt$own.alignment, p(job_directory, opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))		

		print("######################")
		print("Running primer.selection.rscript.R")
		print("######################")
		print("")
		source(primer.selection.script)



		print("######################")
		print("Running primer.evaluation.rscript.R")
		print("######################")
		print("")
		source("./scripts/primer.evaluation.rscript.R")

		print("######################")
		print("Running get.product.sequences.R")
		print("######################")
		print("")
		source("./scripts/get.product.sequences.R")
	} else {

		if(opt$only.primer.selection == F){
				system(p("./scripts/create.folder.structure.sh ", opt$sequence.name, " ", base_directory))
				system(p("cp ./scripts/clean.folders.sh ./", "jobs/", opt$sequence.name, "/primers/"))
				system(p("cp ./scripts/run.primer3.sh ./", "jobs/", opt$sequence.name, "/primers/"))

				source("./scripts/make.fasta.indexes.rscript.R")				

				print("")
				print("######################")
				print("Running perform.blast.rscript.R")
				print("######################")		
				print("")
				source("./scripts/perform.blast.rscript.R")

				print("")
				print("######################")
				print("Running blast.scaffold.parser.rscript.R")
				print("######################")
				print("")
				source("./scripts/blast.scaffold.parser.rscript.R")

				print("")
				print("######################")
				print("Running muscle to generate alignment")
				print("######################")
				print("")
				system(p("muscle -in ", "jobs/", opt$sequence.name, "/seq/extended/seqs/all.fa -out ", "jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))
				source("./scripts/rearrange.muscle.rscript.R")

				browser()

				print("######################")
				print("Running primer.selection.rscript.R")
				print("######################")
				print("")
				source(primer.selection.script)

				print("######################")
				print("Running primer.evaluation.rscript.R")
				print("######################")
				print("")
				source("./scripts/primer.evaluation.rscript.R")
				print("######################")
				print("Running get.product.sequences.R")
				print("######################")
				print("")
				source("./scripts/get.product.sequences.R")
			} else {
				print("######################")
				print("Running primer.selection.rscript.R")
				print("######################")
				print("")
				source(primer.selection.script)

				print("######################")
				print("Running primer.evaluation.rscript.R")
				print("######################")
				print("")
				source("./scripts/primer.evaluation.rscript.R")
				print("######################")
				print("Running get.product.sequences.R")
				print("######################")
				print("")
				source("./scripts/get.product.sequences.R")
			}



		}

	
	
}

}
