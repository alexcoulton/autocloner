#!/usr/bin/Rscript
#load libraries and required scripts

#### INIT ####

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
		make_option(c("--product.full.gene", "-P"), action = "store_true", default = F),
		make_option(c("--min.product.size", "-m"), type = "integer", default = 800),
		make_option(c("--max.product.size", "-M"), type = "integer", default = 1000),
		make_option(c("--start.buffer", "-s"), type = "integer", default = 2000),
		make_option(c("--end.buffer", "-e"), type = "integer", default = 2000),
		make_option(c("--own.alignment.path", "-a"), type = "character", default = ""),
		# make_option(c("--only.msa", "-q"), action = "store_true", default = F),
		# make_option(c("--only.primer.selection", "-O"), action = "store_true", default = F),
		make_option(c("--run.mode", "-r"), type = "character", default = "full") # specify pipeline - "full", "own.align", "only.primer.selection", "only.msa"
		)

	opt = parse_args(OptionParser(option_list = option_list))

	#### DEBUG OPTIONS ####

	if(exists("autocloner.debug")){
		if(autocloner.debug == T){
			opt = list()		
			opt$fasta.path = "debug_seq.fa"
			opt$sequence.name = "debug_seq"			
			opt$product.full.gene = F
			opt$min.product.size = 400
			opt$max.product.size = 2000
			opt$start.buffer = 2000
			opt$end.buffer = 2000
			opt$run.mode = "only.msa"
			# opt$only.primer.selection = F
			opt$own.alignment.path = "" # alignment path
			# opt$only.msa = F
		}
	} 
	
	muscle.debug = T
	blast.debug = T
	

	#### DEFINE FUNCTIONS ####
	print.pipeline.stage = function(x){
		print("")
		print("######################")
		print(paste0("Running ", x))
		print("######################")		
		print("")
	}

	init.new.job = function(){
			system(p("./scripts/create.folder.structure.sh ", opt$sequence.name, " ", base_directory))
			system(p("cp ./scripts/clean.folders.sh ./", "jobs/", opt$sequence.name, "/primers/"))
			system(p("cp ./scripts/run.primer3.sh ./", "jobs/", opt$sequence.name, "/primers/"))
	}

	run.pipeline.own.alignment = function(){
			init.new.job()

			source("./scripts/make.fasta.indexes.rscript.R")		

			file.copy(opt$own.alignment.path, p(job_directory, opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))		
				
			remove.hyphens1 = readLines(p(job_directory, opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))		
			remove.hyphens1 = gsub("-", "", remove.hyphens1)
			writeLines(remove.hyphens1, p(job_directory, opt$sequence.name, "/seq/extended/alignments/all.align.rev.nohyphen.fa"))

			pipeline.stages = c("perform.muscle.R", "rearrange.muscle.rscript.R", "primer.selection.rscript.R",
			  "primer.evaluation.rscript.R", "get.product.sequences.R")

			for(i in pipeline.stages){
				print.pipeline.stage(i)
				source(paste0("./scripts/", i))
			}		
	}

	run.full.pipeline = function(){
		#full run of pipeline
			init.new.job()

			pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
			 "perform.muscle.R", "rearrange.muscle.rscript.R", "primer.selection.rscript.R",
			  "primer.evaluation.rscript.R", "get.product.sequences.R")

			for(i in pipeline.stages){
				print.pipeline.stage(i)
				source(paste0("./scripts/", i))
			}
	}

	run.only.primer.selection = function(){
		pipeline.stages = c("primer.selection.rscript.R",
		  "primer.evaluation.rscript.R", "get.product.sequences.R")

		for(i in pipeline.stages){
			print.pipeline.stage(i)
			source(paste0("./scripts/", i))
		}
	}

	run.only.msa = function(){
		#full run of pipeline
			init.new.job()

			pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
			 "perform.muscle.R", "rearrange.muscle.rscript.R")

			for(i in pipeline.stages){
				print.pipeline.stage(i)
				source(paste0("./scripts/", i))
			}
	}


	#### RUN PIPELINE ####

	if(is.null(opt$sequence.name) == T | ((is.null(opt$fasta.path) == T & opt$run.mode != "own.align") & opt$run.mode != "only.primer.selection")) {
		print("Usage: -n sequence.name -f fasta.path (or -a multiple.alignment.path)")
	} else {
		#load other packages after initial checks 
		suppressMessages(library(Biostrings))
		suppressMessages(library(dplyr))
		suppressMessages(library(tibble))

		if(opt$run.mode == "full") run.full.pipeline()
		if(opt$run.mode == "own.align") run.pipeline.own.alignment()
		if(opt$run.mode == "only.primer.selection") run.only.primer.selection()	
		if(opt$run.mode == "only.msa") run.only.msa()	

	}
}