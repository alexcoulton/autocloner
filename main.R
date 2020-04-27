#!/usr/bin/Rscript
#load libraries and required scripts

#### INIT ####

options(show.error.locations = T)
options(stringsAsFactors = F)
options(error = function() { traceback(2); if(!interactive()) quit('no', status = 1, runLast = FALSE) })

suppressMessages(library("optparse"))
source("./scripts/functions.R")
print("current working directory:")
print(getwd())

base_directory = paste(getwd(), "/", sep = "")
job_directory = paste(base_directory, "jobs", "/", sep = "")

setwd(base_directory)

args = commandArgs(trailingOnly = T)

option_list = list(
		make_option(c("--sequence.name", "-n"), type = "character"),
		make_option(c("--fasta.path", "-f"), type = "character"),		
		make_option(c("--product.full.gene", "-P"), action = "store_true", default = F),
		make_option(c("--min.product.size", "-m"), type = "integer", default = 800),
		make_option(c("--max.product.size", "-M"), type = "integer", default = 1000),
		make_option(c("--start.buffer", "-s"), type = "integer", default = 2000),
		make_option(c("--end.buffer", "-e"), type = "integer", default = 2000),
		make_option(c("--own.alignment.path", "-a"), type = "character", default = ""),		
		make_option(c("--run.mode", "-r"), type = "character", default = "run.full.pipeline"), # specify pipeline - "full", "own.align", "only.primer.selection", "only.msa"
		make_option(c("--perform.masking"), action = "store_true", default = F), 
		make_option(c("--mask.bin.size"), type = "integer", default = 10), 
		make_option(c("--mask.threshold"), type = "integer", default = 40), 
		make_option(c("--allow.hyphens.in.mask"), action = "store_true", default = F),
		make_option(c("--number.homologues"), type = "integer", default = 4),
		make_option(c("--all.homologues"), action = "store_true", default = F),
		make_option(c("--allow.hyphens.for.snp.detection"), action = "store_true", default = T),
		make_option(c("--mask.inter.hsp.distances"), action = "store_true", default = F),
		make_option(c("--alignment.method"), type = "character", default = "muscle"),
		make_option(c("--cds.max.intron.size"), type = "integer", default = 3500),
		make_option(c("--alternate.config"), type = "character", default = './config.txt')
		)

opt = parse_args(OptionParser(option_list = option_list))	

config.file = readLines(opt$alternate.config)
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

if(number.genomes < 1) {
	print("./config.txt is empty")
} else {
	
	#### DEBUG OPTIONS ####

	if(exists("autocloner.debug")){
		if(autocloner.debug == T){
			opt = list()		
			opt$fasta.path = "debug_seq.fa"
			# opt$sequence.name = "debug_seq"			
			opt$sequence.name = "test7"		
			opt$product.full.gene = F
			opt$min.product.size = 400
			opt$max.product.size = 2000
			opt$start.buffer = 2000
			opt$end.buffer = 2000
			
			#ALL RUN MODES
			opt$run.mode = "run.full.pipeline"
			# opt$run.mode = "run.pipeline.own.alignment"
			# opt$run.mode = "run.only.primer.selection"
			# opt$run.mode = "run.only.msa"
			# opt$run.mode = "run.only.snp.selection"
			# opt$run.mode = "run.adv.primer.select.job"					

			opt$own.alignment.path = "" # alignment path			
			opt$perform.masking = F
			opt$mask.bin.size = 10
			opt$mask.threshold = 40
			opt$allow.hyphens.in.mask = "F"
			opt$number.homologues = 4			
			opt$all.homologues = T
			opt$allow.hyphens.for.snp.detection = T
			opt$alignment.method = "muscle"
			opt$mask.inter.hsp.distances = T
			opt$cds.max.intron.size = 3500		
			# opt$alternate.config = './configs/b.rapa.config.txt'	
			opt$alternate.config = './config.txt'	
		}
	} 
	
	muscle.debug = F
	blast.debug = F
	

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
			
			if(opt$alignment.method == "muscle"){
				pipeline.stages = c("perform.muscle.R", "rearrange.muscle.rscript.R", "trim.alignment.rscript.R", "primer.selection.rscript.R",
			  		"primer.evaluation.rscript.R", "get.product.sequences.R")
			} else {
				pipeline.stages = c("run.dialign.R", "primer.selection.rscript.R",
			  		"primer.evaluation.rscript.R", "get.product.sequences.R")
			}


			for(i in pipeline.stages){
				print.pipeline.stage(i)
				source(paste0("./scripts/", i))
			}		
	}

	run.full.pipeline = function(){
		#full run of pipeline
			init.new.job()

			if(opt$alignment.method == "muscle"){
				pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
					"perform.muscle.R", "rearrange.muscle.rscript.R", "trim.alignment.rscript.R", "primer.selection.rscript.R",
					"primer.evaluation.rscript.R", "get.product.sequences.R")
			} else {
				pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
					"run.dialign.R", "primer.selection.rscript.R",
					"primer.evaluation.rscript.R", "get.product.sequences.R")
			}

			# pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
			# "run.dialign.R")


			for(i in pipeline.stages){
				print.pipeline.stage(i)
				source(paste0("./scripts/", i))
			}
	}

	run.single.pipeline.stage = function(script.name){
		source(paste0('./scripts/', script.name))
	}

	run.adv.primer.select.job = function(){		
		ONLY.PRIMER.SELECTION <<- T
			pipeline.stages = c("primer.selection.rscript.R",
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

	run.only.snp.selection = function(){
		ONLY.SNP.SELECTION <<- T
		pipeline.stages = c("primer.selection.rscript.R")
		for(i in pipeline.stages){
			print.pipeline.stage(i)
			source(paste0("./scripts/", i))
		}
	}

	run.only.msa = function(){
		#full run of pipeline
			init.new.job()			

			if(opt$alignment.method == "muscle"){
				pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
			 		"perform.muscle.R", "rearrange.muscle.rscript.R", "trim.alignment.rscript.R")
			} else {
				pipeline.stages = c("make.fasta.indexes.rscript.R", "perform.blast.rscript.R", "blast.scaffold.parser.rscript.R",
					"run.dialign.R")
			}

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

		if(opt$run.mode == "run.full.pipeline") run.full.pipeline()
		if(opt$run.mode == "run.pipeline.own.alignment") run.pipeline.own.alignment()
		if(opt$run.mode == "run.only.primer.selection") run.only.primer.selection()	
		if(opt$run.mode == "run.only.msa") run.only.msa()	
		if(opt$run.mode == "run.only.snp.selection") run.only.snp.selection()	
		if(opt$run.mode == "run.adv.primer.select.job") run.adv.primer.select.job()	

	}
}
