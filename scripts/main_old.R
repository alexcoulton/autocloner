#!/home/ac14037/bin/Rscript

base_directory = paste(getwd(), "/", sep = "")
job_directory = paste(base_directory, "jobs", "/", sep = "")

source("./functions.R")
setwd(base_directory)

suppressMessages(library("optparse"))

args = commandArgs(trailingOnly = T)

config.file = readLines("./config.txt")
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

if(number.genomes < 1){
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

opt = parse_args(OptionParser(option_list = option_list))

# print(opt)

if(is.null(opt$sequence.name) == T | ((is.null(opt$fasta.path) == T & opt$own.alignment == "") & opt$only.primer.selection == F)) {
	print("Usage: -n sequence.name -f fasta.path (or -a multiple.alignment.path)")
} else {
	if(opt$own.alignment != ""){
		system(p("./create.folder.structure.sh ", opt$sequence.name, " ", base_directory))
		system(p("cp ./clean.folders.sh ./", "jobs/", opt$sequence.name, "/primers/"))
		system(p("cp ./run.primer3.sh ./", "jobs/", opt$sequence.name, "/primers/"))

		system(p("./make.fasta.indexes.rscript.R ", base_directory))

		file.copy(opt$own.alignment, p(job_directory, opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))		

		print("######################")
		print("Running primer.selection.script.R")
		print("######################")
		print("")
		system(p("./primer.selection.rscript.R ", opt$sequence.name, " ", opt$product.full.gene, " ",
		opt$min.product.size, " ", opt$max.product.size, " ", base_directory, " ", opt$start.buffer, " ", opt$end.buffer))



		print("######################")
		print("Running primer.evaluation.rscript.R")
		print("######################")
		print("")
		system(p("./primer.evaluation.rscript.R ", opt$sequence.name, " ", base_directory))	
		print("######################")
		print("Running get.product.sequences.R")
		print("######################")
		print("")
		system(p("./get.product.sequences.R ", opt$sequence.name, " ", opt$base_directory))
	} else {

		if(opt$only.primer.selection == F){
				system(p("./create.folder.structure.sh ", opt$sequence.name, " ", base_directory))
				system(p("cp ./clean.folders.sh ./", "jobs/", opt$sequence.name, "/primers/"))
				system(p("cp ./run.primer3.sh ./", "jobs/", opt$sequence.name, "/primers/"))

				system(p("./make.fasta.indexes.rscript.R ", base_directory))

				print("")
				print("######################")
				print("Running perform.blast.rscript.R")
				print("######################")		
				print("")
				system(p("./perform.blast.rscript.R ", opt$sequence.name, " ", opt$fasta.path, " ", base_directory))	

				print("")
				print("######################")
				print("Running blast.scaffold.parser.rscript.R")
				print("######################")
				print("")
				system(p("./blast.scaffold.parser.rscript.R ", opt$sequence.name, " ", opt$fasta.path, " ", base_directory))	

				print("")
				print("######################")
				print("Running muscle to generate alignment")
				print("######################")
				print("")
				system(p("muscle -in ", "jobs/", opt$sequence.name, "/seq/extended/seqs/all.fa -out ", "jobs/", opt$sequence.name, "/seq/extended/alignments/all.align.rev.fa"))
				system(p("./rearrange.muscle.rscript.R ", opt$sequence.name, " ", base_directory))

				print("######################")
				print("Running primer.selection.script.R")
				print("######################")
				print("")
				system(p("./primer.selection.rscript.R ", opt$sequence.name, " ", opt$product.full.gene, " ",
				opt$min.product.size, " ", opt$max.product.size, " ", base_directory, " ", opt$start.buffer, " ", opt$end.buffer))

				print("######################")
				print("Running primer.evaluation.rscript.R")
				print("######################")
				print("")
				system(p("./primer.evaluation.rscript.R ", opt$sequence.name, " ", base_directory))	
				print("######################")
				print("Running get.product.sequences.R")
				print("######################")
				print("")
				system(p("./get.product.sequences.R ", opt$sequence.name, " ", opt$base_directory))
			} else {
				print("######################")
				print("Running primer.selection.script.R")
				print("######################")
				print("")
				system(p("./primer.selection.rscript.R ", opt$sequence.name, " ", opt$product.full.gene, " ",
				opt$min.product.size, " ", opt$max.product.size, " ", base_directory, " ", opt$start.buffer, " ", opt$end.buffer))

				print("######################")
				print("Running primer.evaluation.rscript.R")
				print("######################")
				print("")
				system(p("./primer.evaluation.rscript.R ", opt$sequence.name, " ", base_directory))	
				print("######################")
				print("Running get.product.sequences.R")
				print("######################")
				print("")
				system(p("./get.product.sequences.R ", opt$sequence.name, " ", base_directory))
			}



		}

	
	
}

}