write('run.dialign.R', p("jobs/", opt$sequence.name, "/pipeline.checkpoint.txt"))
gene.name = opt$sequence.name
fa.path1 = opt$fasta.path

system(p("scripts/run.dialign.sh jobs/", gene.name, "/"))