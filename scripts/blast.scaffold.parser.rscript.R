#BLAST PARSER FOR SCAFFOLDS
write('blast.scaffold.parser.rscript.R', p("jobs/", opt$sequence.name, "/pipeline.checkpoint.txt"))
gene.name = opt$sequence.name
fa.path1 = opt$fasta.path
input_sequence = readDNAStringSet(p("jobs/", gene.name, "/seq/extended/seqs/input_seq.fa"))

#   ____________________________________________________________________________
#   DEFINE FUNCTIONS                                                        ####

sort.blastdf = function(blastdf){
  #first clusters BLAST hits by chromosome, then sorts by the start location of each hit within clusters
  #blastdf - a dataframe containing BLAST output in tabular format
  sorted = newdf(colnames(blastdf), no.rows = T)
  blastdf$orientation = "F"
  rev.coords = which(blastdf$sstart > blastdf$send)
  blastdf$orientation[rev.coords] = "R"

  rev.starts = blastdf$sstart[rev.coords]
  rev.ends = blastdf$send[rev.coords]

  blastdf$sstart[rev.coords] = rev.ends
  blastdf$send[rev.coords] = rev.starts
  
  #do some sorting
  for(i in unique(blastdf[, 2])){
    temp = blastdf[blastdf[, 2] == i, ]
    temp = temp[sort(as.numeric(temp[, 9]), index.return = T)$ix, ]
    sorted = rbind(sorted, temp)
  }
  return(sorted)
}

parse.scaffold.blast = function(blastdf1, dist.threshold){
  #parses a BLAST dataframe of a short query sequence against a genome assembly
  #composed of scaffolds or chromosomes. If the assemblie is chromosomal, the parser will split 
  #the chromosome up into groups of hits where hits are more than dist.threshold bp apart.
  #returns a dataframe containing the best groups of hits (average bitscore higher than 200, individual hits no more than dist.threshold bp apart)
  #args:
  # blastdf1 - a BLAST dataframe imported using read.blast()
  # dist.threshold - Integer; the maximum number of bases between two hits for them to be considered part of the same group
  
  blastdf_orig = blastdf1
  blastdf1 = sort.blastdf(blastdf1) 
  
  unique.groups = convert.to.character.data.frame(unique(blastdf1[, 1:2]))
  
  potential.homeologues = newdf(c("query", "scaffold", "start", "end", "length", "rev.comp", "avg.bitscore"), no.rows = T)
  
  count1 = make.counter()
  
  for(i in 1:nrow(unique.groups)){
    temp.df = filter(blastdf1, qseqid == unique.groups[i, 1], sseqid == unique.groups[i, 2])

    split.numeric.vectorv2 = function(sstart, send, threshold){        
        sstart = sstart[-1]
        send = send[-length(send)]
        
        g = data.frame(send, sstart)
        g.diffs = abs(g$sstart - g$send)
        
        cons1 = function(x){
            #get consecutive integer ranges / integer runs
            diffs = c(1, diff(x))
            start_indexes = c(1, which(diffs > 1))
            end_indexes = c(start_indexes - 1, length(x))            
            g = data.frame(x[start_indexes], x[end_indexes])
            colnames(g) = c("start", "end")
            g
        }
        
        group.coords = which(g.diffs < threshold)
        if(length(group.coords) > 0){
            groups1 = cons1(group.coords)
            groups2 = list()
            for(i in 1:nrow(groups1)){
                groups2 = c(groups2, list(c(groups1[i, 1]:groups1[i, 2], groups1[i, 2] + 1)))
            }

            all.rows = 1:(length(sstart) + 1)
            all.rows = all.rows[-which(all.rows %in% unlist(groups2))]
            all.groups = lapply(all.rows, function(x) x)
            all.groups = c(all.groups, groups2)
        } else {
            all.rows = 1:(length(sstart) + 1)            
            all.groups = lapply(all.rows, function(x) x)            
        }
        
        all.groups        
    }
    
    split.temp.df = function(temp.df, orientation1){
        temp.df.corrected = temp.df[which(temp.df$orientation == orientation1), ]    
        #determine whether there is more than one locus involved in this group of hits
        correct.groups = split.numeric.vectorv2(temp.df.corrected$sstart, temp.df.corrected$send, dist.threshold)        
        for(x in 1:length(correct.groups)){
            temp.df.corrected$sseqid[correct.groups[[x]]] = paste0(temp.df.corrected$sseqid[correct.groups[[x]]], ".!!$", orientation1, x)
        }
        
        temp.df.corrected
    }

    rev.orientation.coords = which(temp.df$orientation == "R")
    forward.orientation.coords = which(temp.df$orientation == "F")
    if(length(rev.orientation.coords) > 0){        
        temp.df[rev.orientation.coords, ] = split.temp.df(temp.df, "R")
    }

    if(length(forward.orientation.coords) > 0){
        temp.df[forward.orientation.coords, ] = split.temp.df(temp.df, "F")
    }
    #if two groups of hits are present in the same scaffold / chromosome,
    #these hits will no longer be ordered by bitscore due to sort.blastdf()
    #at the start of the script. Here we calculate mean bitscore for these newly identified
    #groups of hits, and sort groups of hits by this in descending order.
    
    temp.df.unique.scaffolds = unique(temp.df$sseqid)
    mean.bitscores = unlist(lapply(temp.df.unique.scaffolds, function(x){
      temp.df.filtered = filter(temp.df, sseqid == x)
      mean(as.numeric(temp.df.filtered$bitscore))
    }))
    
    transformation.coords = sort(mean.bitscores, decreasing = T, index.return = T)$ix 
        
    transformation.coords2 = unlist(lapply(temp.df.unique.scaffolds[transformation.coords], function(x){
      which(x == temp.df$sseqid)
    }))
    
    temp.df = temp.df[transformation.coords2, ]        

	check.group.orientation = split(temp.df, factor(temp.df$sseqid, levels = unique(temp.df$sseqid)))

	check.group.orientation = lapply(check.group.orientation, function(x){      
	  #if HSPs are near to each other but in different orientations, separate them into different groups
	  same.orientation = (length(unique(x$orientation)) == 1)
	  if(same.orientation == F){
		group1 = x[which(x$sstart < x$send), ]
		group1$sseqid = paste0(group1$sseqid, "_1")
		group2 = x[which(!x$sstart < x$send), ]
		group2$sseqid = paste0(group2$sseqid, "_2")
		x = bind_rows(group1, group2)
	  }
	  x
	})

	temp.df = bind_rows(check.group.orientation)
    blastdf1[which(blastdf1$qseqid == unique.groups[i, 1] & blastdf1$sseqid == unique.groups[i, 2]), ] = temp.df 
    
    for(i2 in 1:length(unique(temp.df$sseqid))){
      #CONCATENATE GROUPS OF HITS TOGETHER INTO potential.homeologues DATAFRAME
      temp.df2 = filter(temp.df, sseqid == unique(temp.df$sseqid)[i2])
      temp.df2$qseqid = as.character(temp.df2$qseqid)
      temp.df2$sseqid = as.character(temp.df2$sseqid)      
      
      group.orientation = temp.df2$orientation[1]         
      
      #populate potential.homeologues dataframe where average bitscore is higher than 200  
      if(mean(temp.df2$bitscore) > 200){        
          #if this is in normal orientation, do x...  
          potential.homeologues = add_row(potential.homeologues)
          potential.homeologues$query[nrow(potential.homeologues)] = temp.df2[1, 1]
          potential.homeologues$scaffold[nrow(potential.homeologues)] = temp.df2[1, 2]
          potential.homeologues$start[nrow(potential.homeologues)] = min(temp.df2$sstart)
          potential.homeologues$end[nrow(potential.homeologues)] = max(temp.df2$send)
          potential.homeologues$avg.bitscore[nrow(potential.homeologues)] = mean(temp.df2$bitscore)
          potential.homeologues$max.bitscore[nrow(potential.homeologues)] = max(temp.df2$bitscore)
          potential.homeologues$avg.percent.identical[nrow(potential.homeologues)] = mean(temp.df2$percentage.identical)
          if(group.orientation == "F"){
              potential.homeologues$rev.comp[nrow(potential.homeologues)] = F
          } else {
              potential.homeologues$rev.comp[nrow(potential.homeologues)] = T
          }
          
          potential.homeologues$query.start[nrow(potential.homeologues)] = min(temp.df2$qstart)
          potential.homeologues$query.end[nrow(potential.homeologues)] = max(temp.df2$qend)
          potential.homeologues$num_hsp[nrow(potential.homeologues)] = nrow(temp.df2)
      }
      
      
      
      
    }
    
  }

    

  potential.homeologues$length = as.numeric(potential.homeologues$end) - as.numeric(potential.homeologues$start)
  potential.homeologues$groupid = potential.homeologues$scaffold
  potential.homeologues$scaffold = multi.str.split(potential.homeologues$scaffold, "\\.\\!\\!\\$", 1)    
  potential.homeologues$homo_length = potential.homeologues$query.end - potential.homeologues$query.start

  # try and identify the matching genomic sequence to the input sequence - avg.bitscore sometimes fails here
  # e.g. if there is a small exon seperated by an intron from the main sequence, it will bring the avg.bitscore down 
  identi.coord = which.max((potential.homeologues$homo_length / length(input_sequence[[1]])) * (potential.homeologues$avg.percent.identical / 100)) 
  g = 1:nrow(potential.homeologues)
  g = g[-identi.coord]  
  potential.homeologues = potential.homeologues[c(identi.coord, g), ] 
  
  if(length(input_sequence[[1]]) > 1500){
    #this will remove all blast hits for small sequences. need an if statement
    coord_to_rm = which(potential.homeologues$length < 500)
    if(length(coord_to_rm) != nrow(potential.homeologues)){
      if(length(coord_to_rm) > 0) potential.homeologues = potential.homeologues[-coord_to_rm, ]           
    }
    
  }

  existing.homeologue.files = grep('potential_homeologues', list.files(p('jobs/', gene.name, '/blast.results/')))
  if(length(existing.homeologue.files) == 0){
    write.csv(potential.homeologues, p('jobs/', gene.name, '/blast.results/potential_homeologues1.csv'), row.names = F)
  } else {
    write.csv(potential.homeologues, p('jobs/', gene.name, '/blast.results/potential_homeologues', (length(existing.homeologue.files) + 1), '.csv'), row.names = F)
  }


  
  list(potential.homeologues, blastdf1)
}



extract.sequence = function(genome1, blast.df.parsed, row.coords, start.buffer, end.buffer){
  #extracts a related sequence from a genome assembly
  #args:
  # genome1 - DNAStringSet object containing the genome assembly of interest
  # blast.df.parsed - dataframe produced by parse.scaffold.blast()
  # row.coords - Numeric vector; the row coordinates of the sequences to used in blast.df.parsed
  # start.buffer - Integer; how much extra sequence before the start indicated in blast.df.parsed to extract
  # end.buffer - Integer; how much extra sequence after the end indicated in blast.df.parsed to extract
  
  #parse the original scaffold name from blastdf1.parsed (remove the appended position in kb)
  original.scaf.names = multi.str.split(blast.df.parsed$scaffold, ".$!", 1) 

  genome2 = genome1[match(original.scaf.names[row.coords], names(genome1))]  
  
  for(i in 1:length(row.coords)){
    #check if start.buffer reaches before the start of the scaffold, and likewise if end.buffer extends bigger than the total length
    extract.start = (as.numeric(blast.df.parsed$start[i]) - start.buffer)
    if(extract.start < 1) extract.start = 1
    extract.end = (as.numeric(blast.df.parsed$end[i]) + end.buffer)
    if(extract.end > length(genome2[[i]])) extract.end = length(genome2[[i]])    
    print("extract.start")
    print(extract.start)
    genome2[[i]] = genome2[[i]][extract.start:extract.end]
    if(blast.df.parsed$rev.comp[i] == T) genome2[[i]] = reverseComplement(genome2[[i]])
  }
  
  #append the start coordinate (in kb) of the blast hit to the name of the sequence for identification later
  names(genome2) = paste(blast.df.parsed$scaffold[row.coords], ".$!", round((as.numeric(blast.df.parsed$start[row.coords]) / 1000)), sep = "")
  print("genome2")  
  print(genome2) 
    genome2  
}

#   ____________________________________________________________________________
#   BEGIN PROCESSING                                                        ####


number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

main.processing = function(){  


    extract.sequence.w.flanking.regions = function(genome.number){
      #read the configuration file

      config.variables = multi.str.split(config.file, "=", 1)
        #begin by parsing the config file for the genome name, fasta file path and blastdb path
        config.settings = config.file[grep(genome.number, config.variables)]
        
        config.settings.temp = config.settings[[1]]
        genome.name = strsplit(config.settings.temp[1], "=")
        genome.name = genome.name[[1]][2]
        
        config.settings.temp = config.settings[[2]]
        fa.path = strsplit(config.settings.temp[1], "=")
        fa.path = fa.path[[1]][2]
        
        config.settings.temp = config.settings[[3]]
        blastdb.path = strsplit(config.settings.temp[1], "=")
        blastdb.path = blastdb.path[[1]][2]
      
        blast.files = list.files(p("jobs/", gene.name, "/blast.results"), pattern = paste0(genome.number, '.*?.blast'))
        
        blastdf1 = tryCatch(read.blast(p("jobs/", gene.name, "/blast.results/", blast.files[1])), error = function(e){
          write('No BLAST hits for this sequence', p('jobs/', gene.name, '/primers/error.txt'))
          stop('No BLAST hits for this sequence')
        })

        
        blastdf1$sstart.mb = blastdf1$sstart / 1000000
        blastdf1$send.mb = blastdf1$send / 1000000
        blastdf_orig = read.blast(p("jobs/", gene.name, "/blast.results/", blast.files[1]))
          
        #read the fasta index for this particular genome
        fasta.index1 = read.csv(p("./fasta.indexes/", genome.name, ".fa.idx"), stringsAsFactors = F, header = T) 

        print('opt$cds.max.intron.size')
        print(opt$cds.max.intron.size)

        blastdf1.parsed_orig = parse.scaffold.blast(blastdf1, opt$cds.max.intron.size)[[1]]
        blastdf1.parsed = parse.scaffold.blast(blastdf1, opt$cds.max.intron.size)[[1]]
        
        original.scaf.names = multi.str.split(blastdf1.parsed$scaffold, ".$!", 1)   
        
        fasta.index1$offset = as.numeric(fasta.index1$offset)

        genome.assembly.subset.genomic.match <<- readDNAStringSet(fasta.index1[match(original.scaf.names[1], multi.str.split(fasta.index1$desc, " ", 1)), ])  
        names(genome.assembly.subset.genomic.match) = multi.str.split(names(genome.assembly.subset.genomic.match), " ", 1)
        template_sequence_genomic = extract.sequence(genome.assembly.subset.genomic.match, blastdf1.parsed[1, ], 1, opt$start.buffer, opt$end.buffer)  

        opt$fasta.path <<- p("jobs/", gene.name, "/seq/extended/seqs/input_w_flanking.fa")
        query.fa.path = p("jobs/", gene.name, "/seq/extended/seqs/input_w_flanking.fa")
        writeXStringSet(template_sequence_genomic, p("jobs/", gene.name, "/seq/extended/seqs/input_w_flanking.fa"))    
        
        input_sequence <<- readDNAStringSet(p("jobs/", gene.name, "/seq/extended/seqs/input_seq.fa"))
        
        #SECOND BLAST WITH FLANKING REGIONS OF INPUT SEQUENCE INCLUDED 
        source("scripts/perform.blast.rscript.R")
      
    }

    extract.homologues = function(genome.number, number.to.extract, setup.dialign.anchors){
        if(missing(setup.dialign.anchors)) setup.dialign.anchors = T
        if(missing(number.to.extract)) number.to.extract = 'all'

        config.variables = multi.str.split(config.file, "=", 1)
        config.settings = config.file[grep(genome.number, config.variables)]
        
        config.settings.temp = config.settings[[1]]
        genome.name = strsplit(config.settings.temp[1], "=")
        genome.name = genome.name[[1]][2]

        #read the fasta index for this particular genome
        fasta.index1 = read.csv(p("./fasta.indexes/", genome.name, ".fa.idx"), stringsAsFactors = F, header = T)        
        
        #### ITERATION TWO ####
        blast.files = list.files(p("jobs/", gene.name, "/blast.results"), pattern = paste0(genome.number, '.*?.blast'))
        
        blastdf0 = read.blast(p("jobs/", gene.name, "/blast.results/", blast.files[grep(paste0(genome.number, ".*?w_flanking"), blast.files)]))
          
        #SEQUENCE EXTRACTION WITH FULL TEMPLATE (INCLUDING FLANKING REGIONS)
        blastdf1.parsed = parse.scaffold.blast(blastdf0, opt$cds.max.intron.size)  
        
		fasta.index1$offset = as.numeric(fasta.index1$offset)
        genome.assembly.subset.genomic.match <<- readDNAStringSet(fasta.index1[match(unique(blastdf1.parsed[[1]]$scaffold), multi.str.split(fasta.index1$desc, " ", 1)), ])
        names(genome.assembly.subset.genomic.match) = multi.str.split(names(genome.assembly.subset.genomic.match), " ", 1)
        sequences = DNAStringSet()  
        if(number.to.extract == 'all') number.to.extract = nrow(blastdf1.parsed[[1]])       


        #SEQUENCE EXTRACTION
        for(i in 1:number.to.extract){                 
          #MASKING OF INTER-HSP DISTANCES WITHIN THE SAME GROUP WITH Ns      
          rev.comp = blastdf1.parsed[[1]][i, ]$rev.comp
          rchr = blastdf1.parsed[[1]][i, ]$scaffold
            temp.df = blastdf1.parsed[[2]][which(blastdf1.parsed[[2]]$sseqid == blastdf1.parsed[[1]]$groupid[i]), ]
            chr = blastdf1.parsed[[1]]$scaffold[i]  
            #remove any _ concatenations that distinguished groups in orientation check
            # chr = strsplit(chr, "_")
            # chr = chr[[1]][1]    

            if(nrow(temp.df) == 1){      
              #if only 1 HSP, just add it to the list of sequences
              if(rev.comp == F) sequences = c(sequences, DNAStringSet(genome.assembly.subset.genomic.match[[chr]][temp.df$sstart[1]:temp.df$send[1]]))
              if(rev.comp == T) sequences = c(sequences, DNAStringSet(reverseComplement(genome.assembly.subset.genomic.match[[chr]][temp.df$sstart[1]:temp.df$send[1]])))
            } else {
              if(opt$mask.inter.hsp.distances == F | i == 1){
                #extract sequences without masking inter HSP distances with Ns if this option is false
                if(rev.comp == F) sequences = c(sequences, DNAStringSet(genome.assembly.subset.genomic.match[[chr]][min(temp.df$sstart):max(temp.df$send)]))
                if(rev.comp == T) sequences = c(sequences, DNAStringSet(reverseComplement(genome.assembly.subset.genomic.match[[chr]][min(temp.df$sstart):max(temp.df$send)])))
              } else {

              temp.df = temp.df[sort(temp.df$sstart, index.return = T)$ix, ]

              if(rev.comp == T) temp.df = temp.df[sort(temp.df$sstart, index.return = T, decreasing = T)$ix, ]

              group.subsequences = DNAStringSet()
              subseq.differences = as.numeric()        
              for(x in 1:nrow(temp.df)){
                if(all(temp.df$sstart < temp.df$send)){
                  if(nrow(temp.df) == 0) browser()
                  if(rev.comp == F) group.subsequences = c(group.subsequences, DNAStringSet(genome.assembly.subset.genomic.match[[chr]][temp.df$sstart[x]:temp.df$send[x]]))
                  if(rev.comp == T) group.subsequences = c(group.subsequences, DNAStringSet(reverseComplement(genome.assembly.subset.genomic.match[[chr]][temp.df$sstart[x]:temp.df$send[x]])))
                }     

                if(x != nrow(temp.df)){
                  #calculate how far apart the HSPs are
                  subseq.differences = c(subseq.differences, abs(temp.df$send[x] - temp.df$sstart[(x + 1)]))
                }
              }
              
              #combine HSPs with interleaving regions masked by Ns
              subseq.gaps = lapply(subseq.differences, function(x) DNAStringSet(DNAString(paste(rep("N", 100), collapse = ""))))
              subsequences.w.gaps = DNAStringSet()
              for(x in 1:length(group.subsequences)){            
                subsequences.w.gaps = c(subsequences.w.gaps, group.subsequences[x])
                if(x != length(group.subsequences)) subsequences.w.gaps = c(subsequences.w.gaps, subseq.gaps[[x]])
              }

              masked.subsequence = DNAStringSet(DNAString(do.call(paste0, lapply(subsequences.w.gaps, as.character)))) 
              sequences = c(sequences, masked.subsequence)           
            }
            }
        }

        

        if(nrow(blastdf1.parsed[[1]]) <= 1){
          write('No homologues found', p("jobs/", gene.name, "/error.txt"))
          stop('No homologues found')    
        }
        
        names(sequences) = paste0(blastdf1.parsed[[1]]$scaffold, "_", blastdf1.parsed[[1]]$query.start)[1:number.to.extract]  
        

        
        # sequences = c(DNAStringSet(input_sequence), sequences)

        if(setup.dialign.anchors == T){
          #SETUP ANCHOR POINTS FOR DIALIGN
          coord.query.start = c(1, sort(blastdf1.parsed[[1]]$query.start[2:number.to.extract], index.return = T)$ix + 1) #order the sequences by query start position
          blastdf1.parsed[[1]] = blastdf1.parsed[[1]][coord.query.start, ]
        
          sequences = sequences[coord.query.start]  

          dialign.df1 = blastdf1.parsed[[1]]
          dialign.df1$dialign1 = 1 #position of the first sequence to be anchored
          dialign.df1$dialign2 = 1:nrow(dialign.df1) #position of the second sequence to be anchored
          dialign.df1$dialign3 = dialign.df1$query.start #beginning position of the anchor point in sequence 1
          dialign.df1$dialign4 = 1 #beginning position of the anchor point in sequence 2
          dialign.df1$dialign5 = 5 #length of anchor
          dialign.df1$dialign6 = 20 #anchor priority
          dialign.df1 = dialign.df1[-1, ]

          dialign.df1 = dialign.df1[, grep('dialign', colnames(dialign.df1))]
          new_dialign_anchors = dialign.df1

          new_dialign_anchors[, 1] = new_dialign_anchors[, 1]
          new_dialign_anchors[, 2] = new_dialign_anchors[, 2]

        } else {
          new_dialign_anchors = 'No anchors'
        }

        # system(p("scripts/run.dialign.sh jobs/", gene.name, "/"))
        return(list(sequences, new_dialign_anchors))
    }


    sequences = lapply(1:number.genomes, function(genome.number){
        if(genome.number == 1){
          extract.sequence.w.flanking.regions(genome.number)
          return(extract.homologues(genome.number))
        } else {          
          return(extract.homologues(genome.number, number.to.extract = 1, setup.dialign.anchors = F))
        }
    })    

    if(number.genomes > 1){
        #combine genome sequences from the other genomes
        other.genome.seq = lapply(sequences[2:number.genomes], function(x){
          x[[1]]
        })

        other.genome.seq = do.call(c, other.genome.seq)

        all.seq = c(input_sequence, sequences[[1]][[1]][1], other.genome.seq, sequences[[1]][[1]][2:length(sequences[[1]][[1]])])

        dialign_anc = sequences[[1]][[2]]
        dialign_anc$dialign1 = dialign_anc$dialign1 + 1
        dialign_anc$dialign2 = dialign_anc$dialign2 + length(other.genome.seq) + 1
    } else {
        all.seq = c(input_sequence, sequences[[1]][[1]])
        dialign_anc = sequences[[1]][[2]]
        dialign_anc$dialign1 = dialign_anc$dialign1 + 1
        dialign_anc$dialign2 = dialign_anc$dialign2 + 1
    }


    write.table(dialign_anc, p("jobs/", gene.name, "/seq/extended/seqs/all.anc"), quote = F, sep = " ", col.names = F, row.names = F)
    writeXStringSet(all.seq, p("jobs/", gene.name, "/seq/extended/seqs/all.fa"))
    
}


main.processing()
