library(dplyr)
options(stringsAsFactors = F)

options(show.error.locations = T)

options(error = function() { traceback(2); if(!interactive()) quit('no', status = 1, runLast = FALSE) })

##### FUNCTIONS ######

p = function(...){
  #quicker paste function. must supply vector of 
  #elements to paste together (i.e. with c())
  paste(..., collapse = "", sep = "")
}

multi.str.split = function(x, character.split, split.id){
  #An implementation of strsplit() that allows selection of a particular element of the output value of strsplit(), 
  #for all elements of the input vector to which the split is applied.
  #--------------
  #args:
  #x = character vector
  #character.split = character specifying what character to split with 
  #split.id = integer specifying which element of the split you want
  
  unlist(lapply(x, function(q) strsplit(q, character.split)[[1]][split.id]))
}


add.blast.colnames = function(blastdf){
  colnames(blastdf)[1:12] = c("qseqid", "sseqid", "percentage.identical", "alignment.length", "no.mismatch",
                              "no.gap.openings", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  return(blastdf)
}


read.blast = function(filepath){
  g = read.table(filepath)
  g = add.blast.colnames(g)
  g[, 3:12] = lapply(g[, 3:12], as.numeric) #convert number columns to numeric
  return(g)
}

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

#make a new dataframe
newdf = function(..., no.rows){
  #...: an unlimited number of column names
  #no.rows: boolean value, should the dataframe contain zero rows?
  if(missing(no.rows)) no.rows = F
  df=as.data.frame(matrix(nrow=1, ncol=length(c(...))))
  df[is.na(df)]=""
  colnames(df)=c(...)
  
  if(no.rows == T){
    df = df[-1, ]
  }
  
  return(df)
}

convert.to.character.data.frame = function(df){
  df = as.data.frame(df)
  g = df
  g[] = lapply(df, as.character)
  return(g)
}

make.counter = function(){
  counter = 0
  ret.count = function(){
    counter <<- counter + 1
    counter
  }
}

combine.list.of.data.frames = function(list.name){
  g = newdf(colnames(list.name[[1]]), no.rows = T)
  for(i in 1:length(list.name)){
    g = rbind(g, list.name[[i]])
  }
  return(g)
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
  # browser()
  identi.coord = which.max((potential.homeologues$homo_length / input_sequence) * (potential.homeologues$avg.percent.identical / 100)) 
  g = 1:nrow(potential.homeologues)
  g = g[-identi.coord]  
  potential.homeologues = potential.homeologues[c(identi.coord, g), ] 
  
  if(input_sequence > 1500){
    #this will remove all blast hits for small sequences. need an if statement
    coord_to_rm = which(potential.homeologues$length < 500)
    if(length(coord_to_rm) != nrow(potential.homeologues)){
      if(length(coord_to_rm) > 0) potential.homeologues = potential.homeologues[-coord_to_rm, ]           
    }
    
  }

  
#   existing.homeologue.files = grep('potential_homeologues', list.files(p('jobs/', gene.name, '/blast.results/')))
#   if(length(existing.homeologue.files) == 0){
#     write.csv(potential.homeologues, p('jobs/', gene.name, '/blast.results/potential_homeologues1.csv'), row.names = F)
#   } else {
#     write.csv(potential.homeologues, p('jobs/', gene.name, '/blast.results/potential_homeologues', (length(existing.homeologue.files) + 1), '.csv'), row.names = F)
#   }



  list(potential.homeologues, blastdf1)
}



extract.sequences = function(blastdf1.parsed, orig.blastdf, genome.assembly.subset.genomic.match, counter){    
    opt = list()
    opt$mask.inter.hsp.distances = F
    sequences = DNAStringSet()  
  
    #SEQUENCE EXTRACTION
    for(i in 1:nrow(blastdf1.parsed)){       
    #MASKING OF INTER-HSP DISTANCES WITHIN THE SAME GROUP WITH Ns       
    rev.comp = blastdf1.parsed[i, ]$rev.comp
    rchr = blastdf1.parsed[i, ]$scaffold
        temp.df = orig.blastdf[which(orig.blastdf$sseqid == blastdf1.parsed$groupid[i]), ]
        chr = blastdf1.parsed$scaffold[i]
        
        #remove any _ concatenations that distinguished groups in orientation check
        chr = strsplit(chr, "_")
        chr = chr[[1]][1]
        
        if(nrow(temp.df) == 1){      
        #if only 1 HSP, just add it to the list of sequences
        if(rev.comp == F) sequences = c(sequences, DNAStringSet(genome.assembly.subset.genomic.match[[chr]][temp.df$sstart[1]:temp.df$send[1]]))
        if(rev.comp == T) sequences = c(sequences, DNAStringSet(reverseComplement(genome.assembly.subset.genomic.match[[chr]][temp.df$sstart[1]:temp.df$send[1]])))
        } else {
        if(opt$mask.inter.hsp.distances == F | i == 1 ){ # don't perform masking if the sequence is the template / input (i == 1)
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
            tryCatch(group.subsequences[x], error = function(e) browser())
            subsequences.w.gaps = c(subsequences.w.gaps, group.subsequences[x])
            if(x != length(group.subsequences)) subsequences.w.gaps = c(subsequences.w.gaps, subseq.gaps[[x]])
        }

        masked.subsequence = DNAStringSet(DNAString(do.call(paste0, lapply(subsequences.w.gaps, as.character)))) 
        sequences = c(sequences, masked.subsequence)           
        }
        }
    }

    names(sequences) = paste0(blastdf1.parsed$scaffold, "_", blastdf1.parsed$query.start)  

    #SETUP ANCHOR POINTS FOR DIALIGN
    coord.query.start = c(1, sort(blastdf1.parsed$query.start[2:nrow(blastdf1.parsed)], index.return = T)$ix + 1)
    blastdf1.parsed = blastdf1.parsed[coord.query.start, ]
    sequences = sequences[coord.query.start]

    # sequences = c(DNAStringSet(input_sequence), sequences)   

    dialign.df1 = blastdf1.parsed
    dialign.df1$dialign1 = 1 #position of the first sequence to be anchored
    dialign.df1$dialign2 = 1:nrow(dialign.df1) #position of the second sequence to be anchored
    dialign.df1$dialign3 = dialign.df1$query.start #beginning position of the anchor point in sequence 1
    dialign.df1$dialign4 = 1 #beginning position of the anchor point in sequence 2
    dialign.df1$dialign5 = 5 #length of anchor
    dialign.df1$dialign6 = 20 #anchor priority
    dialign.df1 = dialign.df1[-1, ]

    # dialign.df2 = blastdf1.parsed
    # dialign.df2$dialign1 = 1  
    # dialign.df2$dialign2 = 1:nrow(dialign.df2)
    # dialign.df2$dialign3 = dialign.df2$query.end
    # dialign.df2$dialign4 = (dialign.df2$length - 3)
    # dialign.df2$dialign5 = 3
    # dialign.df2$dialign6 = 10
    # dialign.df2 = dialign.df2[-1, ]  

    # dialign_anchors = rbind(dialign.df1, dialign.df2)  
    # dialign_anchors = dialign_anchors[, grep('dialign', colnames(dialign_anchors))]

    dialign.df1 = dialign.df1[, grep('dialign', colnames(dialign.df1))]
    new_dialign_anchors = dialign.df1
    # new_dialign_anchors = matrix(nrow = nrow(dialign_anchors) + 1, ncol = ncol(dialign_anchors))
    # new_dialign_anchors[1:nrow(dialign_anchors), ] = as.matrix(dialign_anchors)
    # ndar = nrow(new_dialign_anchors)
    # ndac = ncol(new_dialign_anchors)   
    new_dialign_anchors[, 1] = new_dialign_anchors[, 1] + 1
    new_dialign_anchors[, 2] = new_dialign_anchors[, 2] + 1

    return(list(sequences, new_dialign_anchors))

}

genomic.genes = readDNAStringSet('~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/flanking.genomic.genes/all.genomic.gene.sequences.w.flanking.regions.1kb.fa')
# genomic.gene.lengths = data.frame(names(genomic.genes), width(genomic.genes))
# colnames(genomic.gene.lengths) = c('name', 'length')

# write.csv(genomic.gene.lengths, '~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/flanking.genomic.genes/genomic.gene.lengths.csv', row.names = F)
genomic.gene.lengths = read.csv('~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/flanking.genomic.genes/genomic.gene.lengths.csv')


g = read.delim('~/project.phd.main/bioinf/blast/probe.vs.genome.blast/results.blast/all.genomic.genes.w.flanking.blast', sep = '\t', header = F, stringsAsFactors = F)

colnames(g) = c("qseqid", "sseqid", "percentage.identical", "alignment.length", "no.mismatch",
                              "no.gap.openings", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

#these BLAST hits are in a different order to the genes and need to be reordered
g2 = split(g, factor(g$qseqid, levels = unique(genomic.gene.lengths$name)))

#BLAST SORTING
g3 = g2[sort(genomic.gene.lengths$length, index.return = T, decreasing = T)$ix]
g4 = g3[which(genomic.gene.lengths.sorted$length < 10000)]

#GENE LENGTH SORTING
genomic.gene.lengths.sorted = genomic.gene.lengths[sort(genomic.gene.lengths$length, index.return = T, decreasing = T)$ix, ]
genomic.gene.lengths.sorted2 = genomic.gene.lengths.sorted[which(genomic.gene.lengths.sorted$length < 10000),]

#GENE SEQUENCE SORTING
g.genes.sort = genomic.genes[sort(genomic.gene.lengths$length, index.return = T, decreasing = T)$ix]
g.genes.sort = g.genes.sort[which(genomic.gene.lengths.sorted$length < 10000)]





library(dplyr)
library(parallel)

all.parsed = mclapply(g4, function(x){        
    input_sequence <<- genomic.gene.lengths$length[which(genomic.gene.lengths$name == x$qseqid[1])]
    tryCatch(parse.scaffold.blast(x, 4000), error = function(e) 'error')
}, mc.cores = 50)


right.lengths = sapply(all.parsed, function(x){    
    qlength = as.numeric(x[[1]]$end[1]) - as.numeric(x[[1]]$start[1])
    slength = as.numeric(x[[1]]$query.end[1]) - as.numeric(x[[1]]$query.start[1])
    qlength == slength
})

num.correct = unlist(Map(function(x, y){
    x[[1]]$query.start[1] == 1 & x[[1]]$query.end[1] == y
}, all.parsed, genomic.gene.lengths.sorted2$length))
    
all.parsed.correct = all.parsed[which(right.lengths & num.correct)]
g.gene.length.sort3 = genomic.gene.lengths.sorted2[which(right.lengths & num.correct), ]
g.genes.sort = g.genes.sort[which(right.lengths & num.correct)]
g5 = g4[which(right.lengths & num.correct)]

library(Biostrings)
genome = readDNAStringSet('~/project.phd.main/genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta')


t1 = extract.sequences(all.parsed.correct[[2]][[1]], all.parsed.correct[[2]][[2]], genome)

ptm <- proc.time()
counter = 1
all.sequences = mclapply(all.parsed.correct, function(x){
    print(counter)
    g = tryCatch(extract.sequences(x[[1]], x[[2]], genome, counter), error = function(e) "error")
    counter <<- counter + 1
    g
}, mc.cores = 20)
proc.time() - ptm


#check lengths of extracted sequences match lengths of inputted sequences
checkextract = sapply(all.sequences, function(x){
    width(x[[1]])[1]
}) == width(g.genes.sort)

all.sequences.worked = all.sequences[-which(checkextract != T)]
g.gene.length.sort4 = g.gene.length.sort3[-which(checkextract != T), ]
all.parsed.worked = all.parsed.correct[-which(checkextract != T)]
g6 = g5[-which(checkextract != T)]

ptm <- proc.time()
all.sequences.w.input = mclapply(all.sequences.worked, function(x){      
    g = tryCatch(DNAStringSet(x[[1]][[1]][1001:(length(x[[1]][[1]]) - 1000)]), error = function(e) "no work")
    if(class(g) == "DNAStringSet"){
        x[[1]] = c(g, x[[1]])
        names(x[[1]])[1] = "input_sequence"
    } 

    x
}, mc.cores = 30)
proc.time() - ptm

test.worked = sapply(all.sequences.w.input, function(x) class(x[[1]]) == 'DNAStringSet')
length(which(test.worked))

save(all.sequences, file = "~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/flanking.genomic.genes/all.sequences.w.homologues.robject")

ptm <- proc.time()
counter = 1
all.sequences.w.input_named = lapply(all.sequences.w.input, function(x){    
    if(class(x[[1]]) == "DNAStringSet"){
        names(x[[1]])[1] = paste0(names(x[[1]])[1], "_", counter)
        names(x[[1]])[2] = paste0(names(x[[1]])[2], "_", g.gene.length.sort3$name[counter])
    } 
    counter <<- counter + 1
    x
})
proc.time() - ptm


lapply(all.sequences.w.input_named, function(x) names(x[[1]])[1:2])[1:20000]


coords_start = seq(1, length(all.sequences.w.input), round(length(all.sequences.w.input) / 500))
coords_end = coords_start - 1
coords_end = c(coords_end[-1],  length(all.sequences.w.input))


alignment.dir = '~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/all.gene.alignments'
# dir.create(alignment.dir)


counter = 1
Map(function(x, y){
    dir.create(paste0(alignment.dir, '/set', counter))
    lapply(all.sequences.w.input_named[x:y], function(z){
        if(class(z[[1]]) == "DNAStringSet"){
            sequence.num = strsplit(names(z[[1]])[1], "_")
            sequence.num = sequence.num[[1]][3]
            sequence.name = names(z[[1]])[2]
            
            #write msa
            writeXStringSet(z[[1]], paste0(alignment.dir, '/set', counter, "/", sequence.num, ".", sequence.name, ".fa"))
            #write anchors
            write.table(z[[2]], paste0(alignment.dir, '/set', counter, "/", sequence.num, ".", sequence.name, ".anc"), quote = F, sep = " ", col.names = F, row.names = F)
        }       
    })
    print(paste0("done ", counter))
    counter <<- counter + 1
    
}, coords_start, coords_end)

setwd('/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/IWGSC/all.alignments.dialign/')

all.sets = list.files(pattern = 'set')
all.files = lapply(all.sets, function(x){
    list.files(paste0(x, '/'), pattern = 'Traes')
})

all.files.full = lapply(all.sets, function(x){
    list.files(paste0(x), pattern = 'Traes', full.names = T)
})


all.numbers = lapply(all.files, function(x){
    as.numeric(multi.str.split(x, "\\.", 1))
})


all.names1 = lapply(all.numbers, function(x){
    g.gene.length.sort4$name[x]
})

length(unlist(all.names1)) == length(unlist(all.files.full))

Map(function(x, y){
    blast.chromo = substr(multi.str.split(multi.str.split(x, "\\.", 2), "_", 1), 4, 5)
    named.chromo = substr(multi.str.split(x, "TraesCS", 2), 1, 2)
    fixed.chromo = substr(multi.str.split(y, "TraesCS", 2), 1, 2)
    fixed.chromo[which(fixed.chromo == "U0")] = 'Un'
    
    blast.chromo[which(blast.chromo != fixed.chromo)]
    fixed.chromo[which(blast.chromo != fixed.chromo)]
    x[which(blast.chromo != fixed.chromo)]
    y[which(blast.chromo != fixed.chromo)]

    length(which(blast.chromo != fixed.chromo))
    head(x)
    
    new.paths = paste0('../all.alignments.fixed/', multi.str.split(x, "_", 1), '_', multi.str.split(x, "_", 2), '_', y)

    lapply(new.paths, dir.create)

    Map(function(z, z1){
        fn1 = list.files(z, full.names = T)        
        file.copy(fn1, z1)
    }, x, new.paths)    
    
    print('done1')


}, all.files.full, all.names1)


setwd('../all.alignments.fixed')
all.sets.new = list.files(pattern = 'set')
all.files = lapply(all.sets.new, function(x){
    list.files(paste0(x, '/'), pattern = 'Traes')
})
all.files.new = lapply(all.sets.new, function(x){
    list.files(paste0(x), full.names = T)
})




all.numbers.new = lapply(all.files, function(x){
    as.numeric(multi.str.split(x, "\\.", 1))
})


all.names1 = lapply(all.numbers, function(x){
    g.gene.length.sort4$name[x]
})

lapply(all.files.new, function(x){
    numbers1 = as.numeric(multi.str.split(multi.str.split(x, "/", 2), "\\.", 1))
    Map(function(z, z1){
        # write.csv(g6[[z1]], paste0(z, '/blast.results.csv'), row.names = F)
        write.csv(all.parsed.worked[[z1]][[1]], paste0(z, '/potential.homologues.csv'), row.names = F)
    }, x, numbers1)
    print('done1')
})

all.parsed.worked[[100]][[1]]
length(all.parsed.worked)
length(g6)


library(parallel)
library(Biostrings)
mclapply(all.files.new, function(x){
    lapply(x, function(y){
        gene.name = multi.str.split(y, "_", 3)
        if(file.exists(paste0(y, '/align.w.primers.fa'))){
            align = readDNAStringSet(paste0(y, '/align.w.primers.fa'))
            old.name = strsplit(names(align)[2], "_")[[1]]
            old.name = old.name[-length(old.name)]
            names(align)[2] = paste(c(old.name, gene.name), collapse = "_")
            writeXStringSet(align, paste0(y, '/align.name.corrected.fa'))    
        }
        
    }) 
    
}, mc.cores = 64)


#write blast files 
lapply(all.numbers, function(x){
    browser()
    g6[x]
})