library(Biostrings)
setwd("/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/IWGSC/all.alignments.dialign")

interleave <- function(v1,v2)                                                                                    
{                                                                                                              
	ord1 <- 2*(1:length(v1))-1
	ord2 <- 2*(1:length(v2))
	c(v1,v2)[order(c(ord1,ord2))]
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

conv.matrix.dnastringset = function(mult.mat, row.names1){
  #converts a matrix to a DNAStringSet object
  #args:
  #	mult.mat - DNAMultipleAlignment object that has been converted to character df w/ convert.to.character.data.frame(as.matrix(dnamultiplealignment.object))
  #	row.names1 - Character vector of names for the sequences
  
  mult.mat3 = apply(mult.mat, 1, function(x) paste(x, collapse = ""))
  
  #convert each object in the list to a DNAStringSet object
  mult.mat4 = lapply(mult.mat3, DNAStringSet)
  
  #make new DNAStringSet object
  mult.mat5 = DNAStringSet()
  
  #add our sequences to the new DNAStringSet object one by one
  for(i in 1:length(mult.mat4)){
    mult.mat5 = c(mult.mat5, mult.mat4[[i]])
  }
  
  names(mult.mat5) = row.names1
  
  mult.mat5
}


calculate.primers = function(alignment.file, set.folder){
    gene.name = "gene1"
    full.gene.product = F
    min.product.size = 800
    max.product.size = 2000
    start.buffer = 1000
    end.buffer = 1000

    #the base path for auto primer picker
    project.path = paste0('/home/ac14037/project.phd.main/rotation1scripts_v4/original.data/IWGSC/all.alignments.dialign/', set.folder)

    
    #parse number of genomes in configuration file
    number.genomes = 1        
    gene.name = strsplit(alignment.file, '/')[[1]]
    gene.name = gene.name[length(gene.name)]    
    mult.align_trim = readDNAStringSet(paste0(alignment.file, '/', gene.name, '.fa.dialign.fa'))
    mult.align1 = mult.align_trim

    main.gene.row = 1
    template.row = 2

    homologue.rows = 3:length(mult.align1)
    # variety.rows = (max(homologue.rows) + 1):nrow(mult.align1)



    #   ____________________________________________________________________________
    #   DEFINE FUNCTIONS                                                        ####

    conv.mult.align.dnastringset = function(mult.align){
        #converts a multiplealignment object to dnastringset
        mult.mat = as.matrix(mult.align)
        
        mult.mat3 = apply(mult.mat, 1, function(x) paste(x, collapse = ""))
        
        #convert each object in the list to a DNAStringSet object
        mult.mat4 = lapply(mult.mat3, DNAStringSet)
        
        #make new DNAStringSet object
        mult.mat5 = DNAStringSet()
        
        #add our sequences to the new DNAStringSet object one by one
        for(i in 1:length(mult.mat4)){
            mult.mat5 = c(mult.mat5, mult.mat4[[i]])
        }
        
        names(mult.mat5) = rownames(mult.align)
        
        mult.mat5
        }

        remove.inserts = function(dna.string){
        #removes inserted bases ("-") from a DNAString object
        raw.seq = as.character(dna.string)
        
        #following code removes "-" from raw.seq
        raw.seq2 = strsplit(raw.seq, "")
        raw.seq2 = raw.seq2[[1]]
        raw.seq3 = raw.seq2[!raw.seq2 == "-"]
        raw.seq4 = paste(raw.seq3, collapse = "")
        
        DNAString(raw.seq4)
    }

    calculate.start.and.end.ranges = function(sequence.to.cut, product.size){
        #calculate coordinates for cutting of sequence into product sizes of ~ 700 bp
        #args:
        # sequence.to.cut - a DNAString object
        # product.size - integer specifying the desired length of the PCR product in bases
        
        if(missing(product.size)) product.size = 700
        
        if(product.size > length(sequence.to.cut)){
            start.range = 1
            end.range = length(sequence.to.cut)
        } else {
            num.divisions = ceiling(length(sequence.to.cut) / product.size)
            iterations = ceiling(length(sequence.to.cut) / num.divisions)
            
            start.range = seq(1, length(sequence.to.cut) + 100, iterations)
            end.range = seq(1, length(sequence.to.cut) + 100, iterations)
            
            start.range = start.range[-length(start.range)]
            end.range = end.range[-1] - 1
            
            end.range[length(end.range)] = length(sequence.to.cut)
            
            start.range[2:length(start.range)] = start.range[2:length(start.range)] - 100
            end.range[1:(length(end.range)-1)] = end.range[1:(length(end.range)-1)] + 100
        }
        
        ranges1 = list(start.range, end.range)
        names(ranges1) = c("start.range", "end.range")
        ranges1  
    }

    examine.primer.homologue.concordance = function(left.primer, right.primer, template.sequence){
        #given two primer sequences, checks how well the primers complement a given template sequence. 
        #returns the percentage sequence match against the template for each primer
        #args:
        # left.primer - character string of left primer sequence
        # right.primer - character string of right primer sequence
        left.align = pairwiseAlignment(DNAString(left.primer), template.sequence)
        right.align = pairwiseAlignment(reverseComplement(DNAString(right.primer)), template.sequence)
        
        alignment.pids = list(pid(left.align), pid(right.align))
        names(alignment.pids) = c("left.align.pid", "right.align.pid")
        alignment.pids
        }

        grab.primer.seq = function(primer.id, left.or.right, primer3.output){
        #extracts primer sequence from primer3 output file
        #args:
        # primer.id - integer (starting at 0)
        # left.or.right - character string; either "LEFT" or "RIGHT"
        primer1 = primer3.output[grep(paste0("PRIMER_", left.or.right, "_", primer.id, "_SEQUENCE"), primer3.output)]
        primer1 = strsplit(primer1, "=")
        primer1 = primer1[[1]][2]
        primer1
    }


    #   ____________________________________________________________________________
    #   DETECT SNPS                                                             ####


    grab.homeologous.snps_new = function(input.row, template.row, homologue.rows, multiple.alignment,
    perform.masking, mask.bin.size, mask.threshold, allow.hyphens.in.mask, allow.hyphens.for.snp.detection){
        #gets homeologous snps when there is only one genome
        #takes a DNAMultipleAlignment object and returns a numeric vector of the column coordinates containing homeologous SNPs
        #args:
        # input.row - Integer; the row of the sequence inputted by the user (usually 1)
        #template.row - Integer; the row of the sequence to design primers from (usually 2)
        # homologue.rows - numeric vector containing the row coordinates of the homologous sequences (either paralogous or homeologous)
        # multiple.alignment - DNAMultipleAlignment class 
        # perform.masking - Boolean, indicates whether sequences of low similarity to template should be masked
        #allow.hyphens.for.snp.detection - Boolean, indicates whether hyphens in homologues (not the template sequence) will be considered SNPs
            if(missing(perform.masking)) perform.masking = F
            if(missing(mask.bin.size)) mask.bin.size = 10
            if(missing(mask.threshold)) mask.threshold = 5
            if(missing(allow.hyphens.in.mask)) allow.hyphens.in.mask = T
            if(missing(allow.hyphens.for.snp.detection)) allow.hyphens.for.snp.detection = T
        
        #NB. Insertions "-" in the template sequence cannot be allowed when classifying SNPs, as these
        #will be subsequently removed by get.coordinates.after.removing.hyphens(), meaning that the
        #program maps the primers to the wrong locations in the final multiple sequence alignment output.
        
        mult.align.mat1 = as.data.frame(as.matrix(multiple.alignment))

        perform.masking.function = function(msa.matrix1, mask.bin.size, mask.threshold, allow.hyphens.in.mask){
            # mask.threshold - integer, underneath what percentage of similarity should masking be performed?
            #            e.g. 40 - when bins have less than 40% nucleotides in common, mask them
                if(missing(mask.bin.size)) mask.bin.size = 10
                if(missing(mask.threshold)) mask.threshold = 40
                if(missing(allow.hyphens.in.mask)) allow.hyphens.in.mask = F
                
                msa.matrix1.orig = msa.matrix1
                start.sequence.bins = seq(1, ncol(msa.matrix1), mask.bin.size)
                end.sequence.bins = c((start.sequence.bins[2:length(start.sequence.bins)] - 1), ncol(msa.matrix1))

                #performing masking of regions with low sequence identity (in bins of 10)
                bin.dfs1 = lapply(3:nrow(msa.matrix1), function(z){ #lapply across rows
                    bin.similarities = unlist(Map(function(x, y){
                        template.bin.seq = msa.matrix1[2, x:y]
                        target.bin.seq = msa.matrix1[z, x:y]

                        if(allow.hyphens.in.mask == F){
                            if((length(which(template.bin.seq == "-")) > 3) | length(which(target.bin.seq == "-")) > 3) return(10) #don't include bins with hyphens in masking
                        }
                        

                        length(which(msa.matrix1[2, x:y] == msa.matrix1[z, x:y]))
                    }, start.sequence.bins, end.sequence.bins))

                    bin.df1 = data.frame(start.sequence.bins, end.sequence.bins, bin.similarities)
                    colnames(bin.df1) = c("sbin", "ebin", "nsim")

                    #mask bins with less than 4 nucleotides in common with the template sequence        
                    mask.threshold2 = (mask.bin.size / 100) * mask.threshold
                    bin.df1 = bin.df1[which(bin.df1$nsim < mask.threshold2), ]

                    unlist(Map(function(x, y){
                            seq(x, y, 1)
                    }, bin.df1$sbin, bin.df1$ebin))            
                })
            
            print('Performing sequence masking')
            #mask regions with low sequence identity
            for(i in 1:length(bin.dfs1)){
                msa.matrix1[(i + 2), bin.dfs1[[i]]] = "U"
            }

            msa.matrix1
            }

        if(perform.masking == T) mult.align.mat1 = perform.masking.function(mult.align.mat1, mask.bin.size, mask.threshold, allow.hyphens.in.mask)
        
        g = lapply(mult.align.mat1, function(x){          
            if(x[2] == "-") return(0) #return 0 if template sequence has an insertion
            
            if(allow.hyphens.for.snp.detection == F){
            if("-" %in% x[3:length(x)]){
                return(0)
            }
            }    

            if(x[template.row] %in% x[homologue.rows]){
            snp = 0
            } else {
            snp = 1
            }  
            

            snp
            })

        unlist(g)    
    }


    get.start.and.end.base = function(multiple.alignment1, main.sequence.row.number){
        #Returns the positions of the start base and the end base of a specified sequence in a multiple alignment
        #args:
        # multiple.alignment1 - a DNAMultipleAlignment object
        # main.sequence.row.number - Integer; row number of the desired sequence to return positions for
        multiple.alignment1 = as.data.frame(as.matrix(multiple.alignment1))
        mult.align.mat.start.base = min(grep("A|C|T|G", multiple.alignment1[main.sequence.row.number, ]))
        mult.align.mat.end.base = max(grep("A|C|T|G", multiple.alignment1[main.sequence.row.number, ]))  
        g = c(mult.align.mat.start.base, mult.align.mat.end.base)
        names(g) = c("start.base", "end.base")
        g
    }

    get.coordinates.after.removing.hyphens = function(mult.align1, sequence.row.number, homologous.snps1, start.base1, end.base1){
        #Returns the coordinates of homologous SNPs, start base and end base for a specified sequence in a multiple alignment file after removing hyphens
        #args:
        #mult.align1 - a DNAMultipleAlignment file
        #sequence.row.number - Integer, the row coordinate of the desired sequence to transform coordinate for in the multiple alignment
        #homologous.snps1 - Numeric vector, obtained using grab.homoelogous.snps()
        #start.base1 - Integer, obtained using get.start.and.end.base()
        #end.base1 - Integer, obtained using get.start.and.end.base()
        mult.align2 = conv.mult.align.dnastringset(mult.align1)  
        #need to transform coordinates of bases after removing hyphens
        seq.mat = as.data.frame(as.matrix(mult.align2[[sequence.row.number]]))

        seq.mat$snps = homologous.snps1

        seq.mat$start.base = ""
        seq.mat$start.base[start.base1] = 1
        seq.mat$end.base = ""
        seq.mat$end.base[end.base1] = 1

        seq.mat2 = seq.mat[-which(seq.mat[, 1] == "-"), ] #remove hyphens
        seq.mat2$new.coords = 1:nrow(seq.mat2)
        snp.coords.after.filter = which(seq.mat2$snps == 1) #these are the coordinates of the homologous snps in this particular sequence after removing hyphens 
        start.coord.after.filter = which(seq.mat2$start.base == 1) #coordinate of the start codon ATG after removing hyphens
        end.coord.after.filter = which(seq.mat2$end.base == 1)

        #primer3 uses coordinates starting from 0. Need to update the SNP coordinates to reflect this
        snp.coords.after.filter = snp.coords.after.filter - 1
        start.coord.after.filter = start.coord.after.filter - 1
        end.coord.after.filter = end.coord.after.filter - 1

        g = list(snp.coords.after.filter, start.coord.after.filter, end.coord.after.filter)
        names(g) = c("snp.coords", "start.coord", "end.coord")
        g
    }

    generate_align_w_snps = function(){
        homologous.snps = grab.homeologous.snps_new(main.gene.row, template.row, homologue.rows, mult.align_trim)  
        
        main.start.end = get.start.and.end.base(mult.align_trim, main.gene.row)

        coords = get.coordinates.after.removing.hyphens(mult.align_trim, template.row, homologous.snps, main.start.end[1], main.start.end[2])

            
        output.list1 = list(coords, homologous.snps, main.start.end)
        names(output.list1) = c("coords", "homologous.snps", "main.start.end")
        output.list1

    }

    detect_snps_wrapper = function(){
        if(number.genomes > 1){
        homologous.snps = grab.homeologous.snps_orig(variety.rows, homologue.rows, mult.align1, F, T)  
        } else {
        homologous.snps = grab.homeologous.snps_new(main.gene.row, template.row, homologue.rows, mult.align1)  
        }

        homologous.snps = grab.homeologous.snps_orig(variety.rows, homologue.rows, mult.align_trim, F, T)  
        main.start.end = get.start.and.end.base(mult.align1, main.gene.row)

        coords = get.coordinates.after.removing.hyphens(mult.align1, template.row, homologous.snps, main.start.end[1], main.start.end[2])

        

        snps.smaller.than.start = length(coords$snp.coords[which(coords$snp.coords < coords$start.coord)])

        snps.bigger.than.end = length(coords$snp.coords[which(coords$snp.coords > coords$end.coord)])

        homologue.names = rownames(mult.align1)[homologue.rows]
        variety.names = rownames(mult.align1)[variety.rows]

        output.list1 = list(coords, homologous.snps, main.start.end)
        names(output.list1) = c("coords", "homologous.snps", "main.start.end")
        output.list1
    }

    auto_remove_sequences_and_detect_snps = function(){
        if(number.genomes > 1){
        homologous.snps = grab.homeologous.snps_orig(variety.rows, homologue.rows, mult.align1, F, T)  
        } else {
        homologous.snps = grab.homeologous.snps_new(main.gene.row, template.row, homologue.rows, mult.align1)  
        }

        main.start.end = get.start.and.end.base(mult.align1, main.gene.row)

        coords = get.coordinates.after.removing.hyphens(mult.align1, template.row, homologous.snps, main.start.end[1], main.start.end[2])


        snps.smaller.than.start = length(coords$snp.coords[which(coords$snp.coords < coords$start.coord)])

        snps.bigger.than.end = length(coords$snp.coords[which(coords$snp.coords > coords$end.coord)])

        homologue.names = rownames(mult.align1)[homologue.rows]
        variety.names = rownames(mult.align1)[variety.rows]

        while.counter = 1
        while(snps.smaller.than.start < 2){ 
        #check if less than 2 coordinates that are smaller than the start sequence
        #if so, one of the sequences included in the multiple alignment could be too divergent.
        #Try examining the multiple alignment, checking all homologue sequences (from IWGSC assembly)
        #for hyphens before the start coordinate, and remove the sequence with the most hyphens

        print("Less than 2 SNPs found before start of sequence, removing sequence with most hyphens before start")

        get.sequence.to.remove = function(multiple.alignment1){    
            #analyses a multiple alignment to check which sequence has the most hyphens before 
            #the start position of the main sequence (row 1). Returns the name of the sequence
            #with the most hyphens as a character string.
            #args:
            # multiple.alignment1 - a DNAMultipleAlignment object
            mult.align1.mat = as.matrix(multiple.alignment1)
            start.base.coord = min(grep("A|G|T|C", mult.align1.mat[1, ]))
            matrix.before.start = mult.align1.mat[, 1:start.base.coord]    
            #get number of insertions (hyphens) in each sequence - exclude input sequence and its closest BLAST hit
            #from the evaulation 
            num.hyphens = apply(matrix.before.start, 1, function(x) length(which(x == "-")))[3:nrow(matrix.before.start)]
            names(num.hyphens[which(num.hyphens == max(num.hyphens))])
        }

        seq.to.rm = get.sequence.to.remove(mult.align1)

        homo.update = which(homologue.names == seq.to.rm)
        var.update = which(variety.names == seq.to.rm)

        if(length(homo.update) > 0){
            homologue.rows = homologue.rows[1:(length(homologue.rows) - length(homo.update))]
            } 
            variety.rows = variety.rows = variety.rows - length(homo.update)
        if(length(var.update) > 0){
            variety.rows = variety.rows[1:(length(variety.rows) - length(var.update))] 
            }


        number.of.removed.sequences = length(seq.to.rm)

        mult.align.1.dnastringset = conv.mult.align.dnastringset(mult.align1)  

        mult.align1.dnastringset2 = mult.align.1.dnastringset[-which(names(mult.align.1.dnastringset) %in% seq.to.rm)]

        writeXStringSet(mult.align.1.dnastringset, paste0(project.path, gene.name, "/seq/extended/alignments/all.align.old", while.counter, ".fa"))

        print(paste0("Alignment: all.align.old", while.counter, ".fa"))
        while.counter = while.counter + 1

        #write multiple alignment without the sequence with the most hyphens
        writeXStringSet(mult.align1.dnastringset2, paste0(project.path, gene.name, "/seq/extended/alignments/all.align.rev.fa"))
        writeXStringSet(mult.align1.dnastringset2, paste0(project.path, gene.name, "/seq/extended/alignments/all.align.rev.fa"))

        mult.align1 = readDNAMultipleAlignment(paste0(project.path, gene.name, "/seq/extended/alignments/all.align.rev.fa"))
        
        #if there is more than one sequence with the same number of hyphens, they will both be removed. if this is the case,
        #we need to assign new homologue rows according to the number of sequences that were removed.
        # homo.end = (nrow(mult.align1) - (number.genomes - number.of.removed.sequences))
        # homologue.rows = 3:homo.end
        #if two varietal sequences are removed at the same time, this can create false assignment of a varietal sequence as a homologue.
        # if(homo.end < 3) homologue.rows = 3

        # variety.rows = (max(homologue.rows) + 1):nrow(mult.align1)

        if(number.genomes > 1){
        homologous.snps = grab.homeologous.snps_orig(variety.rows, homologue.rows, mult.align1, F, T)  
        } else {
            homologous.snps = grab.homeologous.snps_new(main.gene.row, template.row, homologue.rows, mult.align1)  
        }

        main.start.end = get.start.and.end.base(mult.align1, main.gene.row)

        coords = get.coordinates.after.removing.hyphens(mult.align1, template.row, homologous.snps, main.start.end[1], main.start.end[2])


        snps.smaller.than.start = length(coords$snp.coords[which(coords$snp.coords < coords$start.coord)])  
        }

        print("snps.bigger.than.end")
        print(snps.bigger.than.end)

        while.counter = 1
        while(snps.bigger.than.end < 2){ 
        #check if less than 2 coordinates that are smaller than the start sequence
        #if so, one of the sequences included in the multiple alignment could be too divergent.
        #Try examining the multiple alignment, checking all homologue sequences (from IWGSC assembly)
        #for hyphens before the start coordinate, and remove the sequence with the most hyphens

        print("Less than 2 SNPs found after end of sequence, removing sequence with most hyphens after end")

        get.sequence.to.remove = function(multiple.alignment1){    
            #analyses a multiple alignment to check which sequence has the most hyphens before 
            #the start position of the main sequence (row 1). Returns the name of the sequence
            #with the most hyphens as a character string.
            #args:
            # multiple.alignment1 - a DNAMultipleAlignment object
            mult.align1.mat = as.matrix(multiple.alignment1)
            end.base.coord = max(grep("A|G|T|C", mult.align1.mat[1, ]))
            matrix.after.end = mult.align1.mat[, end.base.coord:ncol(mult.align1.mat)]    
            #get number of insertions (hyphens) in each sequence - exclude input sequence and its closest BLAST hit
            #from the evaulation 
            num.hyphens = apply(matrix.after.end, 1, function(x) length(which(x == "-")))[3:nrow(matrix.after.end)]
            names(num.hyphens[which(num.hyphens == max(num.hyphens))])
        }

        seq.to.rm = get.sequence.to.remove(mult.align1)

        homo.update = which(homologue.names == seq.to.rm)
        var.update = which(variety.names == seq.to.rm)

        if(length(homo.update) > 0){
            homologue.rows = homologue.rows[1:(length(homologue.rows) - length(homo.update))]
            } 
            variety.rows = variety.rows = variety.rows - length(homo.update)
        if(length(var.update) > 0){
            variety.rows = variety.rows[1:(length(variety.rows) - length(var.update))] 
            }

        print("seq.to.rm")
        print(seq.to.rm)

        number.of.removed.sequences = length(seq.to.rm)

        mult.align.1.dnastringset = conv.mult.align.dnastringset(mult.align1)  

        mult.align1.dnastringset2 = mult.align.1.dnastringset[-which(names(mult.align.1.dnastringset) %in% seq.to.rm)]

        writeXStringSet(mult.align.1.dnastringset, paste0(project.path, gene.name, "/seq/extended/alignments/all.align.old_end_", while.counter, ".fa"))

        print(paste0("Alignment: all.align.old", while.counter, ".fa"))
        while.counter = while.counter + 1

        #write multiple alignment without the sequence with the most hyphens
        writeXStringSet(mult.align1.dnastringset2, paste0(project.path, gene.name, "/seq/extended/alignments/all.align.rev.fa"))
        writeXStringSet(mult.align1.dnastringset2, paste0(project.path, gene.name, "/seq/extended/alignments/all.align.rev.fa"))

        mult.align1 = readDNAMultipleAlignment(paste0(project.path, gene.name, "/seq/extended/alignments/all.align.rev.fa"))
        
        #if there is more than one sequence with the same number of hyphens, they will both be removed. if this is the case,
        #we need to assign new homologue rows according to the number of sequences that were removed.
        # homo.end = (nrow(mult.align1) - (number.genomes - number.of.removed.sequences))
        # homologue.rows = 3:homo.end
        #if two varietal sequences are removed at the same time, this can create false assignment of a varietal sequence as a homologue.
        # if(homo.end < 3) homologue.rows = 3

        if(number.genomes > 1){
        homologous.snps = grab.homeologous.snps_orig(variety.rows, homologue.rows, mult.align1, F, T)  
        } else {
            homologous.snps = grab.homeologous.snps_new(main.gene.row, template.row, homologue.rows, mult.align1)  
        }

        main.start.end = get.start.and.end.base(mult.align1, main.gene.row)

        coords = get.coordinates.after.removing.hyphens(mult.align1, template.row, homologous.snps, main.start.end[1], main.start.end[2])

        print(coords)

        snps.bigger.than.end = length(coords$snp.coords[which(coords$snp.coords < coords$start.coord)])  
        }
        output.list1 = list(coords, homologous.snps)
        names(output.list1) = c("coords", "homologous.snps")
        output.list1
    }


    #   ____________________________________________________________________________
    #   BEGIN MAIN LOOP                                                         ####

    find.best.primers = function(multiple.alignment, template.sequence.row.number, snp.coords.after.filter, start.coord.after.filter, end.coord.after.filter, product.size.range, span.whole.gene, start.buffer, homologous.snps, coords){
        #Automatically obtains primer sequences
        #args:
        # multiple.alignment - a DNAMultipleAlignment object
        # template.sequence.row.number - Integer; the multiple alignment row of the sequence to use as a template in primer3
        # snp.coords.after.filter - Numeric vector; obtained using grab.homeologous.snps() and then get.coordinates.after.removing.hyphens()
        # start.coord.after.filter - Integer; position of the first base of the start codon after removing hyphens
        # end.coord.after.filter - Integer; position of the final base in the coding sequence after removing hyphens
        # product.size.range - a numeric vector with two elements, the first being the minimum product size, the second the maximum
        # span.whole.gene - Boolean; should the product size span the entire gene with only one set of primers?
        # start.buffer - Integer; how many bases before the start of the gene should be allowed in the product?
        
        if(missing(span.whole.gene)) span.whole.gene = F
        if(missing(start.buffer)) start.buffer = start.coord.after.filter 

        list.best.primer.start.coords = as.numeric()
        list.best.primer.end.coords = as.numeric()
        mult.align2 = conv.mult.align.dnastringset(multiple.alignment)

        #make primers for all possible SNPs, then evaluate
        template.sequence = mult.align2[[template.sequence.row.number]]
        
        template.sequence2 = remove.inserts(template.sequence)
        generate.all.primer.penalites = function(x){        
            
            primer3.forward.input = c(paste0("SEQUENCE_ID=AllForwardPrimers"), 
                        paste0("SEQUENCE_TEMPLATE=", template.sequence2, ""),        
                        "PRIMER_TASK=pick_primer_list",
                        "PRIMER_PICK_RIGHT_PRIMER=0",
                        paste0("PRIMER_NUM_RETURN=", (length(template.sequence2) * 6)),                    
                        "=")
            
            writeLines(primer3.forward.input, paste0(set.folder, '/', gene.name, '/forwardprimers.txt'))

            primer3.reverse.input = c(paste0("SEQUENCE_ID=AllReversePrimers"), 
                paste0("SEQUENCE_TEMPLATE=", template.sequence2, ""),        
                "PRIMER_TASK=pick_primer_list",
                "PRIMER_PICK_LEFT_PRIMER=0",
                paste0("PRIMER_NUM_RETURN=", (length(template.sequence2) * 6)),
                "=")

            writeLines(primer3.reverse.input, paste0(set.folder, '/', gene.name, '/reverseprimers.txt'))
            
            system(paste0('primer3_core < ', set.folder, '/', gene.name, '/forwardprimers.txt > ', set.folder, '/', gene.name, '/forwardprimers.output.txt'))
            system(paste0('primer3_core < ', set.folder, '/', gene.name, '/reverseprimers.txt > ', set.folder, '/', gene.name, '/reverseprimers.output.txt'))
            
            #PARSE LEFT PRIMERS 
            
            primer3.forward.output = readLines(paste0(set.folder, '/', gene.name, '/forwardprimers.output.txt'))
            pos1 = primer3.forward.output[grep("PRIMER_LEFT_[0-9]*=", primer3.forward.output)]              

            pos2 = multi.str.split(pos1, "=", 2)
            pos3.1 = multi.str.split(pos2, ",", 1)
            pos3.2 = multi.str.split(pos2, ",", 2)
            pos4 = as.numeric(pos3.1) + (as.numeric(pos3.2)) - 1 #translate primer3 coordinate to SNP coordinate (add length to starting pos - 1)
            
            PRIMER_LEFT_X_PENALTY = as.numeric(multi.str.split(primer3.forward.output[grep("PRIMER_LEFT_[0-9]*_PENALTY=", primer3.forward.output)], "=", 2))
            PRIMER_LEFT_X_SEQUENCE = multi.str.split(primer3.forward.output[grep("PRIMER_LEFT_[0-9]*_SEQUENCE=", primer3.forward.output)], "=", 2)
            PRIMER_LEFT_X = pos2
            PRIMER_LEFT_X_TM = multi.str.split(primer3.forward.output[grep("PRIMER_LEFT_[0-9]*_TM=", primer3.forward.output)], "=", 2)
            PRIMER_LEFT_X_GC_PERCENT = multi.str.split(primer3.forward.output[grep('PRIMER_LEFT_[0-9]*_GC_PERCENT', primer3.forward.output)], "=", 2)
            PRIMER_LEFT_X_SELF_ANY_TH = multi.str.split(primer3.forward.output[grep('PRIMER_LEFT_[0-9]*_SELF_ANY_TH', primer3.forward.output)], "=", 2)
            PRIMER_LEFT_X_SELF_END_TH = multi.str.split(primer3.forward.output[grep('PRIMER_LEFT_[0-9]*_SELF_END_TH', primer3.forward.output)], "=", 2)
            PRIMER_LEFT_X_HAIRPIN_TH = multi.str.split(primer3.forward.output[grep('PRIMER_LEFT_[0-9]*_HAIRPIN_TH', primer3.forward.output)], "=", 2)
            PRIMER_LEFT_X_END_STABILITY = multi.str.split(primer3.forward.output[grep('PRIMER_LEFT_[0-9]*_END_STABILITY', primer3.forward.output)], "=", 2)
            
            
            
            left.parsed = data.frame("name1", pos4, PRIMER_LEFT_X_PENALTY, PRIMER_LEFT_X_SEQUENCE, PRIMER_LEFT_X, PRIMER_LEFT_X_TM, PRIMER_LEFT_X_GC_PERCENT, PRIMER_LEFT_X_SELF_ANY_TH, PRIMER_LEFT_X_SELF_END_TH, PRIMER_LEFT_X_HAIRPIN_TH, PRIMER_LEFT_X_END_STABILITY)
            colnames(left.parsed)[1:2] = c("p.name", "pos")
            left.parsed = left.parsed[sort(left.parsed$pos, index.return = T)$ix, ]
            left.parsed = left.parsed[which(left.parsed$pos %in% snp.coords.after.filter), ]
            nrow(left.parsed)
            left.parsed2 = split(left.parsed, left.parsed$pos)
            left.parsed2 = lapply(left.parsed2, function(x){
            x[which.min(as.numeric(x$PRIMER_LEFT_X_PENALTY)), ]
            })

            left.parsed3 = bind_rows(left.parsed2)

            left.parsed3$MSA.pos = which(homologous.snps == 1)[which(snp.coords.after.filter %in% left.parsed3$pos)]
            left.parsed3 = left.parsed3[, c(1, 2, 12, 3:11)]
            colnames(left.parsed3)[4] = "pen"
            colnames(left.parsed3)[5:ncol(left.parsed3)] = gsub("X", "0", colnames(left.parsed3)[5:ncol(left.parsed3)])

            #PARSE RIGHT PRIMERS
            
            primer3.reverse.output = readLines(paste0(set.folder, '/', gene.name, '/reverseprimers.output.txt'))
            pos1 = primer3.reverse.output[grep("PRIMER_RIGHT_[0-9]*=", primer3.reverse.output)]              

            pos2 = multi.str.split(pos1, "=", 2)
            pos3.1 = multi.str.split(pos2, ",", 1)
            pos3.2 = multi.str.split(pos2, ",", 2)
            pos4 = (as.numeric(pos3.1) - as.numeric(pos3.2)) + 1 #translate primer3 coordinate to SNP coordinate
            
            PRIMER_RIGHT_X_PENALTY = as.numeric(multi.str.split(primer3.reverse.output[grep("PRIMER_RIGHT_[0-9]*_PENALTY=", primer3.reverse.output)], "=", 2))
            PRIMER_RIGHT_X_SEQUENCE = multi.str.split(primer3.reverse.output[grep("PRIMER_RIGHT_[0-9]*_SEQUENCE=", primer3.reverse.output)], "=", 2)
            PRIMER_RIGHT_X = pos2
            PRIMER_RIGHT_X_TM = multi.str.split(primer3.reverse.output[grep("PRIMER_RIGHT_[0-9]*_TM=", primer3.reverse.output)], "=", 2)
            PRIMER_RIGHT_X_GC_PERCENT = multi.str.split(primer3.reverse.output[grep('PRIMER_RIGHT_[0-9]*_GC_PERCENT', primer3.reverse.output)], "=", 2)
            PRIMER_RIGHT_X_SELF_ANY_TH = multi.str.split(primer3.reverse.output[grep('PRIMER_RIGHT_[0-9]*_SELF_ANY_TH', primer3.reverse.output)], "=", 2)
            PRIMER_RIGHT_X_SELF_END_TH = multi.str.split(primer3.reverse.output[grep('PRIMER_RIGHT_[0-9]*_SELF_END_TH', primer3.reverse.output)], "=", 2)
            PRIMER_RIGHT_X_HAIRPIN_TH = multi.str.split(primer3.reverse.output[grep('PRIMER_RIGHT_[0-9]*_HAIRPIN_TH', primer3.reverse.output)], "=", 2)
            PRIMER_RIGHT_X_END_STABILITY = multi.str.split(primer3.reverse.output[grep('PRIMER_RIGHT_[0-9]*_END_STABILITY', primer3.reverse.output)], "=", 2)
            
            
            
            right.parsed = data.frame("name1", pos4, PRIMER_RIGHT_X_PENALTY, PRIMER_RIGHT_X_SEQUENCE, PRIMER_RIGHT_X, PRIMER_RIGHT_X_TM, PRIMER_RIGHT_X_GC_PERCENT, PRIMER_RIGHT_X_SELF_ANY_TH, PRIMER_RIGHT_X_SELF_END_TH, PRIMER_RIGHT_X_HAIRPIN_TH, PRIMER_RIGHT_X_END_STABILITY)
            colnames(right.parsed)[1:2] = c("p.name", "pos")
            right.parsed = right.parsed[sort(right.parsed$pos, index.return = T)$ix, ]
            right.parsed = right.parsed[which(right.parsed$pos %in% snp.coords.after.filter), ]
            nrow(right.parsed)
            right.parsed2 = split(right.parsed, right.parsed$pos)
            right.parsed2 = lapply(right.parsed2, function(x){
            x[which.min(as.numeric(x$PRIMER_RIGHT_X_PENALTY)), ]
            })

            right.parsed3 = bind_rows(right.parsed2)

            right.parsed3$MSA.pos = which(homologous.snps == 1)[which(snp.coords.after.filter %in% right.parsed3$pos)]
            right.parsed3 = right.parsed3[, c(1, 2, 12, 3:11)]
            colnames(right.parsed3)[4] = "pen"
            colnames(right.parsed3)[5:ncol(right.parsed3)] = gsub("X", "0", colnames(right.parsed3)[5:ncol(right.parsed3)])
            
            #######################################
            
            write.csv(left.parsed3, paste0(set.folder, '/', gene.name, '/forward.all.pen.csv'), row.names = F)
            write.csv(right.parsed3, paste0(set.folder, '/', gene.name, '/reverse.all.pen.csv'), row.names = F)

            return(list(left.parsed3, right.parsed3))
        }

        generate.best.primer.set = function(forward.all.pen, reverse.all.pen, forward.coord.used, reverse.coord.used, iteration){
            if(missing(forward.coord.used)) forward.coord.used = as.numeric()
            if(missing(reverse.coord.used)) reverse.coord.used = as.numeric()
            if(missing(iteration)) iteration = 1
            
            if(length(forward.coord.used) > 0) forward.all.pen = forward.all.pen[-which(forward.all.pen$pos %in% forward.coord.used), ]
            if(length(reverse.coord.used) > 0) reverse.all.pen = reverse.all.pen[-which(reverse.all.pen$pos %in% reverse.coord.used), ]

            best.primer.file.end.coord = as.numeric(start.coord.after.filter)
            
            forward.primer.positions = list() 
            reverse.primer.positions = list()

            minimum.snp.coord = 0
            maximum.snp.coord = start.coord.after.filter 

            #run loop to get pairs of primers

            while(best.primer.file.end.coord < end.coord.after.filter){

                
                valid.fwd.coords = which(forward.all.pen$pos < (maximum.snp.coord - 1) & forward.all.pen$pos > minimum.snp.coord)

                while(length(valid.fwd.coords) == 0){
                    maximum.snp.coord = maximum.snp.coord + 10
                    valid.fwd.coords = which(forward.all.pen$pos < (maximum.snp.coord - 1) & forward.all.pen$pos > minimum.snp.coord)
                    # print("No SNPs found for forward primer, expanding start buffer")

                    #add stop condition
                    if(maximum.snp.coord > max(coords$snp.coords)){
                    if(iteration == 1){
                        return(is.numeric())
                    } else {
                        return(is.numeric())
                    }
                    } 
                }

                f.primer.candidates = forward.all.pen[which(forward.all.pen$pos < (maximum.snp.coord - 1) & forward.all.pen$pos > minimum.snp.coord), ] 
                f.primer.candidates = f.primer.candidates[which(f.primer.candidates$pen == min(f.primer.candidates$pen)), ]

                valid.rev.coords = which(reverse.all.pen$pos > (f.primer.candidates$pos + product.size.range[1]) & reverse.all.pen$pos > best.primer.file.end.coord & reverse.all.pen$pos < (f.primer.candidates$pos + product.size.range[2]))

                while(length(valid.rev.coords) == 0){
                    product.size.range[2] = product.size.range[2] + 10
                    valid.rev.coords = which(reverse.all.pen$pos > (f.primer.candidates$pos + product.size.range[1]) & reverse.all.pen$pos > best.primer.file.end.coord & reverse.all.pen$pos < (f.primer.candidates$pos + product.size.range[2]))
                    # print("No SNPs found for reverse primer, expanding maximum product size")

                    #add stop conditions                
                    if(product.size.range[2] > length(template.sequence)){                  
                    if(iteration == 1){
                        return(is.numeric())
                    } else {                    
                        return(is.numeric())
                    }
                    } 
                }



                r.primer.candidates = reverse.all.pen[which(reverse.all.pen$pos > (f.primer.candidates$pos + product.size.range[1]) & reverse.all.pen$pos > best.primer.file.end.coord & reverse.all.pen$pos < (f.primer.candidates$pos + product.size.range[2])), ]
                r.primer.candidates = r.primer.candidates[which(r.primer.candidates$pen == min(r.primer.candidates$pen)), ]
                
                forward.primer.positions = c(forward.primer.positions, list(f.primer.candidates))
                reverse.primer.positions = c(reverse.primer.positions, list(r.primer.candidates))

                best.primer.file.start.coord = as.numeric(f.primer.candidates$pos)
                best.primer.file.end.coord = as.numeric(r.primer.candidates$pos)

                minimum.snp.coord = best.primer.file.start.coord
                maximum.snp.coord = best.primer.file.end.coord #The name of this variable originates from the first iteration. Of cause this is not truly the start.coord.after.filter on subsequent iterations
            }

            forward.primer.positions = lapply(forward.primer.positions, function(x) x[1, ])
            forward.primer.positions = bind_rows(forward.primer.positions)
            forward.primer.positions$orient = "F"


            reverse.primer.positions = lapply(reverse.primer.positions, function(x) x[1, ])
            reverse.primer.positions = bind_rows(reverse.primer.positions)
            reverse.primer.positions$orient = "R"

            list.best.primer.start.coords = c(list.best.primer.start.coords, forward.primer.positions$pos) 
            list.best.primer.end.coords = c(list.best.primer.end.coords, reverse.primer.positions$pos)
            
            return(list(forward.primer.positions, reverse.primer.positions))

        }  
        
        penalties1 = generate.all.primer.penalites(1)
        forward.all.pen = penalties1[[1]]
        reverse.all.pen = penalties1[[2]]
        
        best.primers.all = generate.best.primer.set(forward.all.pen, reverse.all.pen)
        if(length(best.primers.all) > 0){
            
            align.w.primers = matrix(data = "-", nrow = (length(mult.align1) + 1 + (nrow(best.primers.all[[1]]) * 2)), ncol = length(mult.align1[[1]]))
            align.w.primers[1:length(mult.align1), ] = as.matrix(mult.align1)

            snps.row = homologous.snps
            snps.row[which(snps.row == 0)] = "-"
            snps.row[which(snps.row == 1)] = "S"
            align.w.primers[(length(mult.align1) + 1), ] = snps.row

            primer.index = length(mult.align1) + 2
            primer.file1.multi.forward.pos.vector = as.numeric()
            primer.file1.multi.reverse.pos.vector = as.numeric()
            for(x in 1:nrow(best.primers.all[[1]])){
                
                primer.file1.seq = best.primers.all[[1]]$PRIMER_LEFT_0_SEQUENCE[x]
                
                primer.file1.seq = DNAString(primer.file1.seq)		

                
                primer.file1.multi.forward.pos = best.primers.all[[1]]$MSA.pos[x]
                primer.file1.multi.forward.pos.vector = c(primer.file1.multi.forward.pos.vector, primer.file1.multi.forward.pos)

                primer.file1.seq.rev = best.primers.all[[2]]$PRIMER_RIGHT_0_SEQUENCE[x]		

                primer.file1.seq.rev = reverseComplement(DNAString(primer.file1.seq.rev))

                primer.file1.multi.reverse.pos = best.primers.all[[2]]$MSA.pos[x]		
                primer.file1.multi.reverse.pos.vector = c(primer.file1.multi.reverse.pos.vector, primer.file1.multi.reverse.pos)

                both.seqs = c(DNAStringSet(primer.file1.seq), DNAStringSet(primer.file1.seq.rev))
                names(both.seqs) = c(paste0("forward.primer.", primer.file1.multi.forward.pos), paste0("reverse.primer.", primer.file1.multi.reverse.pos))

                #add our forward primer to the multiple alignment via a matrix		
                f.seq = strsplit(as.character(primer.file1.seq), "")
                f.seq = f.seq[[1]]
                
                align.w.primers[primer.index, (primer.file1.multi.forward.pos - (length(primer.file1.seq) - 1)):primer.file1.multi.forward.pos] = f.seq	

                #add our reverse primer to the multiple alignment via a matrix		
                r.seq = strsplit(as.character(primer.file1.seq.rev), "")
                r.seq = r.seq[[1]]
                align.w.primers[primer.index + 1, primer.file1.multi.reverse.pos:(primer.file1.multi.reverse.pos + (length(primer.file1.seq.rev) - 1))] = r.seq
                primer.index = primer.index + 2
            }

            align.row.names = c(names(mult.align1), 'SNPs', interleave(paste0('forward_', best.primers.all[[1]]$MSA.pos), paste0('reverse_', best.primers.all[[2]]$MSA.pos)))
            align.w.primers2 = conv.matrix.dnastringset(align.w.primers, align.row.names)
            
            writeXStringSet(align.w.primers2, paste0(set.folder, '/', gene.name, '/align.w.primers.fa'))
            write.csv(best.primers.all[[1]], paste0(set.folder, '/', gene.name, '/best.primers.forward.csv'), row.names = F)
            write.csv(best.primers.all[[2]], paste0(set.folder, '/', gene.name, '/best.primers.reverse.csv'), row.names = F)
        }
    }

    snp.info1 = generate_align_w_snps()    
    coords = snp.info1$coords
    find.best.primers(mult.align1, template.row, coords[[1]], coords[[2]], coords[[3]], c(min.product.size, max.product.size), F, coords = snp.info1$coords, homologous.snps = snp.info1$homologous.snps)  


}

options(show.error.locations = T)
options(stringsAsFactors = F)
options(error = function() { traceback(2); if(!interactive()) quit('no', status = 1, runLast = FALSE) })
library(dplyr)
library(parallel)

folders.sets = list.files("/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/IWGSC/all.alignments.dialign", full.names = T)
files.all = lapply(folders.sets, function(x){
    list.files(x, full.names = T)
})


library(stringr)
folders1 = sapply(files.all, function(x) gsub("set", "", str_extract(x[1], 'set[0-9]*')))
folders1 = folders1[2:length(folders1)]

files.all = files.all[-1]

Map(function(files.all, folders1){
    set.folder = paste0('set', folders1)    
    counter = 1
    mclapply(files.all, function(x){                
        tryCatch(calculate.primers(x, set.folder), error = function(e) "Error")
        counter <<- counter + 1        
        print(paste0('done ', counter))        
    }, mc.cores = 64)

    print(paste0('COMLETED SET ', folders1))
}, files.all, folders1)

