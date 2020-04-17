#written by Alex Coulton - alex.coulton@bristol.ac.uk

#writing a function to parse BLAST information
#converts a list of markers into percentage of physical distance coverage

#SETUP

# setwd("C:/Users/ac14037/project.phd.main/")
# setwd("~/project.phd.main/")

#FUNCTIONS

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