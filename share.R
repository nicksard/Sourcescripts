#################
### alf.freqs()
#################

#Description
# This fuction was originally written by Mark Christie, but has been modified to add additional
#parameters.

#Input Parameters:

# population - a table with sample.name in column one, and loci after that
# three.char - Default is F, which means the missing data is 0
#   if T then missing data is coded as something else, which is defined by missing.data
# missing.data - Default is "000", only envoked when three.char is true

alf.freq <- function(population, missing.data = "000", one.pop=F, opt = 1) {
  
  if(one.pop==T){
    
    x <- colnames(population)
    population <- data.frame("Pop",population,stringsAsFactors=F)
    colnames(population) <- c("Pops",x)
    
    pops <- unique(population[,1])
    pops <- as.character(pops)
    population <- population[,c(-2)]
    
  } else {
   
    pops <- unique(population[,1])
    population <- population[,c(-2)]
    
  }
  
  OUT1 <- NULL
  i <- 1
  i <- NULL
  for(i in 1:length(pops)){
    
    population1 <- population[population[,1] == pops[i], ]
    head(population1)
  
    #removing the sample names
    population1 <- population1[,-1]
  
  #if the data are in three character form then it replaces whatever missing data is defined as 
  if(missing.data != "000"){  
    population1[population1 == missing.data] <- "0"
    population1 <- data.frame(lapply(population1, as.numeric))
  } else {
    population1 <- data.frame(lapply(population1, as.numeric))
  }
  
  #this bit grabs the locus names for the population1
  L <- 1:ncol(population1)
  locus.locs <- L[seq(1,ncol(population1),2)]
  locus.names <- colnames(population1)
  
  OUT <- NULL 
  
  #this for loop grabs the two columns of alleles for each locus
  #and uses table to count the number of times each allele is observed
  for (x in locus.locs) {
    
    alleles <- c(population1[,x],population1[,x+1])
    alleles2 <- as.data.frame(table(alleles))
    missing.alleles <- alleles2[which(alleles2[,1]==0),2]
    
    #removes the missing values  
    if(length(which(alleles2[,1]==0))>0) {
      alleles2=alleles2[-which(alleles2[,1]==0),]
    }
    
    #this calculates the frequency and puts it all in a table
    alleles4 <- cbind(alleles2,alleles2[,2]/sum(alleles2[,2])) 
    output <- cbind(pops[i],locus.names[x],alleles4)                      
    OUT <- as.data.frame(rbind(OUT,output),stringsAsFactors=F)
  }

    OUT1 <- as.data.frame(rbind(OUT1,OUT),stringsAsFactors=F)
    OUT1  
  }
  
  #adds the column names
  colnames(OUT1) <- c("Pops","Locus","allele","count","frequency")
  Allelefreqs <- OUT1
  
  if(opt == 1){ 
    Allelefreqs$Pops <- as.character(Allelefreqs$Pops)
    Allelefreqs$Locus <- as.character(Allelefreqs$Locus)
    Allelefreqs$allele <- as.character(Allelefreqs$allele)
    return(Allelefreqs) }

  if(opt == 2){
    
    alfs <- Allelefreqs[,c("Locus","frequency")]
    alfs$Locus <- as.character(alfs$Locus)
    locs <- unique(alfs$Locus)
    
    OUT <- NULL
    for(i in 1:length(locs)){
      num <- nrow(alfs[alfs$Locus == locs[i],])
      OUT1 <- rep(as.numeric(paste(i)),times=num)
      OUT <- c(OUT,OUT1)
    }
    alfs$Locus <- OUT
    return(alfs)
  }

}


#################
### gt.summary()
#################

#Description
# This fuction was originally written by Mark Christie, but has been modified to add additional
# lines to calculated observed and expected heterozygosity. In addition it calculates the
# percent genotyped

#Input Parameters:

# population - a table with sample.name in column one, and loci after that
# three.char - Default is F, which means the missing data is 0
#   if T then missing data is coded as something else, which is defined by missing.data
# missing.data - Default is "000", only envoked when three.char is true
# pop.sum - Default is F, which means it will not summarize pop gene stats accross loci
# one.pop - Default is T, which means that it will automatically assume you only have one population
# round.num - Default set to 4, meaning it will round the population summary to 4 decimals



gt.summary <- function(population, missing.data = "000", pop.sum=F, one.pop=T, round.num=4) {
  
  if(one.pop==T){
    
    x <- colnames(population)
    population <- data.frame("Pop",population,stringsAsFactors=F)
    colnames(population) <- c("Pops",x)
    
    pops <- unique(population[,1])
    pops <- as.character(pops)
    population <- population[,c(-2)]
    
  } else{
    
    pops <- unique(population[,1])
    pops <- as.character(pops)
    population <- population[,c(-2)]
    
  }
  

  i <- NULL
  OUT1 <- NULL
  for(i in 1:length(pops)){
    
    population1 <- population[population[,1] == pops[i], ]

    #removing the sample names
    population1 <- population1[,-1]
    
    #if the data are in three character form then it replaces whatever missing data is defined as 

    population1[population1 == missing.data] <- "0"
    population1 <- data.frame(lapply(population1, as.numeric))
    
    #this bit grabs the locus names for the population1
    L <- 1:ncol(population1)
    locus.locs <- L[seq(1,ncol(population1),2)]
    locus.names <- colnames(population1)
    
    OUT <- NULL 

    #counting the number of individuals in the sample
    x <- NULL  
    for (x in locus.locs) {
      tot <- nrow(population1)
      
      #calculating percent genotyped
      missing <- nrow(population1[population1[,x]==0,])
      missing1 <- nrow(population1[population1[,x]==0,])
      Percent.gt <- round((((tot-missing)/tot)*100),2)
      
      #calculating expected heterozygosity
      alleles <- c(population1[,x],population1[,x+1])
      num.alleles <- length(unique(alleles))
      alleles2 <- as.data.frame(table(alleles))
      missing <- alleles2[which(alleles2[,1]==0),2]
      
      #removing zeros
      if(length(which(alleles2[,1]==0))>0) {
        alleles2 <- alleles2[-which(alleles2[,1]==0),]
      }
      
      #calculating p for all loci and then squaring them for all alleles, and calculating He using Nei 1987 eq. 8.4
      N <- tot- as.numeric(missing1)
      alleles4 <- data.frame(alleles2,alleles2[,2]/sum(alleles2[,2])) 
      alleles4$p2 <- alleles4[,3]^2
      bias.he <- 1-sum(alleles4$p2)
      bias.he
      correction <- (2*N)/(2*N-1)
      correction
      exp.het <- round(correction*bias.he,2)
      exp.het
      
      #calculating observed heterozygosity
      alleles <- data.frame(cbind(population1[,x],population1[,x+1]))
      alleles$truth <- ifelse(alleles$X1 == alleles$X2,T,F)
      missing <- nrow(subset(alleles, X1 == 0 & X2 == 0))
      tot2 <- nrow(alleles)-missing
      het <- nrow(subset(alleles, truth ==  F))
      obs.het <- round(het/tot2,2)
      fis <- round(1- (obs.het/exp.het),2)
      
      output <- data.frame(pops[i],locus.names[x], num.alleles, exp.het, obs.het, fis, tot, missing, Percent.gt)
      OUT <- as.data.frame(rbind(OUT,output),stringsAsFactors=F)
    }
    
    OUT1 <- rbind(OUT1,OUT)
    
  }
  
  colnames(OUT1) <- c("Pop","Locus","A","He","Ho","Fis","N","N.Missing","P.GT")
  OUT1 <- as.data.frame(OUT1,stringsAsFactors=F)  
  OUT1$Pop <- as.character(OUT1$Pop)
  OUT1$Locus <- as.character(OUT1$Locus)
  
  if(pop.sum==F){
    return(OUT1)
  } else {
    
    #getting unique pop identifiers
    pops <- as.character(unique(OUT1[,1]))
    
    #
    df <- OUT1
    i <- 1
    i <- NULL
    OUT <- NULL
    
    for(i in 1:length(pops)){
      #getting just the pop of interest
      df1 <- df[df[,1] == pops[i],]
      df1
      
      #getting mean and standard error for number of alleles
      avg.a <- round(mean(df1[,"A"]),round.num)
      se.a <- round(st.err(df1[,"A"]),round.num)
      
      #getting mean and standard error for expected heterozygosity
      avg.he <- round(mean(df1[,"He"]),round.num)
      se.he <- round(st.err(df1[,"He"]),round.num)  
      
      #getting mean and standard error for observed heterozygosity
      avg.ho <- round(mean(df1[,"Ho"]),round.num)
      se.ho <- round(st.err(df1[,"Ho"]),round.num)  
      
      #getting mean and standard error for Fis
      avg.fis <- round(mean(df1[,"Fis"]),round.num)
      se.fis <- round(st.err(df1[,"Fis"]),round.num)  
      
      #getting mean and standard error for n.missing
      avg.nm <- round(mean(df1[,"N.Missing"]),round.num)
      se.nm <- round(st.err(df1[,"N.Missing"]),round.num)  
      
      #getting mean and standard error for percent genotyped
      avg.pgt <- round(mean(df1[,"P.GT"]),round.num)
      se.pgt <- round(st.err(df1[,"P.GT"]),round.num)  
      
      #putting it all together and making a nice neat DF
      OUT1 <- cbind(df1[1,"N"], avg.a, avg.he, avg.ho, avg.fis, avg.nm, avg.pgt, se.a, se.he, se.ho, se.fis, se.nm, se.pgt)
      OUT1 <- as.data.frame(OUT1, stringsAsFactors=F)
      colnames(OUT1) <- c("N","A","He","Ho","Fis","N.Missing","P.GT","A.se","He.se","Ho.se","Fis.se","N.Missing.se","P.GT.se")
      OUT1$Pop <- pops[i]
      OUT1 <- move.me(OUT1,"Pop","first")
      OUT <- rbind(OUT,OUT1)
    }
    return(OUT)
  }
}


#################
### six.char()
#################

#Description
# This fuction takes genotypes in three character form and puts them into six character form

#Input Parameters:

# population - a table with sample.name in column one, and loci after that. All loci must be in three character form. Examples 000, 111, 999
six.char <- function(population, missing.data = "0", one.pop=T) {
  
  missing.data <- as.numeric(missing.data)
  missing.data <- paste0(missing.data,missing.data)
  missing.data1 <- paste0(missing.data,missing.data,missing.data) 
  cnames <- colnames(population)
  
  if(one.pop==T){  
    pop.names <- population[,1]
    population=population[,-1]  
  } else{
    pop.names <- population[,1]
    id.names <- population[,2]
    population=population[,-c(1,2)]
  }

  L<-1:ncol(population)
  locus.locs <- L[seq(1,ncol(population),2)]
  locus.names <- colnames(population)
  
  head(population)
  OUT <- NULL
  for (x in locus.locs) {
    output=paste(population[,x],population[,x+1],sep="")
    OUT <- cbind(OUT,output)
  }
  
  locus.names <- locus.names[seq(1,length(locus.names),2)]
  
  head(OUT)
  if(one.pop==T){
    OUT <- data.frame(pop.names,OUT)
    colnames(OUT) <- c(cnames[1],locus.names)
  } else {
    OUT <- data.frame(pop.names,id.names,OUT)
    head(OUT)
    colnames(OUT) <- c(cnames[1:2],locus.names)
  }

  OUT <- data.frame(lapply(OUT, as.character), stringsAsFactors=FALSE)
  OUT[OUT == missing.data] <- missing.data1
  head(OUT)
  return(OUT)
}

###################
### dup.identifer()
###################

#Description
# This fuction identifies any genotypes that are the same and reports which ones are

#Input Parameters:

# tmp - a table with sample.name in column one, and loci after that
# If the remove.dups = T, then the first genotype in the file is kept and all others that
# have the same genotype are removed.


dup.identifer <- function(tmp, one.pop = T){
  
  if(one.pop == T){
    #this pastes all the genotypes for each locus together  
    tmp$gt<- apply(tmp[,-1],1,paste,collapse="")  
  } else {
    #this pastes all the genotypes for each locus together  
    tmp$gt<- apply(tmp[,-c(1,2)],1,paste,collapse="")      
  }

  
  #the duplicated function is used to duplicated() is used to identify any duplicated genotypes
  tmp$logic <- duplicated(tmp$gt)
  
  if(one.pop == T){
    #this grabs the duplicated ids
    dup.ids <- tmp[tmp$logic == T, 1]
    dup.ids  
  } else {
    #this grabs the duplicated ids
    dup.ids <- tmp[tmp$logic == T, 2]
    dup.ids      
  }

  head(tmp)
  if(length(dup.ids) == 0){
    
    #cleaning up df
    tmp$gt <- NULL
    tmp$logic <- NULL
    print("There were no duplicated genotypes present in data")
    return(tmp)
    
  } else {
   
    #grabbing the genotypes that were genotypes
    dup.gts <- tmp[tmp$logic == T,]$gt
    dup.gts <- unique(dup.gts)
    ids <- paste0(rep("dup.group."),1:length(dup.gts))
    
    tmp$ids <- "unique"
    for(i in 1:length(dup.gts)){
      
      #now changing the logic to T for all genotypes that were duplicated, which is how to find all
      #genotypes that were duplicated
      tmp$ids[tmp$gt ==  dup.gts[i]]<- ids[i]
    }
    
    #cleaning up df
    tmp$gt <- NULL
    tmp$logic <- NULL
    return(tmp)
  } 
}


#################
### three.char()
#################

#Description
# This fuction takes genotypes in six character form and puts them into three character form
#Input Parameters:

# population - a table with sample.name in column one, and loci after that. All loci must be in six character form. Examples 000000, 111111, 999999

three.char <- function(tmp, one.pop=T){
  
  if(one.pop==T){
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1])
    colnames(OUT) <- colnames(tmp)[1]
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1]  
  } else{
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1:2])
    colnames(OUT) <- colnames(tmp)[1:2]
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1:-2]
  }

  
  #counting the number of loci and saving their names
  L <- ncol(gts)
  Locus.names <- colnames(gts)
  
  #using a small for loop to take each column and split it into two columns and cbinding that to OUT, which originally just had the sample names in it
  for(i in 1:ncol(gts)){
    locus <- data.frame(gts[,i])
    colnames(locus) <- Locus.names[i]
    
    locus$x1 <- substr(locus[,1],1,3)
    locus$x2 <- substr(locus[,1],4,6)
    locus[,1] <- NULL
    
    colnames(locus) <- c(paste(Locus.names[i],".1", sep=""),paste(Locus.names[i],".2",sep=""))
    OUT <- cbind(OUT,locus)
  }
  OUT[,1] <- as.character(OUT[,1])
  return(OUT)
}

####################
### m.data.counter()
####################

#Description
# This fuction counts how many loci are missing data for each sample
# It requires data that are in three character form

# Input parameters
# 1) tmp - table with genotypes in three character form and sample names in column one
# 2) opt - options that can be set to 1 or 2
#             Option 1 - returns a summary of the number of missing loci in dataset
#             Option 2 - returns entire dataset with a new column that counts the number of loci missing data
# 3) missing.data - default "000", but if its different such as "999", this can be changed


m.data.counter <- function(tmp, opt = 1, missing.data = "000", one.pop=T){
  
  #saving original genotypes
  tmp1 <- tmp
  
  if(one.pop==T){
    
    #replacing missing data with 0 and making entire data frame numeric
    tmp[tmp == missing.data] <- "0"
    tmp <- data.frame(tmp[,1],lapply(tmp[,-1], as.numeric))
    colnames(tmp) <- c("Sample.Names",colnames(tmp)[-1])
    
    #saving the names for each sample
    ids <- setNames(data.frame(tmp[,1],stringsAsFactors=F),"Sample.Name")
    
    #removing the names from the original data frame
    tmp <- tmp[,-1]
  
  } else {
    
    #replacing missing data with 0 and making entire data frame numeric
    tmp[tmp == missing.data] <- "0"
    tmp <- data.frame(tmp[,c(1,2)],lapply(tmp[,-c(1,2)], as.numeric))
    colnames(tmp) <- colnames(tmp1)
    
    #saving the names for each sample
    ids <- data.frame(tmp[,c(1,2)],stringsAsFactors=F)
    colnames(ids) <- colnames(tmp1)[1:2]
    
    #removing the names from the original data frame
    tmp <- tmp[,-c(1,2)]
  }
  
  #identifying the missing data
  tmp[tmp == 0] <- -10000 
  tmp[tmp != -10000] <- 0 
  tmp[tmp == -10000] <- 1
  
  #calculating the number of missing loci
  tmp$sum <- rowSums(tmp[,c(1:ncol(tmp))])
  tmp$sum <- tmp$sum/2
  
  tmp1$m.data <- tmp$sum
  tmp1[,1] <- as.character(tmp1[,1])
    
  return(tmp1)
}


#################
### uniq.alleles
#################

#Description
# This fuction takes output from alf.freqs() and The program compares
# the alleles at each locus between the two programs and outputs only the alleles
# that are unique to each population along with their counts and frequencies.

#Input Parameters:

# alfs - requires output from alf.freqs() and the data.frame used in that function
# had only two populations


uniq.alleles <- function(alfs) {
  
  x1 <- NULL
  y1 <- NULL
  
  loci <- as.character(unique(alfs$Locus))
  
  for(i in 1:length(loci)){
    alfs1 <- alfs[alfs$Locus == loci[i],]
    
    pops <- as.character(unique(alfs1$Pops))
    
    pop1 <- alfs1[alfs1$Pops == pops[1],]
    pop1 <- as.numeric(as.character(pop1$allele))
    
    pop2 <- alfs1[alfs1$Pops == pops[2],]
    pop2 <- as.numeric(as.character(pop2$allele))
    
    x <- data.frame(pop1,match(pop1,pop2))
    colnames(x) <- c("alleles","matches")
    x <- x[is.na(x$matches) == T,1]
    x <- alfs[alfs$Pops == pops[1] & alfs$Locus == loci[i] & alfs$allele %in% x, ]
    x1 <- rbind(x1,x)
    
    y <- data.frame(pop2,match(pop2,pop1))
    colnames(y) <- c("alleles","matches")
    y <- y[is.na(y$matches) == T,1]
    y <- alfs[alfs$Pops == pops[2] & alfs$Locus == loci[i] & alfs$allele %in% y, ]
    y1 <- rbind(y1,y)
  }
  
  OUT <- rbind(x1,y1)
  return(OUT)
}

#################
### st.err()
#################

#Description
# This fuction calculates standard error

st.err <- function(x){ sd(x)/sqrt(length(x))}

#################
### move.me()
#################

#Description
#This fuction reorders columns in your dataframe in anyway you want. It was orignally developed by Ananda Mahto and found it
#on stackoverflow at this link: http://stackoverflow.com/questions/18339370/reordering-columns-in-a-large-dataframe. By default
#it takes the column(s) you want and moves them to the end of the file, but if designated can move them to the front or infront of
#any specified column

#Input Parameters:
#data - a dataframe
#tomove - can be single or multiple character string
#ba - tells the program what columns to stick tomove columns in front of

move.me <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}

#################
### sim.pop.gts()
#################

#Description
# This fuction creates possion distributed random data to represent microsatellite alleles. 
#  It will repeat that process for however many individuals, loci, and populations you wish.

# Input parameters
# 1) n - the number of individuals per population
# 2) loci - the number of loci desired
# 3) missing.data - can put some percentage of missing data randomly distributed thoroughout the data
# 4) pops - the number of populations desired
# 5) k - mean number of alleles per loci
# 6) k.sd - standard devation of number of alleles per loci


#making the function
sim.pop.gts <- function(n=1000,loci=10,missing.data=0, pops=2, k=20, k.sd=3){
  
  #Pops
  Pops <- rep("POP",times=n*pops)
  pop.ids <- rep(1:pops,each=n)
  Pops <- paste(Pops, pop.ids, sep= ".")
  
  output <- NULL
  for(j in 1:pops){
    
    #IDs
    nums <-1:n
    IDs <- rep("ID.",times=n)
    IDs <- paste0(IDs,nums)
    

    i <- NULL
    OUT <- NULL
    for(i in 1:loci){
      
      #getting the alleles for each locus
      x <- rpois(n=100000, lambda=115)
      x <- sort(x)
      
      x <- x[x>100 & x < 200]
      xs <- unique(x)
      y <- sample(xs,size=round(rnorm(n=1,mean=k,sd=k.sd)),replace=F)
      x <- x[x %in% y]
      
      
      #sampling those alleles the number of times required to
      #make complete genotypes (allele 1 and allele 2)
      z.1 <- sample(x,size=n,replace=T)
      z.2 <- sample(x,size=n,replace=T)
      
      #combining to make a complete genotype in three character form
      OUT1 <- cbind(z.1,z.2)
      
      #making sure the smaller allele is first
      OUT1 <- as.data.frame(OUT1, stringsAsFactors=F)
      OUT1$logic <- OUT1[,1] > OUT1[,2]
      OUT1.1 <- OUT1[OUT1$logic == F,-ncol(OUT1)]
      OUT1.2 <- OUT1[OUT1$logic != F,-ncol(OUT1)]
      OUT1.2 <- OUT1.2[,c(2,1)]
      colnames(OUT1.2) <- colnames(OUT1.1)
      OUT1 <- rbind(OUT1.1,OUT1.2)
      
      if(missing.data != 0){
        #making missing data
        num.md <- missing.data*n
        md.locs <- sample(1:n,size=num.md,replace=F)
        OUT1[md.locs,1] <- "000" 
        OUT1[md.locs,2] <- "000"
      }
      OUT1 <- as.matrix(OUT1)
      OUT <- cbind(OUT,OUT1)
    } # end of for loci loop 
    
    df <- data.frame(IDs,OUT,stringsAsFactors=F)
    loc.names <- rep("Locus",times=2*loci)
    loc.nums <- rep(1:loci, each=2)
    loc.alleles <- rep(1:2,times=loci)
    loc.names <- paste(loc.names,loc.nums,sep="_")
    loc.names <- paste(loc.names,loc.alleles,sep=".")
    colnames(df) <- c("Sample.name",loc.names)
    
    output <- rbind(output,df)
  } #end of Pops loop
  
  #getting rid of any factors!
  output <- data.frame(Pops,output,stringsAsFactors=F)
  return(output)
} #end of sim.pop.gts


####################
### sim.afreqs() ###
####################

#Description
#This fuction is a modified version of Mark Christie's script to create simulate genotype allele frequencies.
#Note that there are high and low frequencies for each locus, and the alleles and their associated frequencies
#are exactly the same. 

#Input Parameters:
#nlcoi - define the number of loci you want
#n.alleles - define the number of alleles you want


sim.afreqs <- function(nloci=10,n.alleles=10, simple=F){
  
  #making sure nloci and n.alleles are numeric
  nloci <- as.numeric(nloci)
  n.alleles <- as.numeric(n.alleles)
  
  #getting the alleles frequencies for each allele
  a <- c(1:n.alleles)
  b <- rev(a)
  c <- b/a
  d <- sum(c)
  
  #now making actual allele frequencies
  freqs <- c/d
  
  #setting the lowest allele to an arbitrary number
  lowest.allele <- 140
  
  #making afreqs
  
  #getting loci names
  loci <- paste0("Locus.",1:nloci)
  loci <- rep(loci, each=n.alleles)
  
  #for alternative output
  loci.1 <- rep(1:nloci, each=n.alleles)
  
  #getting the allele names
  alleles2 <- cbind(seq(lowest.allele,lowest.allele+n.alleles-1,1),freqs)
  alleles <- rep(alleles2[,1],nloci)
  
  
  #getting their frequencies
  freqs  <-  rep(alleles2[,2],nloci)
  
  #returning finished product
  if(simple == F){
    
    afreqs <- data.frame(loci,alleles,freqs, stringsAsFactors=F)
    return(afreqs)  
    
  } else {
    
    afreqs <- data.frame(loci.1,freqs, stringsAsFactors=F)
    return(afreqs)  
  }
  
  
} # end of sim.afreqs

#####################
### sim.gt.data() ###
#####################

#Description
#This fuction is a modified version of Mark Christie's script to create parents and offspring.
#Note that one can create random genotyping error, add unrelated indviduals, and add full siblings inteh mix

#Input Parameters:
#afreq - this is a data frame with two rows in it, the first is the locus, and the second is the allele frequencies
#Nparents - defining thenumber of parents one wants
#Noffs_perpair - defining the number of offspring per parent one wants
#error - defining the amount of random error in the dataset
#Nunrelated - defining the number of unrelated individuals wanted in the dataset
#Nsibs - defining the number of full siblings one wants

#Create Simulated Data Sets==========================================#

sim.gt.data <-  function(afreq, Nparents= 2, Noffs_perpair= 1, error= 0, Nunrelated= 0, Nsibs=0, write2file = F){
  
  #making sure everything is numeric
  Nparents <- as.numeric(Nparents)
  Noffs_perpair <- as.numeric(Noffs_perpair)
  error <- as.numeric(error)
  Nunrelated <- as.numeric(Nunrelated)
  
  #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
  Nadults <- Nparents*2             
  
  #getting afreqs input
  afreqs <- afreq

  #making a function to simulate genotypes for parents
  
  OUT <- NULL
  sims <- function(sims){
    
    #table allele frequencies
    #getting just the frequencies 
    alleles2 <- afreqs[afreqs[,1] == loci[i],]
    
    #changing generic name to three character genotypes
    alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                  
    
    #create homozygote allele frequencies
    homos <- (alleles3[,2])^2                                                          
    homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
    
    #create heterozygote allele frequencies
    
    #first create all possible combinations of freqs
    hets <- t(combn(alleles3[,2],2))                                                   
   
     #now make expected freqs
    hetfreq <- 2*(hets[,1]*hets[,2])

    #create heterozygote allele names
    #now getting the het. genotype combinations
    hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
    
    #combing them
    hets2 <- cbind(hetvals,hetfreq)
    
    #combine hets and homos and create genotypes
    gfreqs <- rbind(hets2,homos2)                                                      
    
    #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
    n <- 1000000                                    
    
    #create genotypes(by coloumn, 1 for each allele)
    #for the first allele
    
    gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3])))) 
    
    #now the second
    gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))

    
    #combining them
    gtypes <- cbind(gfreqs1,gfreqs2)
    
    #mixing them up
    gtypes <- gtypes[sample(1:length(gtypes[,1]),replace= F),]
    
    #now getting a random sample of the gentoypes for the parents
    sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
    
    OUT<<-cbind(OUT,sg1)
    
  } # end of sim function
  
  #defining the number of loci there are
  loci.length <- length(unique(afreqs[,1]))
  loci <- unique(afreqs[,1])
  
  #doing the simulation process for each locus
  for(i in 1:length(loci)){
    lapply(1,sims)
  } # end for loop

  #saving OUT into object called parents
  parents <- OUT
  
  #next making a matrix to add to parents so that gts in every two columns don't overlap
  #first making a matrix the same size of parents and filling it ever two columns with varying
  #numbers in the 1000s
  c <- c(1:(ncol(OUT)))
  odd <- 2*(unique(round(((c-2))/2)))+1
  l <- length(odd) * 1000
  codes <- seq(from= 1,to= l,by= 1000)
  cols <- sort(rep(codes,2))-1
  Anumbs <- matrix(data= cols,nrow= Nadults,ncol= ncol(OUT),byrow= T)

  #now adding it to parents
  parents <- as.numeric(parents)+Anumbs

  #create full sib families (go down the list in pair) ======================#
  OUT2 <- NULL
  sims <- function(sims){
    
    #first grabbing the genotypes of the parents
    p1 <- parents[i,]
    p2 <- parents[i+1,]
    #defining the order of teh alleles
    als <- rep(1:2,length(p1)/2)
    
    #defining the number of offspring that need to be full sibs
    Noffs <- Noffs_perpair                  
    
    #using a for loop to making full sib genotypes
    OUT2 <- NULL
    for (b in 1:Noffs){
      
      #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
      #sampling alleles from parent one
      pos1 <- sample(als,length(p1)/2,replace <- TRUE)                                      
      
      #sampling alleles from parent two
      pos2 <- sample(als,length(p1)/2,replace <- TRUE)

      #getting the position for those alleles from parent one
      pos11 <- pos1+(seq(0,(length(p1)-1),2))
      
      #getting the position for those alleles from parent two
      pos22 <- pos2+(seq(0,(length(p2)-1),2))
      
      #getting those alleles from parent one
      o1 <- p1[pos11]
      #getting those alleles from parent two
      o2 <- p2[pos22]
      
      #putting them together to make the offsprings genotype
      o3 <- sort(c(o1,o2))
      
      #identifying what parents the offspring came fromm
      o3 <- t(c(i,o3))
      
      #writting to file and appending
      write.table(o3,file <- "SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
    } #end of full sib for loop
  } # end of new sims function
  
  #defining which parents are going to have full sibs
  my.picks <- seq(from=1,to=Nadults-1, by=2)
  
  #making full sib gts for each parent pair
  for(i in my.picks){ 
    lapply(i,sims)
  } # end of full sib genotyping loop
  
  #removing the 1000s from the gts
  Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)                                   
  parents <- as.numeric(parents)-Anumbs
  
  #getting parent genotypes
  #first the moms
  Moms <- parents[seq(from=2,to= Nadults,by= 2),]
  #now the dads
  Dads <- parents[seq(from=1,to= Nadults,by= 2),]
  
  #reading SimOffs.txt back in and removing the id column
  Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Offs  <-  as.matrix(Offs2[,-1])
  
  #removing issue with 1000s
  Anumbs <- matrix(cols,length(Offs[,1]),ncol(OUT),byrow = T)
  Offs <- as.numeric(Offs)-Anumbs
  
  #now making offspring names
  if (Noffs_perpair>1){
    
    Offnames <- ceiling(Offs2[,1]/2)                
    Offnames2 <- paste0(Offnames,".",1:length(Offnames))
    Offs <- cbind(paste("Offspring",Offnames2),Offs)
    
  } else {
    
    Offnames <- ceiling(Offs2[,1]/2)
    Offs <- cbind(paste("Offspring",Offnames),Offs)
  } # end of if statement about offspring
  
  Moms <- cbind(paste("Mom",1:length(Moms[,1])),Moms)
  Dads <- cbind(paste("Dad",1:length(Dads[,1])),Dads)
  
  #add error======================================================================#
  error
  #first for the dads
  
  #getting the gts
  Dadsg <- Dads[,-1]
  
  #figuring out how many gts to make errors in
  ldad <- length(Dadsg)*error
  
  #randomly selecting where to put the errors
  pdad1 <- sample(1:length(Dadsg),ldad,replace= FALSE)
  
  #randomly selecting some gets to replace in the those locations
  pdad2 <- Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
  
  #actually putting in the error
  Dadsg2 <- replace(Dadsg,pdad1,pdad2)
  
  #replacing the old Dads gts
  Dads <- cbind(Dads[,1],Dadsg2)
  
  #doing the same thing for the kids
  Offsg=Offs[,-1]
  loff=length(Offsg)*error
  poff1=sample(1:length(Offsg),loff,replace=FALSE)
  poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
  Offsg2=replace(Offsg,poff1,poff2)
  Offs=cbind(Offs[,1],Offsg2)
  
  #no error rate for moms - why?
  #Momsg <- Moms[,-1]
  #ldad <- length(Momsg)*error
  #pdad1 <- sample(1:length(Momsg),ldad,replace=FALSE)
  #pdad2 <- Momsg[sample(1:length(Momsg),ldad,replace=FALSE)]
  #Momsg2 <- replace(Momsg,pdad1,pdad2)
  #Moms <- cbind(Moms[,1],Momsg2)
  
  #===============================================================================#
  #making generic locus names
  colsnmz <- rep("Locus",ncol(Offs)-1)
  colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
  colsnmz3 <- paste(colsnmz,colsnmz2)
  colsnmz3 <- c("IDs",colsnmz3)
  
  #usig them to make the colnames for the Offs, Dads, and Moms
  colnames(Offs)<- colsnmz3
  colnames(Dads)<- colsnmz3
  colnames(Moms)<- colsnmz3
  
  #Add unrealated individuals <- ====================================================#
  
  if (Nunrelated>0){
    
    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nunrelated*3                                                            
    afreqs <- afreq
    OUT <- NULL
    
    #making the allele frequencies again, this way they are completely random and unrelated
    #doing the same thing as lines 34-113
    sims <- function(sims) {
      alleles2 <- afreqs[which(afreqs[,1] == z),]
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])
      
      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2                                                         
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
      
      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))                                                   
      hetfreq <- 2*(hets[,1]*hets[,2])
      
      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
      hets2 <- cbind(hetvals,hetfreq)
      
      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)                                                      
      
      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000         
      
      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
      
      OUT <<-cbind(OUT,sg1)
      
    } #end of sims function
    
    #now applying function
    z1 <- length(unique(afreqs[,1]))
    for(z in 1:z1) {lapply(z,sims)}
    
    #repeating lines 122-204
    parents <- OUT
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    parents <- as.numeric(parents)-Anumbs
    
    #to called this thing unrelated
    unrelated <- cbind(paste("Individual",1:length(parents[,1])),parents)
    Lid <- length(unrelated[,1])/3
    
    #splitting the unrelated individuals among dads, moms, and kids
    m1 <- unrelated[1:Lid,]
    d1 <- unrelated[(Lid+1):(Lid*2),]
    j1 <- unrelated[((Lid*2)+1):(Lid*3),]
    
    #now for the column names again
    colsnmz <- rep("Locus",ncol(Offs)-1)
    colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
    colsnmz3 <- paste(colsnmz,colsnmz2)
    colsnmz3 <- c("IDs",colsnmz3)
    colnames(Offs)<- colsnmz3
    colnames(Dads)<- colsnmz3
    colnames(Moms)<- colsnmz3
    
    #adding the unrelated individuals in the three files
    Moms <- rbind(Moms,m1)
    Dads <- rbind(Dads,d1)
    Offs <- rbind(Offs,j1)
    
  } #end of unrelated if statement
  
  #removing the file SimOff.txt from working directory
  unlink("SimOffs.txt")
  
  #Begin add siblings=============================================================#
  #This scripts creates a desired number of siblings and splits them between adults and offspring files
  Nsibs <- as.numeric(Nsibs)
  if (Nsibs > 0) {
    #could change this, but as stands only evaluating 1 pair of siblings per parent-pair
    Noffs_per_pair <- 2                                                                
    
    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nsibs*2                                                                 
    afreqs <- afreq
    
    #repeating allele frequency simulation from above
    OUT <- NULL
    sims <- function(sims){
      alleles2 <- afreqs[which(afreqs[,1]==z),]
      
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   
      
      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2                                                          
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
      
      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))                                                   
      hetfreq <- 2*(hets[,1]*hets[,2])
      
      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
      hets2 <- cbind(hetvals,hetfreq)
      
      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)                                                      
      
      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000                                                                       
      
      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
      
      OUT <<-cbind(OUT,sg1)
      
    } #end of allele frequency simulation
    
    #applyingn to making genotypes of parengs
    z <- length(unique(afreqs[,1]))
    for(z in 1:z) {lapply(z,sims)}
    parents <- OUT
    
    #repeating the part where he addes 1000 to each line
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    
    
    #create full sib families (go down the list in pair)============================#
    OUT2 <- NULL
    sims <- function(sims)  {
      N <- 1:length(parents[,1])
      u <- sample(N,1)
      u2 <- sample(N,1)
      p1 <- parents[u,]
      p2 <- parents[u2,]
      als <- rep(1:2,length(p1)/2)
      
      #number of offspring per pair
      Noffs <- Noffs_per_pair                                                            
      OUT2 <- NULL
      for (b in 1:Noffs){
        
        #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
        pos1 <- sample(als,length(p1)/2,replace=TRUE)                                      
        pos2 <- sample(als,length(p1)/2,replace=TRUE)
        pos11 <- pos1+(seq(0,(length(p1)-1),2))
        pos22 <- pos2+(seq(0,(length(p2)-1),2))
        o1 <- p1[pos11]
        o2 <- p2[pos22]
        o3 <- sort(c(o1,o2))
        o3 <- t(c(z,o3))
        write.table(o3,file="SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
      } #end of for loop
    } # end of simulation for sibs
    
    #applying it to the parents
    z <- length(parents[,1])
    
    #used to move down list, now sample randomly so is just used to produce wanted number of offspring
    C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,sims)                 
    Dads <- parents[seq(from=1,to=length(parents[,1]),by=2),]
    Moms <- parents[seq(from=2,to=length(parents[,1]),by=2),]
    
    #see code before functions for adding 1000s (here am removing 1000s)
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)                                   
    parents <- as.numeric(parents)-Anumbs
    Dads <- parents[seq(from=1,to=length(parents[,1]),by=2),]
    Moms <- parents[seq(from=2,to=length(parents[,1]),by=2),]
    
    #reading in the kids again
    Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    Offs  <-  as.matrix(Offs2[,-1])
    Anumbs <- matrix(cols,length(Offs[,1]),ncol(OUT),byrow <- T)
    Offs <- as.numeric(Offs)-Anumbs
    
    #naming offspring
    if (Noffs_per_pair>1){
      #naming of offpsring
      Offnames <- ceiling(Offs2[,1]/2)                                                   
      Offnames2 <- paste(Offnames,".",1:length(Offnames))
      Offs <- cbind(paste("Sibling",Offnames2),Offs)
    } else {
      Offnames <- ceiling(Offs2[,1]/2)
      Offs <- cbind(paste("Sibling",Offnames),Offs)
    } # end of offspring naming if statement
    
    Moms <- cbind(paste("Mom",1:length(Moms[,1])),Moms)
    Dads <- cbind(paste("Dad",1:length(Dads[,1])),Dads)
    
    
    #add error======================================================================#
    Offsg <- Offs[,-1]
    loff <- length(Offsg)*error
    poff1 <- sample(1:length(Offsg),loff,replace=FALSE)
    poff2 <- Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
    Offsg2 <- replace(Offsg,poff1,poff2)
    Offs <- cbind(Offs[,1],Offsg2)
    
    #===============================================================================#
    colsnmz <- rep("Locus",ncol(Offs)-1)
    colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
    colsnmz3 <- paste(colsnmz,colsnmz2)
    colsnmz3 <- c("IDs",colsnmz3)
    colnames(Offs)<- colsnmz3
    
    #calculate shared alleles among pairs of siblings===============================#
    sib1 <- Offs[seq(from=1,to=length(Offs[,1]),by=2),]
    sib2 <- Offs[seq(from=2,to=length(Offs[,1]),by=2),]
    write.table(sib1,file="Sib1.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(sib2,file="Sib2.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    
    #appending siblings to file 
    #getting dad siblings
    dadsibs <- read.table("Dads_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    sib1 <- read.table("Sib1.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    dadsibs <- rbind(dadsibs,sib1)
    #write.table(dadsibs,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)
    
    #putting them in the kid file
    juvsibs <- read.table("Juveniles_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    sib2 <- read.table("Sib2.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    juvsibs <- rbind(juvsibs,sib2)
    #write.table(juvsibs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)
    
    #adding them to the end of the parent files
    Dads <- rbind(Dads,dadsibs)
    Offs <- rbind(Offs,juvsibs)
    
    #removing files
    unlink("Sib1.txt")
    unlink("Sib2.txt")
    unlink("SimOffs.txt")
  } # end of sibling function
  
  
  if(write2file == T){
    print("Mom, Dad, and Offspring genotypes have been saved to current working directory")  
    write.table(Dads,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(Moms,file="Moms_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(Offs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
  } else {
    
    Parents <- rbind(Moms,Dads)
    Parents <- data.frame("Parents",Parents,stringsAsFactors=F)
    colnames(Parents) <- c("Type", colnames(Parents)[-1])
    
    Offspring <- data.frame("Offspring",Offs,stringsAsFactors=F)
    colnames(Offspring) <- c("Type",colnames(Offspring)[-1])
    
    output <- data.frame(rbind(Parents,Offspring),stringsAsFactors=F)
    #making the genotypes numeric
    output <- data.frame(output[,1:2],sapply(output[,3:ncol(output)], as.numeric),stringsAsFactors=F)
    return(output)
  }
} # end of sim.gt.data function


########################
### power.analysis() ###
########################

#Description
#This fuction is a modified version of Mark Christie's script to calculate the exclusion power within a dataset
# It will calculates the expected number of false parent offspring pairs (POPs) within the dataset and
# it also calculates the exclusion power per locus


#Input Parameters:
#adults - this is a data frame with parents ids and genotypes
#offspring - this is a data frame with offspring ids and genotypes
#loci_number - define as the number of loci you want to analyze
#             NOTE: from what I understand it randomly selects which ones to analyze
#opt - has three versions
#       opt 1 - writes both Exp. Number of false POPs and Locus.Power to file
#       opt 2 - returns Exp. Number of false POPs
#       opt 3 - reutnrs Locus.Power

power.analysis <- function(adults,offspring, loci_number=0,opt=1){
  
  if(length(loci_number)==1){
    if(loci_number == 0){print("Need to define the number of loci to analyze before I can continue..."); stop}
    multi <- T
  } else {
    multi <- F
  }
  
  if(multi==F){
    
    #getting the adults and offspring gts without their sample names
    Adults1 <- adults
    a <- ncol(Adults1)
    lociset <- seq(from= 2,to= a,by= 2)
    loci <- as.numeric(loci_number)
    s1 <- sample(lociset,loci,replace= FALSE)
    s2 <- sort(c(s1,s1+1))
    Adults <- Adults1[,s2]
    Offspring1 <- offspring
    Offspring <- Offspring1[,s2]
    
    head(Offspring)
    head(Adults)
    a <- 1
    #Begin calcualtion of Pr(Z) and Pr(delta)=======================================#
    massive <- function(massive){
      
      #making a table with unique alleles, their counts, and frequencies for parents
      AAT40 <- c(Adults[,i],Adults[,(i+1)])
      locus <- as.data.frame(table(AAT40))
      Frequency <- locus[,2]/sum(locus[,2])
      Data <- cbind(locus, Frequency)
      Data
      
      #removing missing data from the table
      if (Data[1,1]==0) {
        
        Data <- Data[-1,]
        
      } else {
        
        Data <- Data
      } #end of missing data remover
      
      
      #calculating the expected number of each homozygote allele
      n1 <- ((length(AAT40))/2)
      A1 <- (Data[,3]*(2*n1))
      A2 <- (Data[,3]^2)*n1
      AA <- A1-A2
      AAA <- cbind(Data,AA)
      AAA
      
      #making a table with unique alleles, their counts, and frequencies for offspring
      AAT402 <- c(Offspring[,i],Offspring[,(i+1)])
      locus2 <- as.data.frame(table(AAT402))
      Frequency2 <- locus2[,2]/sum(locus2[,2])
      Data2 <- cbind(locus2, Frequency2)
      Data2
      
      #removing missing data again
      if (Data2[1,1]==0){
        
        Data2=Data2[-1,]
        
      } else {
        
        Data2=Data2 
      } #end of missing data remover
      
      #calculating the expected number of each homozygote allele
      n2 <- ((length(AAT402))/2)
      B1 <- (Data2[,3]*(2*n2))
      B2 <- (Data2[,3]^2)*n2
      BB <- B1-B2
      BBB <- cbind(Data2,BB)
      BBB
      
      #getting the alleles in the parents that are found in the kids
      AAA1 <- match(BBB[,1],AAA[,1])
      g <- which(AAA1>0)
      g1 <- AAA1[g]
      AAA1 <- AAA[g1,]
      AAA1
      
      #doing the same thing for the kids
      BBB1 <- match(AAA[,1],BBB[,1])
      j <- which(BBB1>0)
      j1 <- BBB1[j]
      BBB1 <- BBB[j1,]
      
      #multiplying the number of expected from parents and kids and summing 
      AB <- AAA1[,4]*BBB1[,4]
      AB <- sum(AB)
      AB
      
      #getting all possible combinations of alleles and removing homozygotes for the parents
      y <- AAA1[,1]
      Agenotypes <- expand.grid(y,y)
      remove <- which(Agenotypes[,1]==Agenotypes[,2])
      Agenotypes <- Agenotypes[-remove,]
      Agenotypes
      
      #doing the same for the kids
      z <- BBB1[,1]
      Bgenotypes <- expand.grid(z,z)
      remove <- which(Bgenotypes[,1]==Bgenotypes[,2])
      Bgenotypes <- Bgenotypes[-remove,]
      Bgenotypes
      
      #getting the expected number of parents that would have each heterozygous genotype
      one <- match(Agenotypes[,1],AAA1[,1])
      two <- match(Agenotypes[,2],AAA1[,1])
      one <- AAA1[one,3]
      two <- AAA1[two,3]
      Agfreq <- one*two*n1*2
      Agfreq
      
      #getting the expected number of kids that would have each heterozygous genotype
      oneb <- match(Bgenotypes[,1],BBB1[,1])
      twob <- match(Bgenotypes[,2],BBB1[,1])
      oneb <- BBB1[oneb,3]
      twob <- BBB1[twob,3]
      Ogfreq <- oneb*twob*n2*2
      Ogfreq
      
      #combining the two and calculating the probability of being false at that locus
      Gfreqs <- Agfreq*Ogfreq
      Gfreqs <- floor(Gfreqs)
      Gfreq <- sum(Gfreqs)/2
      PrB <- (AB-Gfreq)/(n1*n2)
      write.table(PrB,file="Output_Per_Locus_Exclusion_Probabilities.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
      
    } # end of massive function
    
    a <- ncol(Adults)
    loci.locs <- seq(from=1,to=ncol(Adults)-1, by=2)
    loci.locs
    
    #applying massive function to all loci
    for(i in loci.locs) {
      lapply(i,massive)
    } # end of for loop calculating exclusion probs
    
    #reading the file with each locus' exclusion prob in
    PrBs <- read.table("Output_Per_Locus_Exclusion_Probabilities.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    
    #getting the number of parents and offspring
    nn1 <- length(Adults[,1])
    nn2 <- length(Offspring[,1])
    
    #getting the product of all the independent probs
    Pr_delta <- prod(PrBs[,1])
    
    #getting exclusion power per locus
    Locus.Power <- data.frame(colnames(Offspring)[loci.locs],PrBs,stringsAsFactors=F)
    colnames(Locus.Power) <- c("Loci","Exc.Power")
    Locus.Power
    
    #making table with information in it
    Expected_Number_of_False_Pairs <- Pr_delta*(nn1*nn2)
    pvalue1 <- data.frame(loci_number,Expected_Number_of_False_Pairs,Pr_delta,stringsAsFactors=F)
    
    #removing the saved file
    unlink("Output_Per_Locus_Exclusion_Probabilities.txt")
    
  } else {
    
    pvalue1 <- NULL
    for(i in 1:length(loci_number)){
      
      # getting the adults and offspring gts without their sample names
      Adults1 <- adults
      a <- ncol(Adults1)
      lociset <- seq(from= 2,to= a,by= 2)
      loci <- as.numeric(loci_number[i])
      s1 <- sample(lociset,loci,replace= FALSE)
      s2 <- sort(c(s1,s1+1))
      Adults <- Adults1[,s2]
      Offspring1 <- offspring
      Offspring <- Offspring1[,s2]
      
      head(Offspring)
      head(Adults)
      a <- 1
      #Begin calcualtion of Pr(Z) and Pr(delta)=======================================#
      massive <- function(massive){
        
        #making a table with unique alleles, their counts, and frequencies for parents
        AAT40 <- c(Adults[,i],Adults[,(i+1)])
        locus <- as.data.frame(table(AAT40))
        Frequency <- locus[,2]/sum(locus[,2])
        Data <- cbind(locus, Frequency)
        Data
        
        #removing missing data from the table
        if (Data[1,1]==0) {
          
          Data <- Data[-1,]
          
        } else {
          
          Data <- Data
        } #end of missing data remover
        
        
        #calculating the expected number of each homozygote allele
        n1 <- ((length(AAT40))/2)
        A1 <- (Data[,3]*(2*n1))
        A2 <- (Data[,3]^2)*n1
        AA <- A1-A2
        AAA <- cbind(Data,AA)
        AAA
        
        #making a table with unique alleles, their counts, and frequencies for offspring
        AAT402 <- c(Offspring[,i],Offspring[,(i+1)])
        locus2 <- as.data.frame(table(AAT402))
        Frequency2 <- locus2[,2]/sum(locus2[,2])
        Data2 <- cbind(locus2, Frequency2)
        Data2
        
        #removing missing data again
        if (Data2[1,1]==0){
          
          Data2=Data2[-1,]
          
        } else {
          
          Data2=Data2 
        } #end of missing data remover
        
        #calculating the expected number of each homozygote allele
        n2 <- ((length(AAT402))/2)
        B1 <- (Data2[,3]*(2*n2))
        B2 <- (Data2[,3]^2)*n2
        BB <- B1-B2
        BBB <- cbind(Data2,BB)
        BBB
        
        #getting the alleles in the parents that are found in the kids
        AAA1 <- match(BBB[,1],AAA[,1])
        g <- which(AAA1>0)
        g1 <- AAA1[g]
        AAA1 <- AAA[g1,]
        AAA1
        
        #doing the same thing for the kids
        BBB1 <- match(AAA[,1],BBB[,1])
        j <- which(BBB1>0)
        j1 <- BBB1[j]
        BBB1 <- BBB[j1,]
        
        #multiplying the number of expected from parents and kids and summing 
        AB <- AAA1[,4]*BBB1[,4]
        AB <- sum(AB)
        AB
        
        #getting all possible combinations of alleles and removing homozygotes for the parents
        y <- AAA1[,1]
        Agenotypes <- expand.grid(y,y)
        remove <- which(Agenotypes[,1]==Agenotypes[,2])
        Agenotypes <- Agenotypes[-remove,]
        Agenotypes
        
        #doing the same for the kids
        z <- BBB1[,1]
        Bgenotypes <- expand.grid(z,z)
        remove <- which(Bgenotypes[,1]==Bgenotypes[,2])
        Bgenotypes <- Bgenotypes[-remove,]
        Bgenotypes
        
        #getting the expected number of parents that would have each heterozygous genotype
        one <- match(Agenotypes[,1],AAA1[,1])
        two <- match(Agenotypes[,2],AAA1[,1])
        one <- AAA1[one,3]
        two <- AAA1[two,3]
        Agfreq <- one*two*n1*2
        Agfreq
        
        #getting the expected number of kids that would have each heterozygous genotype
        oneb <- match(Bgenotypes[,1],BBB1[,1])
        twob <- match(Bgenotypes[,2],BBB1[,1])
        oneb <- BBB1[oneb,3]
        twob <- BBB1[twob,3]
        Ogfreq <- oneb*twob*n2*2
        Ogfreq
        
        #combining the two and calculating the probability of being false at that locus
        Gfreqs <- Agfreq*Ogfreq
        Gfreqs <- floor(Gfreqs)
        Gfreq <- sum(Gfreqs)/2
        PrB <- (AB-Gfreq)/(n1*n2)
        write.table(PrB,file="Output_Per_Locus_Exclusion_Probabilities.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
        
      } # end of massive function
      
      a <- ncol(Adults)
      loci.locs <- seq(from=1,to=ncol(Adults)-1, by=2)
      loci.locs
      
      #applying massive function to all loci
      for(i in loci.locs) {
        lapply(i,massive)
      } # end of for loop calculating exclusion probs
      
      #reading the file with each locus' exclusion prob in
      PrBs <- read.table("Output_Per_Locus_Exclusion_Probabilities.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
      
      #getting the number of parents and offspring
      nn1 <- length(Adults[,1])
      nn2 <- length(Offspring[,1])
      
      #getting the product of all the independent probs
      Pr_delta <- prod(PrBs[,1])
      
      #getting exclusion power per locus
      Locus.Power <- data.frame(colnames(Offspring)[loci.locs],PrBs,stringsAsFactors=F)
      colnames(Locus.Power) <- c("Loci","Exc.Power")
      Locus.Power
      
      #making table with information in it
      Expected_Number_of_False_Pairs <- Pr_delta*(nn1*nn2)
      pvalue1.1 <- data.frame(Expected_Number_of_False_Pairs,Pr_delta,stringsAsFactors=F)
      pvalue1 <- rbind(pvalue1, pvalue1.1)
      
      #removing the saved file
      unlink("Output_Per_Locus_Exclusion_Probabilities.txt")
    }
    pvalue1 <- cbind(loci_number,pvalue1)
  }
  
  #writing both to file
  if(opt == 1){
    
    write.table(Locus.Power,file="Locus power.txt",row.names=FALSE,col.names=T,sep="\t",append=T)
    write.table(pvalue1,file="Expected Number of False Pairs.txt",row.names=FALSE,col.names=T,sep="\t",append=T)
    
    print("Locus power.txt saved to current working directory")
    print("Expected Number of False Pairs.txt saved to current working directory")
  } #end of option 1
  
  #returning only the expected number of false pairs
  if(opt == 2){
    
    return(pvalue1)
  } #end of option 2
  
  #returning the exc. power per locus
  if(opt == 3){
    
    return(Locus.Power)
  } #end of option 3
  
} #end of power analysis function

#################
### NEP.calc()
#################

#Description
#This funciton calculations non-exclusion power for any given set of loci with any given set of alleles
#originally from Jamieson, A. 1997

#Input Parameters:
#tmp - a data.frame with only two columns, the first is the name of the loci and the second the frequencies
#      of the alleles at that locus - in long form 

NEP.calc <- function(tmp, opt=1){
  loci <- unique(tmp[,1])
  
  i <- NULL
  OUT <- NULL
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- -2*sum(tmp1^2)
    p1
    
    p2 <- sum(tmp1^3)
    p2
    
    p3 <- 2*sum(tmp1^4)
    p3
    
    p4 <- -3*sum(tmp1^5)
    p4
    
    p5 <- -2*(sum(tmp1^2)^2)
    p5
    
    p6 <- 3*sum(tmp1^2)*sum(tmp1^3)
    p6
    
    P <- 1+(sum(p1,p2,p3,p4,p5,p6))
    
    OUT <- c(OUT,P)
  }
  
  OUT
  
  i<-1
  i <- NULL
  OUT1 <- NULL
  
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    p1 <- -4*sum(tmp1^2)
    p1
    
    p2 <- 2*(sum(tmp1^2)^2)
    p2
    
    p3 <- 4*sum(tmp1^3)
    p3
    
    p4 <- -3*sum(tmp1^4)
    p4
    
    
    P <- 1+(sum(p1,p2,p3,p4))
    
    OUT1 <- c(OUT1,P)
  }
  
  OUT1
  
  
  
  i<-1
  i <- NULL
  OUT2 <- NULL
  
  for(i in 1:length(loci)){
    tmp1 <- tmp[tmp[,1] == loci[i],2]
    
    p1 <- 4*sum(tmp1^4)
    p1
    
    p2 <- -4*sum(tmp1^5)
    p2
    
    p3 <- -3*sum(tmp1^6)
    p3
    
    p4 <- -8*(sum(tmp1^2)^2)
    p4
    
    p5 <- 8*sum(tmp1^2)*sum(tmp1^3)
    p5
    
    p6 <- 2*(sum(tmp1^3)^2)
    p6
    
    P <- 1+(sum(p1,p2,p3,p4,p5,p6))
    
    OUT2 <- c(OUT2,P)
  }
  
  
  #parparing for first for loops
  i <-  NULL
  PId <- NULL
  
  for(i in 1:length(loci)){
    
    #grabbing each locus individually and ordering some smallest to largest
    tmp1 <- tmp[tmp[,1] == loci[i],]
    tmp1 <- tmp1[order(tmp1[,2]),]
    tmp1
    
    #raising the allele frequencies to the fourth power
    p14 <-  sum(tmp1[,2]^4)
    p14 #makes test in excel
    
    #parparing for second for loops
    j <- NULL
    calc <- NULL
    
    for(j in 1:nrow(tmp1)){
      
      #take each allele frequencies one at a time
      tmp1.f <- tmp1[j,2]
      tmp1.f
      
      #and taking all allele frequencies greater than that and (2*Pi*Pj)^2 and summing
      tmp2 <- tmp1[-1:-j,2]    
      calc1 <- sum((2*tmp1.f*tmp2)^2)
      
      #adding to calc
      calc <- c(calc,calc1)
    }
    
    #summing all calcs
    calc <- sum(calc)
    calc
    
    #adding p14 and calc
    PId1 <- p14 +calc
    
    #adding to PId
    PId <- c(PId,PId1)
  }
  
  
  
  #parparing for first for loops
  i <-  NULL
  PrSib <- NULL
  
  for(i in 1:length(loci)){
    
    #grabbing each locus individually and ordering some smallest to largest
    tmp1 <- tmp[tmp[,1] == loci[i],]
    tmp1 <- tmp1[order(tmp1[,2]),]
    tmp1
    
    tmp1$p2 <- tmp1[,2]^2
    tmp1$p4 <- tmp1[,2]^4
    
    tmp1
    
    0.5*sum(tmp1$p2)
    (0.5*sum(tmp1$p2))^2
    -0.25*sum(tmp1$p4)
    
    p1 <- 0.25
    p2 <- 0.5*sum(tmp1$p2)
    p3 <- 0.5*(sum(tmp1$p2)^2)
    p4 <- -0.25*(sum(tmp1$p4))
    c(p1,p2,p3,p4)
    ps1 <- sum(p1,p2,p3,p4)
    PrSib <- c(PrSib,ps1)
  }
  
  
  #getting the power calcs accross all loci and adding it to the individual locus calcs
  Locus <- c(loci,"Combined")
  Locus
  
  NE_2P <- c(OUT , prod(1-OUT))
  NE_2P
  
  NE_1P <- c(OUT1 , prod(1-OUT1))
  NE_1P
  
  NE_PP <- c(OUT2 , prod(1-OUT2))
  NE_PP
  
  PrId <- c(PId,prod(PId))
  PrId
  
  PrIdSib <- c(PrSib,prod(PrSib))
  PrIdSib
  
  Output <- data.frame(Locus,NE_1P,NE_2P,NE_PP, PrId,PrIdSib, stringsAsFactors=F)
  
  if(opt==1){
    return(Output)
  }
  
  if(opt==2){
    Output <- Output[Output$Locus == "Combined",]
    return(Output)
  }
}


#################
### exclusion()
#################

#Description
#This funciton takes a data.frame with parents and a data.frame with offspring and identify all parent offspring
#relationships with up to the number of mismatches you describe and it has two alternative options for output
#take from SOLOMON

#Input Parameters:
#adults - contains the sample names and genotypes of the parents
#offspring - contains the sample names and genotypes of the offspring
#mismatch - the maximum number of loci allowed to mistmatch between parent and offspring
#opt - has two options for output
#   1) parents and offspring with genotyeps, and number of mismatching loci
#   2) parents, offspring, and number of mismatching loci

exclusion <-function(adults, offspring, mismatch=0, opt= 1){
  
  # making sure MM is numeric and saving adults and offspring to other objects
  mismatch <- as.numeric(mismatch)                                                 
  Adults1 <- adults
  Offspring1 <- offspring
  
  #getting just the number of columns so gts and ids can be seperated
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]
  
  #getting the names for adults and offspring
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]
  
  #getting the number of genotypes, and how many adults and juveniles there are 
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])
  
  #making all possible comparisons between adults and offspring using generic numbers
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  
  #setting up two dfs with the gts of interest in the appropriate order
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
    
  #getting the IDs
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]

  #writing names to file
  write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  
  #reading them back in and combining them
  IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  
  #making a function to identify all possible matches between parent and offspring pairs (POPs) at a given locus
  matches <- function(matches){
    #comparing the genotypes of adult to offspring
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]
    
    #multipling them together and getting ride of NAs
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    #replacing those greater than 0 with 1, this is to ID those that do not match
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)
    
    #writing this to file and appending, I see the trick to avoid for loops
    write.table(f,file="Sort.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=T)
  }
  
  #getting the columns to apply this to
  z <- ncol(Ads)
  
  #applying the match function to all loci
  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
  
  #reading results back in
  Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  
  #converting from long form to wide form
  a <- unique(Observed[,1])
  i<- NULL
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  } # end of wide form convert
  
  #adding up the riws for each column
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  
  #identifying the matches
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)
  
  #identifying the locations of the names where matches occurred
  IDS <- which(stuff<(mismatch+1))
  
  #getting those ids for offspring and parents
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]
  
  #getting the number of putative pairs and filling a DF with them
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  Putativepairs
  
  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"
  
  #getting the gts for each parent and offspring drawing from f, which appeared earlier
  sorts <- function(sorts) {
    tell <- rbind(PAdults[f,],POffspring[f,])
    write.table(tell,file="Output_genotypes_tmp.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
  } #end of sort loop
  
  f <- length(PAdults[,1])
  
  #applying to all parents
  for(f in 1:f) lapply(f,sorts)
  
  #removing sort IdnamesA and IdnaemsO
  unlink("Sort.txt")
  unlink("IdnamesA.txt")
  unlink("IdnamesO.txt")
  
  #Outputfile processing==========================================================#
  putative<- read.table("Output_genotypes_tmp.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
  unlink("Output_genotypes_tmp.txt")  
  if(opt == 1){
    
    return(putative)
    
  } #end option 1
  
  if(opt == 2){
    
    #seperating parents and offspring again
    parents <- putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
    offspring <- putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]

    #now making two columns with IDs for parents and offspring and their number of MMs
    out1 <- cbind(parents,offspring)
    out2 <- out1[order(out1[,1]),-2]
    colnames(out2) <- c("Parent","Offspring","MM")

    
    return(out2)
  } #end option 2
  
} #end of exclusion


#################
### exclusion.ns()
#################

#Description
#This funciton takes a data.frame with parents and a data.frame with offspring and identify all parent offspring
#relationships with up to the number of mismatches you describe and it has two alternative options for output
#take from SOLOMON

#Input Parameters:
#adults - contains the sample names and genotypes of the parents
#offspring - contains the sample names and genotypes of the offspring
#mismatch - the maximum number of loci allowed to mistmatch between parent and offspring
#opt - has two options for output
#   1) parents and offspring with genotyeps, and number of mismatching loci
#   2) parents, offspring, and number of mismatching loci

exclusion.ns <-function(adults, offspring, mismatch=0, opt= 1){
  
  # making sure MM is numeric and saving adults and offspring to other objects
  mismatch <- as.numeric(mismatch)                                                 
  Adults1 <- adults
  Offspring1 <- offspring
  
  #getting just the number of columns so gts and ids can be seperated
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]
  
  #getting the names for adults and offspring
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]
  
  #getting the number of genotypes, and how many adults and juveniles there are 
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])
  
  #making all possible comparisons between adults and offspring using generic numbers
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  G
  
  #setting up two dfs with the gts of interest in the appropriate order
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
  head(Ads)
  head(Offs)
  
  #getting the IDs
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]
  IdnamesA
  IdnamesO
  
  getwd()
  #writing names to file
  write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  
  #reading them back in and combining them
  IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  Names
  
  z <- 1
  z <- NULL
  #making a function to identify all possible matches between parent and offspring pairs (POPs) at a given locus
  matches <- function(matches){
    #comparing the genotypes of adult to offspring
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]
    
    #multipling them together and getting ride of NAs
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    #replacing those greater than 0 with 1, this is to ID those that do not match
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)
    
    #writing this to file and appending, I see the trick to avoid for loops
    write.table(f,file="Sort.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=T)
  }
  
  #getting the columns to apply this to
  z <- ncol(Ads)
  
  #applying the match function to all loci
  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
  
  #reading results back in
  Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Observed
  
  #converting from long form to wide form
  a <- unique(Observed[,1])
  i<- NULL
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  } # end of wide form convert
  
  #adding up the riws for each column
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  head(Sorted)
  
  #identifying the matches
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)
  
  #identifying the locations of the names where matches occurred
  IDS <- which(stuff<(mismatch+1))
  
  #getting those ids for offspring and parents
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]
  
  #getting the number of putative pairs and filling a DF with them
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  Putativepairs
  
  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  PAdults
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  POffspring
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"
  
  PAdults$sorter <- as.numeric(row.names(PAdults))
  POffspring$sorter <- as.numeric(row.names(POffspring))
  putative <- rbind(PAdults,POffspring)
  putative <- putative[order(putative$sorter),]
  putative$sorter <- NULL
  
  #getting the gts for each parent and offspring drawing from f, which appeared earlier
  #sorts <- function(sorts) {
  #tell <- rbind(PAdults[f,],POffspring[f,])
  #  write.table(tell,file="Output_genotypes_tmp.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
  #} #end of sort loop
  
  #f <- length(PAdults[,1])
  
  #applying to all parents
  #for(f in 1:f) lapply(f,sorts)
  
  #removing sort IdnamesA and IdnaemsO
  unlink("Sort.txt")
  unlink("IdnamesA.txt")
  unlink("IdnamesO.txt")
  
  #Outputfile processing==========================================================#
  #putative<- read.table("Output_genotypes_tmp.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
  #putative  
  #unlink("Output_genotypes_tmp.txt")  
  if(opt == 1){
    
    return(putative)
    
  } #end option 1
  
  if(opt == 2){
    
    #seperating parents and offspring again
    parents <- putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
    offspring <- putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]
    parents
    offspring
    
    #now making two columns with IDs for parents and offspring and their number of MMs
    out1 <- cbind(parents,offspring)
    out2 <- out1[order(out1[,1]),-2]
    colnames(out2) <- c("Parent","Offspring","MM")
    head(out2)
    
    return(out2)
  } #end option 2
  
} #end of exclusion.ns

##################
### no.par.bayes()
##################

#Description
#This funciton takes a data.frame with parents and a data.frame with offspring and identify all parent offspring
#relationships with up to the number of mismatches bayes on bayesian priors and gives posterior probabilities for each
#take from SOLOMON - Mark Christie's library

#Input Parameters:
#adults - contains the sample names and genotypes of the parents
#offspring - contains the sample names and genotypes of the offspring
#reps - set the number of simulated datasets desired
#Ngts - set the number of simulated genotypes desired

no.par.bayes <- function(adults,offspring,reps=1000,Ngts=50000000)     {
  
  #loading in adults and offspring genotypes
  Adults1 <- adults
  Offspring1 <- offspring
  
  #making sure the number of reps is numeric and counting the number of loci
  reps <- as.numeric(reps)
  loci <- ncol(Adults1)
  
  #removing ID column
  Adults <- Adults1[,c(2:loci)]                                                    
  Offspring <- Offspring1[,c(2:loci)]
  
  #getting the number loci, again for Progress bar
  total <- ncol(Adults)/2                                                       
  
  #putting in a progress bar
  
  
  #Begin Master simulation function
  afreqs <- function(afreqs){
    
    #making locus names? *** need to figure out ** looks like used in progress bar
    locus_name <- L
    
    #setting up progress bar again
    
    
    #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
    vect <- c(Adults[,L],Adults[,L+1])                                                 
    alleles <- data.frame(table(vect))
    alleles <- alleles[order(alleles[,1]),]
    if (as.numeric(as.character(alleles[1,1]))==0) {alleles <- alleles[-1,]}           #removes missing data
    if(length(alleles[,1])==1) {                                                    #deals with monomorphic loci by adding 1 very strange allele (799)
      alleles <- cbind(vect[1],alleles[2])
      alleles[2,]<-c(799,1)}
    alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
    homos <- (alleles2[,3])^2                                                          #create homozygote allele frequencies
    homos2 <- cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
    hets <- t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
    hetfreq <- 2*(hets[,1]*hets[,2])
    hetvals <- t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
    hets2 <- cbind(hetvals,hetfreq)
    gfreqs <- rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
    csum <- cumsum(as.numeric(gfreqs[,3]))
    gfreqs1 <- cbind(gfreqs,csum)
    Nadults <- length(Adults[,1])
    Noffs <- length(Offspring[,1])
    
    
    #===============================================================================#end locus-specific HWE genotype frequency calculations
    alength <- length(alleles2[,1])
    for (Y in 1:reps) {
      positions <- 1:length(gfreqs[,1])
      sg3 <- sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                      #sample the repeated genotype positions, by the number of adults
      sadults <- gfreqs[sg3,1:2]                                                         #index gfreqs to create genotypes
      og3 <- sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                        #create juvenile genotyes
      soffs <- gfreqs[og3,1:2]
      soffs <- cbind(as.numeric(soffs[,1]),as.numeric(soffs[,2]))
      asp <- cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
      asp <- rbind(cbind(locus_name,0,0),asp)
      for (i in 1:Nadults) {
        parent1 <- as.numeric(sadults[i,1])                                                #first allele in parent
        parent2 <- as.numeric(sadults[i,2])                                                #second allele in parent
        p1 <- soffs-parent1
        p2 <- soffs-parent2
        pp1 <- which(p1[,1]==0)
        pp2 <- which(p1[,2]==0)
        allele1 <- unique(c(pp1,pp2))
        p21 <- which(p2[,1]==0)
        p22 <- which(p2[,2]==0)
        allele2 <- unique(c(p21,p22))
        Out51 <- cbind(parent1,length(allele1))
        Out52 <- cbind(parent2,length(allele2))
        Out53 <- cbind(0,Noffs-length(unique(c(allele1,allele2))))
        Out5 <- rbind(Out51,Out52,Out53)
        if(parent2==parent1) {Out5 <- Out5[-1,]}                                           #remove 1 of alleles count if homozygous
        if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
          diffs <- sum(Out5[,2])-Noffs                                                     #comment out to be more conservative!
          maxa <- max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
          pos <- which(Out5[,2]==maxa)
          Out5[pos,2]<-Out5[pos,2]-diffs}
        m1 <- match(Out5[,1],asp[,2])
        m2 <- asp[m1,3]+as.numeric(Out5[,2])
        asp[m1,3]<-m2
        asp<-asp
      }
      write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  L <- ncol(Adults)
  C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
  
  #Bayes P-value==================================================================#
  OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  locname <- unique(OUT[,1])                                                         #compile calculations for each locus
  OUTALL <- NULL
  for (z in locname) {
    Loc1 <- OUT[which(OUT[,1]==z),]
    allfreqs <- unique(Loc1[,2])
    OUT2 <- NULL
    for (x in allfreqs) {
      a1<-Loc1[which(Loc1[,2]==x),]
      a2 <- sum(a1[,3])
      a3 <- cbind(x,a2)
      OUT2 <- rbind(OUT2, a3)
    }
    OUT3 <- cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
    OUTALL <- rbind(OUTALL, cbind(z,OUT3))
  }
  #Create multilocus genotypes, 50000 at a time, calculate number of loci that mismatch, and calculate freqs of shared alleles
  NL <- length(unique(OUTALL[,1]))
  ngtypes <- 1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
  Ngts <- as.numeric(Ngts)
  inreps <- 50000     #was tested as 10000 for SNPS
  repnumber <- round(Ngts/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  asp <- cbind(0:NL,rep(0,length(0:NL)))
  
  for (n in 1:repnumber) {
    
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    distm <- apply(OUT, 1, function(x)sum(x == 1))
    distm2 <- data.frame(table(distm))
    m1 <- match(distm2[,1],asp[,1])
    m2 <- asp[m1,2]+distm2[,2]
    asp[m1,2]<-m2
    asp<-asp
  }
  
  #tabulate number of multilocus genotypes with x mismatches
  d2 <- asp
  d3 <- cbind(d2,d2[,2]/sum(d2[,2]))
  #===============================================================================# Create plot of exclusionary power.  Also, some necessary data formatting
  
  
  Adults<- adults
  Offs<- offspring
  Nads <- length(Adults[,1])
  Noffs <- length(Offs[,1])
  asp <- d3
  asp <- cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
  asp <- cbind(asp,cumsum(asp[,2]))
  asp <- cbind(asp,asp[,5]/max(asp[,5]))
  #find minimum Nloci to mismatch (could modify this) and perform exclusion with decided-upon mismatches==#
  distm <- cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
  #mismatch=min(which(round(distm[,6],1)==.9))                                    #deprecated
  mismatch <- which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
  Adults1 <- Adults                                                                   #begin exclusion
  Offspring1 <- Offs
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]
  write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  matches <- function(matches)
  {
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)
    write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  }
  z <- ncol(Ads)
  C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
  Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)
  IDS <- which(stuff<(mismatch+1))
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"
  sorts <- function(sorts)
  {
    tell <- rbind(PAdults[f,],POffspring[f,])
    write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  }
  f <- length(PAdults[,1])
  C1 <-  for(f in 1:f) lapply(f,sorts)
  unlink("Sort.txt")
  unlink("IdnamesA.txt")
  unlink("IdnamesO.txt")
  #Calculate phi for each number mismatching loci=================================#
  
  Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  observed <- data.frame(table(Putative[,2]))
  observed <- cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
  zerom <- 0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
  zerom2 <- which(is.na(match(zerom,observed[,1])))
  if (length(zerom2>0)) {observed <- observed[,3]
                         for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}  #not really 0, to prevent divide by 0
  }   else {observed=observed[,3]}
  expected <- distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
  #expected=distm[1:(mismatch+1),3]*Nads*Noffs                                     #not using cumulative sum
  phi <- expected/observed
  phi <- replace(phi,which(phi>=1),1)
  Offs<- offspring
  actualTrue <- length(grep("Off",Offs[,1]))
  info <- cbind(actualTrue,expected,observed,phi)
  #calculate phi and index values ================================================#
  
  phibase <- phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
  observed <- observed[min(which(phi[]==1))-1]                                       #do the same for observed
  testob <- which(observed==0.000000001)
  phi2 <- cbind(1:length(phi),phi)
  if (length(testob>0)) {phi4 <- phi2[-testob,]} else {phi4 <- phi2}
  nmismatch <- min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
  index <- phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
  index <- index[which(index>-1)]
  if (length(index)>1) {
    if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
  phi <- phi[index]
  index <- index-1
  #Create Plot ===================================================================#
  pdf(file <- "Output_Dataset_Power.pdf")
  x <- 0:(length(info[,1])-1)
  y1 <- info[,3]
  y2 <- info[,2]
  p1 <- which(y2==0)
  y2 <- replace(y2,p1,.000000001)
  p1 <- which(y1<y2)
  y3 <- replace(y1,p1,y2[p1])
  y2 <- y2+1
  y3 <- y3+1
  par(mar = c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
  ats <- c(0,1,2,3,4,5,6,7,8,9)
  labs <- c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
  axis(side=2,ats,labs)
  lines(x,log10(y3),lwd=2)
  lines(x,log10(y2),lwd=2)
  points(x,log10(y3),pch=21,bg="green",cex=2)
  points(x,log10(y2),pch=21,bg="blue",cex=2)
  legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
  yphi <- y2/y3
  par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
  lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
  points(x,yphi,cex=2,pch=21,bg="gray")
  mtext("Number of Mismatching Loci",side=1,line=1.94)
  dev.off()
  info2 <- cbind(x,info[,4])
  colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
  write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  #===============================================================================#
  ngtypes <- 20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
  inreps <- 10000
  repnumber <- round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  #writes values at all loci, to be analyzed further below
  for (n in 1:repnumber) {
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      alleles3 <- alleles3[-which(alleles3[,2]==0),]
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    for (i in index) {                                                              #loop over numbers of mismatched loci
      if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                    #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distp <- cbind(i,DIST)
      write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  
  #Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
  OUT9 <- NULL
  for (n in index){
    Putative2 <- Putative[which(Putative[,2]==n),]
    write.table(Putative2, file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                                #need at least 1000 values , else assigned 0 (may be redundant now)
    OBS <- NULL
    afreqs <- function(afreqs) {
      PutL <- Putative2[,c(L+2,L+3)]
      PutLadults <- PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
      PutLoffs <- PutL[seq(from=2,to=length(PutL[,1]),by=2),]
      Puts <- cbind(PutLadults,PutLoffs)
      c1 <- Puts[,1]-Puts[,3]                                                    #find the matching alleles
      c2 <- Puts[,1]-Puts[,4]
      c3 <- Puts[,2]-Puts[,3]
      c4 <- Puts[,2]-Puts[,4]
      Puts2 <- cbind(Puts,c1,c2,c3,c4)
      P5 <- replace(Puts2[,5],which(Puts2[,5]!=0),-10)
      P6 <- replace(Puts2[,6],which(Puts2[,6]!=0),-10)
      P7 <- replace(Puts2[,7],which(Puts2[,7]!=0),-10)
      P8 <- replace(Puts2[,8],which(Puts2[,8]!=0),-10)
      P5 <- replace(P5,which(P5==0),1)
      P6 <- replace(P6,which(P6==0),1)
      P7 <- replace(P7,which(P7==0),1)
      P8 <- replace(P8,which(P8==0),1)
      Puts3 <- cbind(Puts,P5,P6,P7,P8)
      Puts4 <- cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
      alleles2 <- OUTALL[which(OUTALL[,1]==L),]
      alfreq1 <- alleles2[match(Puts4[,1],alleles2[,2]),4]
      alfreq2 <- alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
      alfreq3 <- alleles2[match(Puts4[,3],alleles2[,2]),4]
      alfreq4 <- alleles2[match(Puts4[,4],alleles2[,2]),4]
      Puts5 <- cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
      R1 <- replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
      R2 <- replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
      R3 <- replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
      R4 <- replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
      Puts6 <- cbind(R1,R2,R3,R4)
      Put_share <- apply(Puts6, 1, min)                                          #find row minimum
      Put_share2 <- apply(Puts4, 1, max)                                         #find shared allele name
      Put_share3 <- c(Put_share,Put_share2)
      OBS <<- cbind(OBS,Put_share3)
    }
    L <- ncol(Putative2)-2
    C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths <- length(OBS[,1])/2
    if (lengths==1) {OBA3 <- t(OBS[(lengths+1):(2*lengths),])} else {OBA3 <- OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file="True_shared_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS <- t(OBS[1:lengths,])} else {OBS <- OBS[1:lengths,]}   #formatting for if there is only a single pair
    obsp <- apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
    OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
  }
  #calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#
  
  OBA3<- read.table("True_shared_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
      vect <- c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
      alleles <- data.frame(table(vect))
      alleles <- alleles[order(alleles[,1]),]
      if (as.numeric(as.character(alleles[1,1]))<=0) {alleles <- alleles[-1,]}
      alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))
      gtrue <- sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
      
      if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
        mm <- alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
        write.table(cbind(r,mm),file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
      }
      }
    }
    for (i in index) {
      if (i==0) {DIST <- as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                   #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distt<-cbind(i,DIST)
      write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  #Calculate lamdaphi=============================================================#
  
  Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Putadults <- Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
  Putoffs <- Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
  Names <- cbind(as.character(Putadults),as.character(Putoffs))
  empirical <- cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
  distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  P2 <- NULL
  for (i in index) {                                                              #loop over numbers of mismatched loci
    empirical2 <- empirical[which(as.numeric(as.character(empirical[,3]))==i),]
    if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
    if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
      a3 <- distp[which(distp[,1]==i),2]
      DIST<-as.data.frame(a3)
      P <- NULL
      for (b in 1:length(empirical2[,1])) {
        p1 <- length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
        if (p1==0) {p1=0.00001}
        p2 <- cbind(empirical2[b,1],empirical2[b,2],p1)
        p3 <- cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
        P <- rbind(P,p3)
      }
    }
    P2<-rbind(P2,cbind(i,P))
  }
  lamdaphi <- as.numeric(as.character(P2[,5]))
  #Calculate lamda|phic ==========================================================#
  
  lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  mm1 <- which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
  Observed[mm1,2] <- 1                                                            #replace NAs with 1
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  lamdaphic <- apply(U, 1, prod)
  l1 <- length(which(OUT9[,2]==0))
  if (l1>0) lamdaphic <- c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
  P3 <- cbind(P2,lamdaphic)
  for (i in index) {                                                              #loop over numbers of mismatched loci
    e2 <- P3[which(as.numeric(as.character(P3[,1]))==i),]
    if(length(e2)==6) {e2 <- t(e2)}                                                 #deals with one indiviudal formatting
    if (e2[1,5]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
      a3 <- lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
      DIST<-as.data.frame(a3)
      for (b in 1:length(e2[,1])) {                                                 #calculate p values
        p1 <- length(which(DIST[,1] <= e2[b,6]))
        p2 <- p1/length(DIST[,1])
        write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      }
    }
  }
  lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  #Put it all together with Bayes theorem!========================================#
  
  vals <- cbind(P2,lamdaphic[,1])
  philength <- cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
  phis <- rep(philength[,2],philength[,3])
  vals <- cbind(vals,phis)
  colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
  phi <- as.numeric(as.character(vals[,7]))
  lamdaphi <- as.numeric(as.character(vals[,5]))
  lamdaphic <- as.numeric(as.character(vals[,6]))
  lamdaphi <- replace(lamdaphi,which(lamdaphi==0),1)
  lamdaphic <- replace(lamdaphic,which(lamdaphic==0),1)
  pval <- (lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
  pval <- cbind(vals[,2],vals[,3],vals[,1],pval)
  pval <- pval[order(as.numeric(pval[,4])),]
  colnames(pval) <- c("Adult","Offspring","NL_mismatch","Probability of pair being false given frequencies of shared alleles")
  #write.table(vals, file="Posterior_Components_Bayes.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  
  unlink("False_allele_freqs.txt")                                                #clear all sims files
  unlink("True_allele_freqs.txt")
  unlink("Shared_allele_freqs.txt")
  unlink("lamdaphic.txt")
  unlink("Putative.txt")
  unlink("True_shared_freqs.txt")
  unlink("Output_genotypes.txt")
  unlink("*.sims")
  rm(list=ls())
} #end of exclusion


#########################
### no.par.bayes.ns() ###
#########################

#function the same as SOLOMON's no.par.bayes(), but with the slight change to the script to make 
#computation time faster


no.par.bayes.ns <- function(adults,offspring,reps=1000,Ngts=50000000)     {
  
  #loading in adults and offspring genotypes
  Adults1 <- adults
  Offspring1 <- offspring
  
  #making sure the number of reps is numeric and counting the number of loci
  reps <- as.numeric(reps)
  loci <- ncol(Adults1)
  
  #removing ID column
  Adults <- Adults1[,c(2:loci)]                                                    
  Offspring <- Offspring1[,c(2:loci)]
  
  #getting the number loci, again for Progress bar
  total <- ncol(Adults)/2                                                       
  
  #putting in a progress bar
  
  
  #Begin Master simulation function
  afreqs <- function(afreqs){
    
    #making locus names? *** need to figure out ** looks like used in progress bar
    locus_name <- L
    
    #setting up progress bar again
    
    
    #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
    vect <- c(Adults[,L],Adults[,L+1])                                                 
    alleles <- data.frame(table(vect))
    alleles <- alleles[order(alleles[,1]),]
    if (as.numeric(as.character(alleles[1,1]))==0) {alleles <- alleles[-1,]}           #removes missing data
    if(length(alleles[,1])==1) {                                                    #deals with monomorphic loci by adding 1 very strange allele (799)
      alleles <- cbind(vect[1],alleles[2])
      alleles[2,]<-c(799,1)}
    alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
    homos <- (alleles2[,3])^2                                                          #create homozygote allele frequencies
    homos2 <- cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
    hets <- t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
    hetfreq <- 2*(hets[,1]*hets[,2])
    hetvals <- t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
    hets2 <- cbind(hetvals,hetfreq)
    gfreqs <- rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
    csum <- cumsum(as.numeric(gfreqs[,3]))
    gfreqs1 <- cbind(gfreqs,csum)
    Nadults <- length(Adults[,1])
    Noffs <- length(Offspring[,1])
    
    
    #===============================================================================#end locus-specific HWE genotype frequency calculations
    alength <- length(alleles2[,1])
    for (Y in 1:reps) {
      positions <- 1:length(gfreqs[,1])
      sg3 <- sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                      #sample the repeated genotype positions, by the number of adults
      sadults <- gfreqs[sg3,1:2]                                                         #index gfreqs to create genotypes
      og3 <- sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                        #create juvenile genotyes
      soffs <- gfreqs[og3,1:2]
      soffs <- cbind(as.numeric(soffs[,1]),as.numeric(soffs[,2]))
      asp <- cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
      asp <- rbind(cbind(locus_name,0,0),asp)
      for (i in 1:Nadults) {
        parent1 <- as.numeric(sadults[i,1])                                                #first allele in parent
        parent2 <- as.numeric(sadults[i,2])                                                #second allele in parent
        p1 <- soffs-parent1
        p2 <- soffs-parent2
        pp1 <- which(p1[,1]==0)
        pp2 <- which(p1[,2]==0)
        allele1 <- unique(c(pp1,pp2))
        p21 <- which(p2[,1]==0)
        p22 <- which(p2[,2]==0)
        allele2 <- unique(c(p21,p22))
        Out51 <- cbind(parent1,length(allele1))
        Out52 <- cbind(parent2,length(allele2))
        Out53 <- cbind(0,Noffs-length(unique(c(allele1,allele2))))
        Out5 <- rbind(Out51,Out52,Out53)
        if(parent2==parent1) {Out5 <- Out5[-1,]}                                           #remove 1 of alleles count if homozygous
        if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
          diffs <- sum(Out5[,2])-Noffs                                                     #comment out to be more conservative!
          maxa <- max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
          pos <- which(Out5[,2]==maxa)
          Out5[pos,2]<-Out5[pos,2]-diffs}
        m1 <- match(Out5[,1],asp[,2])
        m2 <- asp[m1,3]+as.numeric(Out5[,2])
        asp[m1,3]<-m2
        asp<-asp
      }
      write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  L <- ncol(Adults)
  C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
  
  #Bayes P-value==================================================================#
  OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  locname <- unique(OUT[,1])                                                         #compile calculations for each locus
  OUTALL <- NULL
  for (z in locname) {
    Loc1 <- OUT[which(OUT[,1]==z),]
    allfreqs <- unique(Loc1[,2])
    OUT2 <- NULL
    for (x in allfreqs) {
      a1<-Loc1[which(Loc1[,2]==x),]
      a2 <- sum(a1[,3])
      a3 <- cbind(x,a2)
      OUT2 <- rbind(OUT2, a3)
    }
    OUT3 <- cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
    OUTALL <- rbind(OUTALL, cbind(z,OUT3))
  }
  #Create multilocus genotypes, 50000 at a time, calculate number of loci that mismatch, and calculate freqs of shared alleles
  NL <- length(unique(OUTALL[,1]))
  ngtypes <- 1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
  Ngts <- as.numeric(Ngts)
  inreps <- 50000     #was tested as 10000 for SNPS
  repnumber <- round(Ngts/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  asp <- cbind(0:NL,rep(0,length(0:NL)))
  
  for (n in 1:repnumber) {
    
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    distm <- apply(OUT, 1, function(x)sum(x == 1))
    distm2 <- data.frame(table(distm))
    m1 <- match(distm2[,1],asp[,1])
    m2 <- asp[m1,2]+distm2[,2]
    asp[m1,2]<-m2
    asp<-asp
  }
  
  #tabulate number of multilocus genotypes with x mismatches
  d2 <- asp
  d3 <- cbind(d2,d2[,2]/sum(d2[,2]))
  #===============================================================================# Create plot of exclusionary power.  Also, some necessary data formatting
  
  
  Adults<- adults
  Offs<- offspring
  Nads <- length(Adults[,1])
  Noffs <- length(Offs[,1])
  asp <- d3
  asp <- cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
  asp <- cbind(asp,cumsum(asp[,2]))
  asp <- cbind(asp,asp[,5]/max(asp[,5]))
  #find minimum Nloci to mismatch (could modify this) and perform exclusion with decided-upon mismatches==#
  distm <- cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
  #mismatch=min(which(round(distm[,6],1)==.9))                                    #deprecated
  mismatch <- which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
  Adults1 <- Adults                                                                   #begin exclusion
  Offspring1 <- Offs
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]
  write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
  IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  matches <- function(matches)
  {
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)
    write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  }
  z <- ncol(Ads)
  C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
  Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)
  IDS <- which(stuff<(mismatch+1))
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"
  head(PAdults)
  head(POffspring)
  PAdults$sorter <- as.numeric(row.names(PAdults))
  POffspring$sorter <- as.numeric(row.names(POffspring))
  Putative <- rbind(PAdults,POffspring)
  Putative <- Putative[order(Putative$sorter),]
  Putative$sorter <- NULL
  #head(Putative)
  #sorts <- function(sorts){
  #  tell <- rbind(PAdults[f,],POffspring[f,])
  #  write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  #}
  #f <- length(PAdults[,1])
  #C1 <-  for(f in 1:f) lapply(f,sorts)
  
  unlink("Sort.txt")
  unlink("IdnamesA.txt")
  unlink("IdnamesO.txt")
  #Calculate phi for each number mismatching loci=================================#
  
  #Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  
  observed <- data.frame(table(Putative[,2]))
  observed <- cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
  zerom <- 0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
  zerom2 <- which(is.na(match(zerom,observed[,1])))
  if (length(zerom2>0)) {observed <- observed[,3]
                         for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}  #not really 0, to prevent divide by 0
  }   else {observed=observed[,3]}
  expected <- distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
  #expected=distm[1:(mismatch+1),3]*Nads*Noffs                                     #not using cumulative sum
  phi <- expected/observed
  phi <- replace(phi,which(phi>=1),1)
  Offs<- offspring
  actualTrue <- length(grep("Off",Offs[,1]))
  info <- cbind(actualTrue,expected,observed,phi)
  #calculate phi and index values ================================================#
  
  phibase <- phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
  observed <- observed[min(which(phi[]==1))-1]                                       #do the same for observed
  testob <- which(observed==0.000000001)
  phi2 <- cbind(1:length(phi),phi)
  if (length(testob>0)) {phi4 <- phi2[-testob,]} else {phi4 <- phi2}
  nmismatch <- min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
  index <- phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
  index <- index[which(index>-1)]
  if (length(index)>1) {
    if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
  phi <- phi[index]
  index <- index-1
  #Create Plot ===================================================================#
  pdf(file <- "Output_Dataset_Power.pdf")
  x <- 0:(length(info[,1])-1)
  y1 <- info[,3]
  y2 <- info[,2]
  p1 <- which(y2==0)
  y2 <- replace(y2,p1,.000000001)
  p1 <- which(y1<y2)
  y3 <- replace(y1,p1,y2[p1])
  y2 <- y2+1
  y3 <- y3+1
  par(mar = c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
  ats <- c(0,1,2,3,4,5,6,7,8,9)
  labs <- c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
  axis(side=2,ats,labs)
  lines(x,log10(y3),lwd=2)
  lines(x,log10(y2),lwd=2)
  points(x,log10(y3),pch=21,bg="green",cex=2)
  points(x,log10(y2),pch=21,bg="blue",cex=2)
  legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
  yphi <- y2/y3
  par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
  lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
  points(x,yphi,cex=2,pch=21,bg="gray")
  mtext("Number of Mismatching Loci",side=1,line=1.94)
  dev.off()
  info2 <- cbind(x,info[,4])
  colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
  write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  #===============================================================================#
  ngtypes <- 20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
  inreps <- 10000
  repnumber <- round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  #writes values at all loci, to be analyzed further below
  for (n in 1:repnumber) {
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      alleles3 <- alleles3[-which(alleles3[,2]==0),]
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    for (i in index) {                                                              #loop over numbers of mismatched loci
      if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                    #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distp <- cbind(i,DIST)
      write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  
  #Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
  OUT9 <- NULL
  for (n in index){
    Putative2 <- Putative[which(Putative[,2]==n),]
    write.table(Putative2, file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                                #need at least 1000 values , else assigned 0 (may be redundant now)
    OBS <- NULL
    afreqs <- function(afreqs) {
      PutL <- Putative2[,c(L+2,L+3)]
      PutLadults <- PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
      PutLoffs <- PutL[seq(from=2,to=length(PutL[,1]),by=2),]
      Puts <- cbind(PutLadults,PutLoffs)
      c1 <- Puts[,1]-Puts[,3]                                                    #find the matching alleles
      c2 <- Puts[,1]-Puts[,4]
      c3 <- Puts[,2]-Puts[,3]
      c4 <- Puts[,2]-Puts[,4]
      Puts2 <- cbind(Puts,c1,c2,c3,c4)
      P5 <- replace(Puts2[,5],which(Puts2[,5]!=0),-10)
      P6 <- replace(Puts2[,6],which(Puts2[,6]!=0),-10)
      P7 <- replace(Puts2[,7],which(Puts2[,7]!=0),-10)
      P8 <- replace(Puts2[,8],which(Puts2[,8]!=0),-10)
      P5 <- replace(P5,which(P5==0),1)
      P6 <- replace(P6,which(P6==0),1)
      P7 <- replace(P7,which(P7==0),1)
      P8 <- replace(P8,which(P8==0),1)
      Puts3 <- cbind(Puts,P5,P6,P7,P8)
      Puts4 <- cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
      alleles2 <- OUTALL[which(OUTALL[,1]==L),]
      alfreq1 <- alleles2[match(Puts4[,1],alleles2[,2]),4]
      alfreq2 <- alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
      alfreq3 <- alleles2[match(Puts4[,3],alleles2[,2]),4]
      alfreq4 <- alleles2[match(Puts4[,4],alleles2[,2]),4]
      Puts5 <- cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
      R1 <- replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
      R2 <- replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
      R3 <- replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
      R4 <- replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
      Puts6 <- cbind(R1,R2,R3,R4)
      Put_share <- apply(Puts6, 1, min)                                          #find row minimum
      Put_share2 <- apply(Puts4, 1, max)                                         #find shared allele name
      Put_share3 <- c(Put_share,Put_share2)
      OBS <<- cbind(OBS,Put_share3)
    }
    L <- ncol(Putative2)-2
    C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths <- length(OBS[,1])/2
    if (lengths==1) {OBA3 <- t(OBS[(lengths+1):(2*lengths),])} else {OBA3 <- OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file="True_shared_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS <- t(OBS[1:lengths,])} else {OBS <- OBS[1:lengths,]}   #formatting for if there is only a single pair
    obsp <- apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
    OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
  }
  #calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#
  
  OBA3<- read.table("True_shared_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
      vect <- c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
      alleles <- data.frame(table(vect))
      alleles <- alleles[order(alleles[,1]),]
      if (as.numeric(as.character(alleles[1,1]))<=0) {alleles <- alleles[-1,]}
      alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))
      gtrue <- sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
      
      if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
        mm <- alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
        write.table(cbind(r,mm),file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
      }
      }
    }
    for (i in index) {
      if (i==0) {DIST <- as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                   #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distt<-cbind(i,DIST)
      write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  #Calculate lamdaphi=============================================================#
  
  Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Putadults <- Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
  Putoffs <- Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
  Names <- cbind(as.character(Putadults),as.character(Putoffs))
  empirical <- cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
  distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  P2 <- NULL
  for (i in index) {                                                              #loop over numbers of mismatched loci
    empirical2 <- empirical[which(as.numeric(as.character(empirical[,3]))==i),]
    if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
    if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
      a3 <- distp[which(distp[,1]==i),2]
      DIST<-as.data.frame(a3)
      P <- NULL
      for (b in 1:length(empirical2[,1])) {
        p1 <- length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
        if (p1==0) {p1=0.00001}
        p2 <- cbind(empirical2[b,1],empirical2[b,2],p1)
        p3 <- cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
        P <- rbind(P,p3)
      }
    }
    P2<-rbind(P2,cbind(i,P))
  }
  lamdaphi <- as.numeric(as.character(P2[,5]))
  #Calculate lamda|phic ==========================================================#
  
  lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  mm1 <- which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
  mm1
  Observed[mm1,2] <- 1                                                            #replace NAs with 1
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  lamdaphic <- apply(U, 1, prod)
  l1 <- length(which(OUT9[,2]==0))
  if (l1>0) lamdaphic <- c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
  P3 <- cbind(P2,lamdaphic)
  for (i in index) {                                                              #loop over numbers of mismatched loci
    e2 <- P3[which(as.numeric(as.character(P3[,1]))==i),]
    if(length(e2)==6) {e2 <- t(e2)}                                                 #deals with one indiviudal formatting
    if (e2[1,5]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
      a3 <- lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
      DIST<-as.data.frame(a3)
      for (b in 1:length(e2[,1])) {                                                 #calculate p values
        p1 <- length(which(DIST[,1] <= e2[b,6]))
        p2 <- p1/length(DIST[,1])
        write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      }
    }
  }
  lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  #Put it all together with Bayes theorem!========================================#
  
  vals <- cbind(P2,lamdaphic[,1])
  philength <- cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
  phis <- rep(philength[,2],philength[,3])
  vals <- cbind(vals,phis)
  colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
  phi <- as.numeric(as.character(vals[,7]))
  lamdaphi <- as.numeric(as.character(vals[,5]))
  lamdaphic <- as.numeric(as.character(vals[,6]))
  lamdaphi <- replace(lamdaphi,which(lamdaphi==0),1)
  lamdaphic <- replace(lamdaphic,which(lamdaphic==0),1)
  pval <- (lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
  pval <- cbind(vals[,2],vals[,3],vals[,1],pval)
  pval <- pval[order(as.numeric(pval[,4])),]
  colnames(pval) <- c("Adult","Offspring","NL_mismatch","Probability of pair being false given frequencies of shared alleles")
  #write.table(vals, file="Posterior_Components_Bayes.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  
  unlink("False_allele_freqs.txt")                                                #clear all sims files
  unlink("True_allele_freqs.txt")
  unlink("Shared_allele_freqs.txt")
  unlink("lamdaphic.txt")
  unlink("Putative.txt")
  unlink("True_shared_freqs.txt")
  unlink("Output_genotypes.txt")
  unlink("*.sims")
  rm(list=ls())
}

#################
### sim.gt.data2()
#################

#Description
# This fuction creates genotypic data for parents and offspring using a breeding matrix

# Input parameters
# 1) alfs - a simple data frame first column 1:number of loci and 2 the frequecies of alleles at each locus
# 2) breeding matrix - three column data frame with mom id, dad id, and RS for each pair
# 3) error - defining the amount of random error in the dataset
# 4) Nunrelated - defining the number of unrelated individuals wanted in the dataset
# 5) write2file - defualt is set to FALSE and returns data to console, if TRUE writes data to file in current working directory

#Defining function
sim.gt.data2 <-  function(alfs, breeding.matrix, error= 0, Nunrelated= 0, write2file = F){
  
  #breeding.matrix names
  breeding.matrix [,1] <- paste("Mom",breeding.matrix[,1])
  breeding.matrix [,2] <- paste("Dad",breeding.matrix[,2])
  
  #figuring our how many parents are here
  moms <- unique(breeding.matrix[,1])
  dads <- unique(breeding.matrix[,2])
  par.names <- c(moms,dads)
  Nadults <- sum(length(moms),length(dads))
  
  #getting afreqs input
  afreqs <-alfs
  
  #making a function to simulate genotypes for parents
  
  OUT <- NULL
  sims <- function(sims){
    
    #table allele frequencies
    #getting just the frequencies changing generic name to three character genotypes
    alleles2 <- afreqs[afreqs[,1] == loci[i],]
    alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                  
    
    #create homozygote allele frequencies
    homos <- (alleles3[,2])^2                                                          
    homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
    
    #create heterozygote allele frequencies
    #first create all possible combinations of freqs now make expected freqs
    hets <- t(combn(alleles3[,2],2))                                
    hetfreq <- 2*(hets[,1]*hets[,2])
    
    #create heterozygote allele names
    #now getting the het. genotype combinations and combining them
    hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
    hets2 <- cbind(hetvals,hetfreq)
    
    #combine hets and homos and create genotypes
    gfreqs <- rbind(hets2,homos2)                                                      
    
    #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
    n <- 1000000                                    
    
    #create genotypes(by coloumn, 1 for each allele)
    #for the first and second allele
    gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
    gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
    
    #combining them and mixing them up
    gtypes <- cbind(gfreqs1,gfreqs2)
    gtypes <- gtypes[sample(1:length(gtypes[,1]),replace= F),]
    
    #now getting a random sample of the gentoypes for the parents
    sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
    OUT<<-cbind(OUT,sg1)
    
  } # end of sim function
  
  #defining the number of loci there are
  loci.length <- length(unique(afreqs[,1]))
  loci <- unique(afreqs[,1])
  loci
  
  #doing the simulation process for each locus
  for(i in 1:length(loci)){
    lapply(1,sims)
  } # end for loop
  OUT
  #saving OUT into object called parents
  parents <- OUT
  
  parents <- as.data.frame(parents,stringsAsFactors = F)
  for(i in 1:ncol(parents)){parents[,i] <- as.numeric(parents[,i])}
  
  #adding names of parents 
  parents$names <- par.names
  
  i <- NULL
  out <- NULL
  for(i in 1:nrow(breeding.matrix)){
    
    #getting the parent ids of interest
    mom <- breeding.matrix[i,1]
    mom.id <-as.numeric(gsub("Mom ","",breeding.matrix[i,1]))
    
    dad <- breeding.matrix[i,2]
    dad.id <-as.numeric(gsub("Dad ","",breeding.matrix[i,2]))
    
    #now picking those parents gts
    mom.gts <- parents[parents$names == mom,-ncol(parents)]
    dad.gts <- parents[parents$names == dad,-ncol(parents)]
    
    #combining parents gts
    parent.gts <- rbind(mom.gts,dad.gts)
    
    #identifying where a new locus starts in the DF
    loci.locs <- seq(from= 1, to=ncol(mom.gts),by = 2)
    
    #making kids genotype
    j<- NULL
    out1 <- NULL 
    for(j in loci.locs){
      
      #mom gts sampled
      alleles <- parent.gts[1,j]
      alleles <- c(alleles, parent.gts[1,j+1])
      
      mom.allele <- sample(alleles,size = breeding.matrix[i,3],replace = T)
      
      #dad gts sampled
      alleles <- parent.gts[2,j]
      alleles <- c(alleles, parent.gts[2,j+1])
      dad.allele <- sample(alleles,size = breeding.matrix[i,3],replace = T)
      
      #combining parents alleles
      alleles <- cbind(mom.allele,dad.allele)
      
      #saving to out1
      out1 <- cbind(out1,alleles)
    }
    
    off.id <-paste("Offspring",paste(mom.id,dad.id, sep="."))
    off.id <- paste(off.id,1:breeding.matrix[i,3],sep=".")  
    out1 <- cbind(off.id,out1)  
    out <- rbind(out,out1)
  }
  
  Offs <- data.frame(out,stringsAsFactors = F)
  for(i in 2:ncol(Offs))(Offs[,i] <- as.numeric(Offs[,i]))
  
  #seperating the parents
  parents <- parents[,c(ncol(parents),1:(ncol(parents)-1))]
  
  Dads <- parents[grepl("Dad",parents[,1])==T,]
  Moms <- parents[grepl("Mom",parents[,1])==T,]
  
  #add error======================================================================#
  #first for the dads
  #error=0
  #getting the gts
  Dadsg <- Dads[,-1]
  
  #figuring out how many gts to make errors in
  ldad <- length(Dadsg)*error
  
  #randomly selecting where to put the errors
  pdad1 <- sample(1:length(Dadsg),ldad,replace= FALSE)
  
  #randomly selecting some gets to replace in the those locations
  pdad2 <- Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
  
  #actually putting in the error
  Dadsg2 <- replace(Dadsg,pdad1,pdad2)
  
  #replacing the old Dads gts
  Dads <- cbind(Dads[,1],Dadsg2)
  
  #doing the same thing for the kids
  Offsg=Offs[,-1]
  loff=length(Offsg)*error
  poff1=sample(1:length(Offsg),loff,replace=FALSE)
  poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
  Offsg2=replace(Offsg,poff1,poff2)
  Offs=cbind(Offs[,1],Offsg2)
  
  
  #############################
  ### unrelated individuals ###
  #############################
  
  #Add unrealated individuals <- ====================================================#
  
  
  if (Nunrelated>0){
    
    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nunrelated*3                                                            
    afreqs <- alfs
    OUT <- NULL
    
    #making the allele frequencies again, this way they are completely random and unrelated
    #doing the same thing as lines 34-113
    sims <- function(sims) {
      alleles2 <- afreqs[which(afreqs[,1] == z),]
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])
      
      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2                                                         
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
      
      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))                                                   
      hetfreq <- 2*(hets[,1]*hets[,2])
      
      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
      hets2 <- cbind(hetvals,hetfreq)
      
      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)                                                      
      
      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000         
      
      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
      
      OUT <<-cbind(OUT,sg1)
      
    } #end of sims function
    
    #now applying function
    z1 <- length(unique(afreqs[,1]))
    for(z in 1:z1) {lapply(z,sims)}
    
    #repeating lines 122-204
    parents <- OUT
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    parents <- as.numeric(parents)-Anumbs
    
    #to called this thing unrelated
    unrelated <- cbind(paste("Individual",1:length(parents[,1])),parents)
    Lid <- length(unrelated[,1])/3
    
    #splitting the unrelated individuals among dads, moms, and kids
    m1 <- unrelated[1:Lid,]
    d1 <- unrelated[(Lid+1):(Lid*2),]
    j1 <- unrelated[((Lid*2)+1):(Lid*3),]
    
    #adding the unrelated individuals to the moms, dads, and kids
    m1 <- data.frame(m1)
    colnames(m1) <- colnames(Moms)
    m1[,1] <- gsub(pattern = "Individual", replacement = "Ind.mom",x = m1[,1])
    Moms <- rbind(Moms,m1)
    for(i in 2:ncol(Moms)){Moms[,i] <- as.numeric(Moms[,i])}
    
    d1 <- data.frame(d1)
    colnames(d1) <- colnames(Dads)
    d1[,1] <- gsub(pattern = "Individual", replacement = "Ind.dad",x = d1[,1])
    Dads <- rbind(Dads,d1)
    for(i in 2:ncol(Dads)){Dads[,i] <- as.numeric(Dads[,i])}
    
    j1 <- data.frame(j1)
    colnames(j1) <- colnames(Offs)
    j1[,1] <- gsub(pattern = "Individual", replacement = "Ind.off",x = j1[,1])
    Offs <- rbind(Offs,j1)
    for(i in 2:ncol(Offs)){Offs[,i] <- as.numeric(Offs[,i])}
  }
  ### final touches
  
  #making generic locus names
  colsnmz <- rep("Locus",ncol(Offs)-1)
  colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
  colsnmz3 <- paste(colsnmz,colsnmz2)
  colsnmz3 <- c("IDs",colsnmz3)
  
  #usig them to make the colnames for the Offs, Dads, and Moms
  colnames(Offs)<- colsnmz3
  colnames(Dads)<- colsnmz3
  colnames(Moms)<- colsnmz3
  
  df <- rbind(Moms,Dads,Offs)
  
  p.count <- nrow(Moms)+nrow(Dads)
  k.count <- nrow(Offs)
  Type <- c(rep("Parents",times = p.count),rep("Kids",times = k.count))
  
  df <- data.frame(Type,df,stringsAsFactors = F)
  
  if(write2file==F){
    
    return(df)
    
  } else{
    
    write.table(Moms, "moms.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
    write.table(Dads, "dads.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
    write.table(Offs, "kids.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  }
}

####################
### genepop.format()
####################

#Description
#This fuction takes three character genotypic data with a pop and sample id column in the beginning and puts it into 
#a format acceptable to genepop

#Input Parameters:
#genos - a dataframe with two columns - "pop" and "sample.name", then the three character genotypic data
#title - the name of the data frame you are putting into genepop - Default: "Genotypic data"

genepop.format <- function(genos, title="Genotypic data"){
  
  #gts names for genepop
  genos[,1] <- paste(genos[,1], ",")
  head(genos)
  
  #now for the genepop prep
  my.sep <- genos
  my.sep[1,] <- c("Pop",rep("",(ncol(genos)-1)))
  my.sep <- my.sep[1,]
  my.sep
  
  #getting the unique pop identifiers
  my.ids <- unique(genos[,1])
  
  tmp2 <- NULL
  for(i in 1:length(my.ids)){
    x <- genos[genos[,1] == my.ids[i],]
    tmp2 <- rbind(tmp2,my.sep,x)
  }
  
  #just a check to make sure I have eight pops
  nrow(tmp2[tmp2[,1] == "Pop",])
  
  #making into six.char()
  tmp2 <- six.char(tmp2,missing.data="000",one.pop=F)
  head(tmp2)
  
  #removing the Pop ID
  tmp2 <- tmp2[,-2]
  head(tmp2)
  
  
  #prepping for the header for the genepop file
  loci <- c(paste0(colnames(tmp2)[-c(1, length(colnames(tmp2)))],", "),paste0(colnames(tmp2)[length(colnames(tmp2))]), "")
  loci
  
  head(tmp2)
  ncol(tmp2)
  my.sep <- tmp2
  my.sep[1,] <- c(title,rep("",(ncol(tmp2)-1)))
  my.sep <- my.sep[1,]
  my.sep
  
  tmp2 <- rbind(my.sep,loci,tmp2)
  
  
  print("Writing Genepop formatted data to current working directory")
  #writing to file now
  write.table(tmp2,"genepop.gts.txt",append=F,quote=F,sep="\t",row.names=F,col.names=F)
}

#######################
### all.freq.compar ###
#######################

#Description
#This fuction compares the allele frequencies of two groups to see if alleles are lost or retained or new. I created to compare
#genetic variation retrained in offspring compared to parents

#Input Parameters:
#gts - a dataframe with two columns - "Pops" and "sample.name", then the three character genotypic data
#grp1 - a character string to select the one group in the first column Pops, e.g. "Group 1"
#grp2 - a character string to select the two group in the second column Pops, e.g. "Group 2"

afreq.compar <- function(gts, grp1, grp2){
  
  #getting allele frequencies for all populations combined and individual populations
  ind.alfs <- alf.freq(population = gts,one.pop = F)
  all.alfs <- alf.freq(population = gts[,-2],one.pop=T)
  
  #remocinv some columns and adding acolumn
  head(all.alfs)
  all.alfs$frequency <- NULL
  all.alfs$count <- NULL
  all.alfs$parents <- 0
  
  
  #now doing some matching with the first group
  print("Analyzing genetic data from the first group")
  df1 <- filter(ind.alfs, Pops == grp1)
  head(df1)
  for(i in 1:nrow(all.alfs)){
    for(j in 1:nrow(df1)){
      if(all.alfs$Locus[i] == df1$Locus[j] & all.alfs$allele[i] == df1$allele[j]){
        all.alfs$parents[i] = df1$frequency[j]
      }
    }
  }
  
  #adding another column for the second group and repeating the project
  print("Analyzing genetic data from the second group")
  all.alfs$offspring <- 0
  df1 <- filter(ind.alfs, Pops == grp2)
  head(df1)
  for(i in 1:nrow(all.alfs)){
    for(j in 1:nrow(df1)){
      #    print(paste(i,j))
      if(all.alfs$Locus[i] == df1$Locus[j] & all.alfs$allele[i] == df1$allele[j]){
        all.alfs$offspring[i] = df1$frequency[j]
      }
    }
  }
  head(all.alfs)
  head(df1)
  
  all.alfs$alf.diff <- all.alfs$parents-all.alfs$offspring
  head(all.alfs)
  
  all.alfs$lost <- ifelse(all.alfs$parents > 0 & all.alfs$offspring == 0, "Lost","Retained")
  all.alfs$lost[all.alfs$parents == 0 & all.alfs$offspring >0]<-"New"
  table(all.alfs$lost)
  return(all.alfs)
}


#################
### obey.mendel()
#################

#Description
#This fuction evaluted if an given assignment based on SOLOMON obeys mendelian inheritance laws. It counts the number of mismatches (mm) between
#each parent and the offspring indivdiually and then together after taking into account mendelian inheritance, mean that both parents can't assign to a single
#allele at any given locus. In addition, it does the same process, but disregards missing data.

#Input Parameters:
#tmp - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs - must have all there IDs for a given assignment on one row
#gts - 
#opt - dictates what is returned to user
# opt 1 - currently the only option, which just returns the original tmp with 3 additional columns:
#      mom.mm - number of mms because mother and offspring, after account for missing data
#      dad.mm - number of mms betweeen father and offspring, after account for missing data
#      trio.mm - number of mms between offspring and both parents - takes into account mendelian inheritance, after account for missing data

assign <- NULL

obey.mendel <- function(assign,gts,opt=1){
  
  
  # definining three.char first
  #################
  ### three.char()
  #################
  
  #Description
  # This fuction takes genotypes in six character form and puts them into three character form
  
  #Input Parameters:
  
  # population - a table with sample.name in column one, and loci after that. All loci must be in six character form. Examples 000000, 111111, 999999
  
  three.char <- function(assign){
    
    #saving sample names into a seperate file
    OUT <- data.frame(assign[,1])
    colnames(OUT) <- "Sample.Name"
    
    #just grabbing all the genotypes with no sample names
    gts <- assign[,-1]
    
    #counting the number of loci and saving their names
    L <- ncol(gts)
    Locus.names <- colnames(gts)
    
    #using a small for loop to take each column and split it into two columns and cbinding that to OUT, which originally just had the sample names in it
    for(i in 1:ncol(gts)){
      locus <- data.frame(gts[,i])
      colnames(locus) <- Locus.names[i]
      
      locus$x1 <- substr(locus[,1],1,3)
      locus$x2 <- substr(locus[,1],4,6)
      locus[,1] <- NULL
      
      colnames(locus) <- c(paste(Locus.names[i],".1", sep=""),paste(Locus.names[i],".2",sep=""))
      OUT <- cbind(OUT,locus)
    }
    return(OUT)
  }
  
  j <- 3
  j <- NULL
  assign.out <- NULL
  
  for(j in 1:nrow(assign)){
    #for(j in x){  
    #getting each row seperately
    assign1 <- assign[j,]
    assign1
    
    #preparing data to merge with genotype file
    genos <- as.data.frame(t(assign1))
    colnames(genos) <- "Sample.Name"
    genos$order <- 1:nrow(genos)
    genos
    
    #merging genotypes
    genos <- merge(genos,gts,by="Sample.Name")
    genos <- genos[order(genos$order),]
    genos <- genos[,-2]
    
    #putting in a check to make sure that you have the correct number of genotypes to do the rest of the script
    if(nrow(genos) != 3) { stop(paste("You are missing the genotypes for one individual in assignment",j,"- check the names to see if they are spelled differently")) }
    
    three.char(genos)
    genos <- t(genos)
    geno.names <- genos[1,]
    geno.names 
    genos <- genos[-1,]
    genos <- as.data.frame(genos)
    colnames(genos) <- geno.names
    genos
    
    genos$mom.mm <- 0
    genos$dad.mm <- 0
    genos$trio.mm <- 0
    
    genos$mom.mID <- 0
    genos$dad.mID <- 0
    genos$trio.mID <- 0
    genos
    
    loci <- nrow(genos)
    loci
    
    #counting the number of matches and will use later to calculate the number of mismatches
    mom.counter <- 0
    dad.counter <- 0
    trio.counter <- 0
    
    genos
    
    #this will be in a for loop
    #i <- 6
    i <- NULL
    assign1
    genos
    
    for(i in 1:nrow(genos)){
      
      off.gt <- genos[i,1]
      mom.gt <- genos[i,2]
      dad.gt <- genos[i,3]
      
      #splitting the genotype into its alleles
      off.gt <- c(substr(off.gt,1,3),substr(off.gt,4,6))
      mom.gt <- c(substr(mom.gt,1,3),substr(mom.gt,4,6))
      dad.gt <- c(substr(dad.gt,1,3),substr(dad.gt,4,6))
      
      #checking them out
      off.gt
      mom.gt
      dad.gt
      
      assign1
      t(assign1)
      
      #identifying if any individuals are homozygotes - if the parents are I remove one of the alleles
      
      #for offspring
      
      #first determining the genotype is missing
      off.gt1 <- off.gt[1] == "000"
      off.gt2 <- off.gt[2] == "000"
      
      #getting an object with T/F 
      off.mID <- off.gt1 == T & off.gt2 == T
      
      #using logicals to define if homo.o is T\F
      if(off.mID == T){
        #making the homozygote identifier FALSE if the data is missing
        homo.o <- F
      } else {
        #if the logic isn't met then I ask the question 
        homo.o <- off.gt[1] == off.gt[2] 
      }
      
      if(homo.o ==  T){off.gt <- off.gt[1]}
      homo.o
      off.gt
      
      #for mothers
      
      #first determining the genotype is missing
      mom.gt1 <- mom.gt[1] == "000"
      mom.gt2 <- mom.gt[2] == "000"
      
      #getting an object with T/F 
      mom.mID <- mom.gt1 == T & mom.gt2 == T
      
      #using logicals to define if homo.m is T\F
      if(mom.mID  == T){
        #making the homozygote identifier FALSE if the data is missing
        homo.m <- F
      } else {
        #if the logic isn't met then I ask the question 
        homo.m <- mom.gt[1] == mom.gt[2] 
      }
      
      if(homo.m ==  T){mom.gt <- mom.gt[1]}
      homo.m
      mom.gt
      
      #for fathers
      
      #first determining the genotype is missing
      dad.gt1 <- dad.gt[1] == "000"
      dad.gt2 <- dad.gt[2] == "000"
      
      #getting an object with T/F 
      dad.mID <- dad.gt1 == T & dad.gt2 == T
      #using logicals to define if homo.d is T\F
      if(dad.mID == T){
        #making the homozygote identifier FALSE if the data is missing
        homo.d <- F
      } else {
        #if the logic isn't met then I ask the question 
        homo.d <- dad.gt[1] == dad.gt[2] 
      }
      
      if(homo.d == T){dad.gt <- dad.gt[1]}  
      homo.d
      dad.gt
      
      
      #identifying if the offspring and the parents have the same genotypes
      off.mom.check <- length(mom.gt[mom.gt %in% off.gt])
      off.mom.check
      
      off.dad.check <- length(dad.gt[dad.gt %in% off.gt])
      off.dad.check    
      
      
      homo.d
      homo.m
      homo.o
      
      ######################
      ### No Homozygotes ###
      ######################
      if(homo.o == F & homo.m == F & homo.d == F){
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        off.gt
        mom.gt
        dad.gt
        
        #if the mother and offspring have the same genotype I compare the father first
        if(off.mom.check == 2 & off.dad.check != 2){
          off.gt2 <- off.gt[(off.gt %in% dad.gt) == F]
          off2.mom.check <- length(mom.gt[(mom.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.mom.check > 0) {trio.counter <- trio.counter + 1}
        } #end off.mom.check loop
        
        #if the father and offspring have the same genotype I compare the mother first
        if(off.mom.check != 2 & off.dad.check == 2){
          off.gt2 <- off.gt[(off.gt %in% mom.gt) == F]
          off2.dad.check <- length(dad.gt[(dad.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.dad.check > 0) {trio.counter <- trio.counter + 1}
        }#end of off.dad.check loop
        
        off.mom.check
        off.dad.check
        #if neither of the genotypes are identical then I randomize who I pick first
        if(off.mom.check != 2 & off.dad.check != 2){
          
          #if mom nor dad match do nothing
          if(off.mom.check == 0 & off.dad.check == 0){trio.counter <- trio.counter}
          
          #if just dad matches do nothing
          if(off.mom.check == 0 & off.dad.check == 1){trio.counter <- trio.counter}
          
          #if just mom matches do nothing
          if(off.mom.check == 1 & off.dad.check == 0){trio.counter <- trio.counter}
          
          #if mom and dad match then do the logic below
          if(off.mom.check == 1 & off.dad.check == 1){
            
            mom.gt
            dad.gt
            off.gt
            trio.counter
            
            #To avoid a bias in heterozgotic individuals I randomize which parent I look at first to remove potential matches in genotypes
            my.pick <- sample(c("mom","dad"),1)
            
            if((my.pick == "mom") == T){
              #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
              off.gt2 <- off.gt[(off.gt %in% mom.gt) == F]
              off.gt2
              off2.dad.check <- length(dad.gt[(dad.gt %in% off.gt2) == T] > 0)
              off2.dad.check
              if(off.mom.check > 0 & off.dad.check > 0 & off2.dad.check > 0) {trio.counter <- trio.counter + 1}
              
            } else {
              
              #if the dad matched the offspring the genotype they matched at is removed and one is added to trio counter
              off.gt2 <- off.gt[(off.gt %in% dad.gt) == F]
              off.gt2
              off2.mom.check <- length(mom.gt[(mom.gt %in% off.gt2) == T] > 0)
              off2.mom.check
              if(off.mom.check > 0 & off.dad.check > 0 & off2.mom.check > 0) {trio.counter <- trio.counter + 1}
            } #end of my.pick logic
            trio.counter
          }#end of randomization loop
        }#end of mom and dad 1 loop
        
        #if all genotypes match then by default the trio.counter goes up by one
        if(off.mom.check == 2 & off.dad.check == 2){
          trio.counter <- trio.counter + 1
        }#end of all match loop
        
      } #end of No Homozygotes
      
      #######################
      ### All Homozygotes ###
      #######################
      if(homo.o == T & homo.m == T & homo.d == T){
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
        mom.check <- off.mom.check > 0
        dad.check <- off.dad.check > 0
        
        if(mom.check == T & dad.check == T) {trio.counter <- trio.counter + 1}
      } #end of All Homozygote
      
      #####################################
      ### Mother and Father Homozygotes ###
      #####################################
      
      if(homo.o == F & homo.m == T & homo.d == T){
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #To avoid a bias in heterozgotic individuals I randomize which parent I look at first to remove potential matches in genotypes
        my.pick <- sample(c("mom","dad"),1)
        
        if(my.pick == "mom"){
          #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
          off.gt2 <- off.gt[(off.gt %in% mom.gt) == F]
          off2.dad.check <- length(dad.gt[(dad.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.dad.check > 0) {trio.counter <- trio.counter + 1}
          
        } else {
          
          #if the dad matched the offspring the genotype they matched at is removed and one is added to trio counter
          off.gt2 <- off.gt[(off.gt %in% dad.gt) == F]
          off2.mom.check <- length(mom.gt[(mom.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.mom.check > 0) {trio.counter <- trio.counter + 1}
        } #end of my.pick logic
      } # end of Mother and Father Homozygotes
      
      ############################
      ### Offspring Homozygote ###
      ############################
      if(homo.o == T & homo.m == F & homo.d == F){
        
        #matching parent genotypes to offspring genotypes
        
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
        mom.check <- off.mom.check > 0
        dad.check <- off.dad.check > 0
        
        if(mom.check == T & dad.check == T) {trio.counter <- trio.counter + 1}
      } #end of Offspring Homozygote
      
      #########################
      ### Father Homozygote ###
      #########################
      homo.o 
      homo.m
      homo.d
      if(homo.o == F & homo.m == F & homo.d == T){
        
        off.mom.check
        off.dad.check
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #if the dad matched the offspring the genotype they matched at is removed and one is added to trio counter
        off.gt2 <- off.gt[(off.gt %in% dad.gt) == F]
        off.gt2
        off2.mom.check <- length(mom.gt[(mom.gt %in% off.gt2) == T] > 0)
        if(off.mom.check > 0 & off.dad.check > 0 & off2.mom.check > 0) {trio.counter <- trio.counter + 1}
        trio.counter
      } #end of Father Homozygote
      
      #########################
      ### Mother Homozygote ###
      #########################
      
      if(homo.o == F & homo.m == T & homo.d == F){
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
        off.gt2 <- off.gt[(off.gt %in% mom.gt) == F]
        off2.dad.check <- length(dad.gt[(dad.gt %in% off.gt2) == T] > 0)
        if(off.mom.check > 0 & off.dad.check > 0 & off2.dad.check > 0) {trio.counter <- trio.counter + 1}
      } #end of Mother Homozygote
      
      ########################################
      ### Offspring and Father Homozygotes ###
      ########################################
      
      if(homo.o == T & homo.m == F & homo.d == T){
        
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #determining if the father and offspring are homozgotes for the same allele
        if((off.gt[1] == dad.gt[1])==T){
          
          #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
          off.gt2 <- off.gt[1]
          off2.mom.check <- length(mom.gt[(mom.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.mom.check > 0) {trio.counter <- trio.counter + 1}
          
        } else {
          
          #if the mom matched the offspring the genotype they matched at is removed and one is added to trio counter
          off.gt2 <- off.gt[(off.gt %in% dad.gt) == F]
          off2.mom.check <- length(mom.gt[(mom.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.mom.check > 0) {trio.counter <- trio.counter + 1}
        } # end of father and offspring same homozygote check
      } #end of Offspring and Father Homozygotes
      
      ########################################
      ### Offspring and Mother Homozygotes ###
      ########################################
      if(homo.o == T & homo.m == T & homo.d == F){
        
        #matching parent genotypes to offspring genotypes
        if(off.mom.check > 0) {mom.counter <- mom.counter + 1}
        if(off.dad.check > 0) {dad.counter <- dad.counter + 1}
        
        #determining if the father and offspring are homozgotes for the same allele
        if((off.gt[1] == mom.gt[1])==T){
          
          #if the dad matched the offspring the genotype they matched at is removed and one is added to trio counter
          off.gt2 <- off.gt[1]
          off2.dad.check <- length(dad.gt[(dad.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.dad.check > 0) {trio.counter <- trio.counter + 1}
          
        } else {
          
          #if the dad matched the offspring the genotype they matched at is removed and one is added to trio counter
          off.gt2 <- off.gt[(off.gt %in% mom.gt) == F]
          off2.dad.check <- length(dad.gt[(dad.gt %in% off.gt2) == T] > 0)
          if(off.mom.check > 0 & off.dad.check > 0 & off2.dad.check > 0) {trio.counter <- trio.counter + 1}
        }# end of mother and offspring same homozygote check
      } # end of Offspring and Mother Homozygotes
      
      genos$mom.mm[i] <- off.mom.check
      genos$dad.mm[i] <- off.dad.check
      genos$trio.mm[i] <- loci - trio.counter
      
      if(mom.mID == T){
        genos$mom.mID[i] <- 1
        genos$trio.mID[i] <- 1
      }
      
      if(dad.mID == T){
        genos$dad.mID[i] <- 1
        genos$trio.mID[i] <- 1
      }
      
      
    } #end of geno loop
    
    
    genos 
    assign[2,]
    assign1
    
    mom.counter
    dad.counter
    trio.counter
    
    ############################  
    mom.mm <-  loci - mom.counter
    dad.mm <- loci - dad.counter
    trio.mm <- loci - trio.counter
    
    mom.mm
    dad.mm
    trio.mm
    
    sum(genos$mom.mID)
    sum(genos$dad.mID)
    sum(genos$trio.mID)
    
    #    assign1$mom.mm <- mom.mm
    #    assign1$dad.mm <- dad.mm
    #    assign1$trio.mm <- trio.mm
    assign1$mom.mmc <- mom.mm- sum(genos$mom.mID)
    assign1$dad.mmc <- dad.mm- sum(genos$dad.mID)
    assign1$trio.mmc <- trio.mm- sum(genos$trio.mID)
    assign1
    
    assign.out <- rbind(assign.out,assign1)
    
  }
  
  if(opt==1){
    return(assign.out)
  }
}


######################
### mendel.gt.compar()
######################

#Description
#This fuction takes output from obey.mm and gets the gts for all there individuals made in the assignment and aligns them in adjacent rows for visual inspection
#Input Parameters:
#tmp - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs - must have all there IDs for a given assignment on one row
#gts - a dataframe containing all the genotypes for all the indivudals involved in the assignment 

mendel.gt.compar <- function(df,gts){
  
  # definining three.char first
  #################
  ### three.char()
  #################
  
  #Description
  # This fuction takes genotypes in six character form and puts them into three character form
  
  #Input Parameters:
  
  # population - a table with sample.name in column one, and loci after that. All loci must be in six character form. Examples 000000, 111111, 999999
  
  three.char <- function(tmp){
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1])
    colnames(OUT) <- "Sample.Name"
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1]
    
    #counting the number of loci and saving their names
    L <- ncol(gts)
    Locus.names <- colnames(gts)
    
    #using a small for loop to take each column and split it into two columns and cbinding that to OUT, which originally just had the sample names in it
    for(i in 1:ncol(gts)){
      locus <- data.frame(gts[,i])
      colnames(locus) <- Locus.names[i]
      
      locus$x1 <- substr(locus[,1],1,3)
      locus$x2 <- substr(locus[,1],4,6)
      locus[,1] <- NULL
      
      colnames(locus) <- c(paste(Locus.names[i],".1", sep=""),paste(Locus.names[i],".2",sep=""))
      OUT <- cbind(OUT,locus)
    }
    return(OUT)
  }
  
  
  #for names, mms, and gts
  
  j<- NULL
  genos.out <- NULL
  for(j in 1:nrow(df)){
    
    assign1 <- df[j,]
    
    #getting mm information
    mm1 <- data.frame(assign1[1,c(4:6)])
    
    #getting the first row of data
    assign1 <- assign1[1,c(1:3)]
    
    #preping column names for rbind and putting in correct order
    colnames(mm1) <- colnames(assign1)
    mm1[2,]<- 1:3
    mm1
    
    #the rbind
    assign1 <- rbind(assign1,mm1)
    
    #preparing data to merge with genotype file
    genos <- t(assign1)
    colnames(genos) <- c("Sample.Name","mms","sorter")
    genos
    
    #merging genotypes and sorter correctly
    genos <- merge(genos,gts,by="Sample.Name")
    genos <- genos[order(genos$sorter,decreasing=T),]
    
    #preping for thre.char()
    genos$sorter <- NULL
    mms <- genos$mms
    genos$mms <- NULL
    
    #running it and adding information back in
    genos <- three.char(genos)
    genos$assignment <- paste(row.names(genos),(j), sep=".")
    genos$mms <- mms
    
    genos.out <- rbind(genos.out,genos)
  }
  
  return(genos.out)
  
}

