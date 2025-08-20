###This is the final code used to calculate the gross growth rates for L2 and L6
rm(list=ls(all=TRUE))

Dir.o = "./output/100-growth-rates/"

#Set working directory and load libraries & scripts:
  # setwd("/Users/blazewicz1/Documents/manuscripts/Heavy_Water_wetup_2/qSIP_Mar18")
  library(reshape2)
  library(VennDiagram)
  source("./scripts/qSIP_repo/sample.vec.R")                      #sample.vec
  source("./scripts/qSIP_repo/WAD.func.R")                        #WAD.func
  source("./scripts/qSIP_repo/fit.norm.func.R")                   #fit.norm.func
  source("./scripts/qSIP_repo/boot.WAD.func.R")                   #boot.WAD.func
  source("./scripts/qSIP_repo/diff.wad.calc.R")                   #diff.wad.calc
  source("./scripts/qSIP_repo/boot.diff.wad.R")                   #boot.diff.wad
  source("./scripts/qSIP_repo/MW.calc.R")                         #MW.calc
  source("./scripts/qSIP_repo/MW.calc.Schildkraut.R")             #MW.calc.Schildkraut
  source("./scripts/qSIP_repo/comparison.message.R")              #comparison.message
  source("./scripts/qSIP_repo/ape.calc.lh.R")                        #ape.calc
  source("./scripts/qSIP_repo/boot.diff.ape.R")                   #boot.diff.ape
  source("./scripts/qSIP_repo/r.calc.R")                          #r.calc
  source("./scripts/qSIP_repo/boot.diff.r.R")                     #boot.diff.r
  source("./scripts/qSIP_repo/boot.TUBE.func.R")                  #boot.TUBE.func
  source("./scripts/qSIP_repo/f.calc.R")                          #f.calc
  source("./scripts/qSIP_repo/boot.diff.f.R")                     #boot.diff.f
  source("./scripts/qSIP_repo/all.taxa.calcs.R")                  #all.taxa.calcs
  source("./scripts/qSIP_repo/boot.pop.R")                        #boot.pop
  source("./scripts/qSIP_repo/boot.TUBE.pop.R")                   #boot.TUBE.pop
  source("./scripts/qSIP_repo/comparison.message.pop.R")          #comparison.message.pop
  source("./scripts/qSIP_repo/r.calc.pop.R")                      #r.calc.pop
  source("./scripts/qSIP_repo/f.calc.pop.R")                      #f.calc.pop
  source("./scripts/qSIP_repo/boot.r.pop.R")                      #boot.r.pop
  source("./scripts/qSIP_repo/boot.f.pop.R")                      #boot.f.pop
  source("./scripts/qSIP_repo/all.taxa.calcs.pop.R")              #all.taxa.calcs.pop
  source("./scripts/qSIP_repo/id.reps.R")                         #id.reps
  source("./scripts/qSIP_repo/select.rep.R")                      #select.rep
  source("./scripts/qSIP_repo/explore.filter.taxa.R")             #explore.filter.taxa
  source("./scripts/qSIP_repo/filter.taxa.R")                     #filter.taxa
  source("./scripts/qSIP_repo/explore.filter.fractions.taxa.R")   #explore.filter.fractions.taxa
  source("./scripts/qSIP_repo/filter.fractions.taxa.R")           #filter.fractions.taxa  
  source("./scripts/qSIP_repo/WAD.by.taxon.func.R")               #WAD.by.taxon.func
  source("./scripts/qSIP_repo/SE.WAD.by.taxon.plot.R")            #SE.WAD.by.taxon.plot
  source("./scripts/qSIP_repo/find.unlabeled.correction.R")       #find.unlabeled.correction
  source("./scripts/qSIP_repo/find.labeled.correction.R")         #find.labeled.correction
  source("./scripts/qSIP_repo/td.pos.resid.R")                    #td.pos.resid
  source("./scripts/qSIP_repo/td.abs.resid.R")                    #td.abs.resid
  source("./scripts/qSIP_repo/bu.abs.resid.R")                    #bu.abs.resid
  source("./scripts/qSIP_repo/select.best.iteration.R")           #select.best.iteration
  source("./scripts/qSIP_repo/find.labeled.correction.plot.R")    #find.labeled.correction.plot
  source("./scripts/qSIP_repo/get.seq.taxa.nums.R")               #get.seq.taxa.nums
  source("./scripts/qSIP_repo/add.lab.WAD.corr.summary.R")        #add.lab.WAD.corr.summary
  source("./scripts/qSIP_repo/apply.unlabeled.correction.R")      #apply.unlabeled.correction
  source("./scripts/qSIP_repo/apply.labeled.correction.R")        #apply.labeled.correction
  source("./scripts/qSIP_repo/U.sens.func.R")                     #U.sens.func
  source("./scripts/qSIP_repo/U.sens.plot.func.R")                #U.sens.plot.func


#APE CALCULATIONS:


#Import raw data and format raw data for analysis:
  #NOTE: treatment code names cannot contain spaces
  data.all <- read.table("./data/ASV_table_relabun_100percent_with_sample_info_total_copies.txt", header=T, sep="\t", stringsAsFactors=T, check.names=F) #Note that this input has total copies of 16S in original binned SIP fractions to account for the different volumes - See line 435 for further explanation
  summary(data.all)
  head(data.all)
  names(data.all)
  dim(data.all)


#Import taxa legend:
  taxa.id <- read.table("./data/taxonomy_table_reads.txt", header=F, sep="\t", stringsAsFactors=T, check.names=F, na.strings = "", col.names = c("taxon","code"))
  summary(taxa.id)
  head(taxa.id)
  names(taxa.id)
  dim(taxa.id)


#Create a working copy of the data:
  data <- data.all




#Create a column of unique tube IDs (where a unique tube is a combination of Hour, Isotope, Replicate, and I added treatment for 4th wedge)
  #unique.tube <- paste(data$Replicate, data$Hour, data$Isotope, sep="_")
  unique.tube <- paste(data$Replicate, data$Hour, data$Isotope, sep="_")
  data <- data.frame(data[1:8], unique.tube, data[9:dim(data)[2]])

#Create a column of unique treatment IDs (where a unique treatment is a combination of Hour, Isotope, and I added treatment for 4th wedge)
  unique.tmt <- paste(data$Hour, data$Isotope, sep="_")
  data <- data.frame(data[1:9], unique.tmt, data[10:dim(data)[2]])

  data[1:13]

#Add a column for the sum of total (normalized for read depth) coverages of all taxa by fraction
  data$sum.abundance <- rowSums(data[,13:ncol(data)])
  data$sum.abundance

#Calculate mass of DNA per fraction, based on relative abundance (relative coverage) and total amount of DNA (ng) per fraction:
  copies <- data$copies.ul*(data[,13:(ncol(data)-1)])
  dim(copies)
  copies <- cbind(data[,1:12], copies)  # add first 12 columns of data to ng.dna
  dim(copies)
  head(copies)
  names(copies)

#Melt data into long format by metadata;
#Do this for ng.dna and for raw coverage number, which is just our data file. Merge these to into 1 masterfile: data.melted
  copies.melted <- melt(copies, id=c("SampleID","tube","fraction","Replicate", "Hour","Isotope","Trt.code","unique.tube","unique.tmt","Density","copies.ul"), measure.vars=as.character(taxa.id$taxon), variable.name="taxon", value.name="t.copies.ul")
  rel.abundance.melted <- melt(data, id=c("SampleID","tube","fraction","Replicate", "Hour","Isotope","Trt.code","unique.tube","unique.tmt","Density","copies.ul"), measure.vars=as.character(taxa.id$taxon), variable.name="taxon", value.name="rel.abundance")
  data.melted <- merge(copies.melted, rel.abundance.melted)
  head(data.melted)
  dim(copies.melted)
  dim(rel.abundance.melted)
  dim(data.melted)

#Merge taxa data and reorder data frame by taxon and tube and fraction
  data.melted <- merge(data.melted, taxa.id)
  data.melted <- data.melted[order(data.melted$taxon, data.melted$tube, data.melted$fraction),]
  row.names(data.melted) <- 1:dim(data.melted)[1]     #rename observations to be sequential
  head(data.melted)

#Look at taxon-fractions per unique tube:
  table(data.melted$unique.tube, data.melted$Replicate)

#Look at taxon-fractions per unique tube:
  table(data.melted$unique.tube, data.melted$fraction)

#Look at taxon-fractions per unique treatment:
  table(data.melted$unique.tmt, data.melted$Replicate)
  

#Create a tube-level copies dataset by summing fractions within tubes in the Time0 samples and by doing the same for the 16O & 18O tubes (assume tube = sum of fractions):
  #Create an empty data frame:
    copies.tube <- data.frame(matrix(NA, nrow=length(levels(factor(copies$unique.tube))), ncol=7+dim(copies)[2]-12))
    names(copies.tube) <- c("tube", "Replicate", "Hour", "Isotope", "unique.tube", "unique.tmt", "copies.ul", names(copies)[13:dim(copies)[2]])

  #Calculate copies per tube:
    tubes <- unique(copies$unique.tube)
    for (k in 1:length(tubes)){
      copies.tube$unique.tube[k] <- as.character(tubes[k])
      #copies.tube$Spin[k] <- as.character(unique(copies$Spin[copies$unique.tube == tubes[k]]))
      copies.tube$Replicate[k] <- as.character(unique(copies$Replicate[copies$unique.tube == tubes[k]]))
      #copies.tube$treatment[k] <- as.character(unique(copies$treatment[copies$unique.tube == tubes[k]]))
      copies.tube$Hour[k] <- as.character(unique(copies$Hour[copies$unique.tube == tubes[k]]))
      copies.tube$Isotope[k] <- as.character(unique(copies$Isotope[copies$unique.tube == tubes[k]]))
      copies.tube$tube[k] <- as.character(unique(copies$tube[copies$unique.tube == tubes[k]]))
      copies.tube$unique.tmt[k] <- as.character(unique(copies$unique.tmt[copies$unique.tube == tubes[k]]))
      copies.tube$copies.ul[k] <- sum(copies$copies.ul[copies$unique.tube == tubes[k]])
      copies.tube[k,8:dim(copies.tube)[2]] <- apply(copies[copies$unique.tube == tubes[k], 13:dim(copies)[2]], 2, sum, na.rm=TRUE)
    }
    #Convert the factor columns to factor
      copies.tube <- as.data.frame(lapply(copies.tube, function(x) if(is.factor(x)) factor(x) else x))
      names(copies.tube)[8:dim(copies.tube)[2]] <- names(copies)[13:dim(copies)[2]]
      copies.tube[,1:20]
  
  #Check that the sum of copies.tube across all taxa is the same as the tube-level estimate of copies.tube:
    data.frame(tube.level=copies.tube$copies.ul, fraction.level=apply(copies.tube[,8:dim(copies.tube)[2]], 1, sum))

  #Melt tube-level copies data into long format by metadata:
    copies.tube.melted <- melt(data=copies.tube, id=c("tube", "Replicate", "Hour", "Isotope", "unique.tube", "unique.tmt", "copies.ul"), measure.vars=as.character(taxa.id$taxon), variable.name="taxon", value.name="t.copies.ul")

  #Merge with taxa data and reorder data frame by taxon and tube:
    copies.tube.melted <- merge(copies.tube.melted, taxa.id)
    copies.tube.melted <- copies.tube.melted[order(copies.tube.melted$taxon, copies.tube.melted$tube),]
    rownames(copies.tube.melted) <- 1:dim(copies.tube.melted)[1]     #rename observations to be sequential
    head(copies.tube.melted)
    dim(copies.tube.melted)


#Import data frame containing the treatment comparisons to perform:

  Tcompare3 <- read.table("./data/SB_TreatmentComparisons_3.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare3)
  Tcompare3

  Tcompare24 <- read.table("./data/SB_TreatmentComparisons_24.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare24)
  Tcompare24
  
  Tcompare48 <- read.table("./data/SB_TreatmentComparisons_48.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare48)
  Tcompare48
  
  
  Tcompare72 <- read.table("./data/SB_TreatmentComparisons_72.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare72)
  Tcompare72
  
  Tcompare168 <- read.table("./data/SB_TreatmentComparisons_168.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare168)
  Tcompare168
  
  TcompareTime0 <- read.table("./data/SB_TreatmentComparisons_Time0_all.txt", header=TRUE, sep="\t", colClasses=c("factor","character","character","numeric","character"))
  summary(TcompareTime0)
  TcompareTime0

#Filter the data so that only taxa with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
  #First back up the complete data frame:
    data.melted.orig <- data.melted

  #For now, we have elected to a simple, straightforward criterion that a taxon must occur in all 3 (of 3) replicates of each labeled 18O treatment AND the taxon must also occur in at least 3 (of 12) replicates of the unlabeled treatments.

  #Filtering for the different comparisons:
  #3_18O:
    Tcompare3$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
  #24_18O:
    Tcompare24$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare24$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    dev.off()
    data.melted.24 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare24$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    #48_18O:
    Tcompare48$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare48$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    dev.off()
    data.melted.48 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare48$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    #72_18O:
    Tcompare72$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare72$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    dev.off()
    data.melted.72 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare72$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
  #168_18O:
    Tcompare168$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare168$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    dev.off()
    data.melted.168 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare168$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
  #3_16O; 24_16O; 48_16O; 72_16O; 168_16O:
    Tcompare3$trt.code.1[1]
    Tcompare3$trt.refs[1]
    #Just look at the occurrence of taxa in the 16O reps from the as-yet-unfiltered data:
      explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.1[1], trt.refs=Tcompare3$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
      dev.off()
    #Filter for the 3 hour comparison:    
      explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.1[1], trt.refs=Tcompare3$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
      dev.off()
      data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.1[1], trt.refs=Tcompare3$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    #Filter for the 24 hour comparison:    
      explore.filter.taxa(DATA=data.melted.24, trt.code.1=NULL, trt.code.2=Tcompare24$trt.code.1[1], trt.refs=Tcompare24$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
      dev.off()
      data.melted.24 <- filter.taxa(DATA=data.melted.24, trt.code.1=NULL, trt.code.2=Tcompare24$trt.code.1[1], trt.refs=Tcompare24$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    #Filter for the 72 hour comparison:    
      explore.filter.taxa(DATA=data.melted.72, trt.code.1=NULL, trt.code.2=Tcompare72$trt.code.1[1], trt.refs=Tcompare72$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
      dev.off()
      data.melted.72 <- filter.taxa(DATA=data.melted.72, trt.code.1=NULL, trt.code.2=Tcompare72$trt.code.1[1], trt.refs=Tcompare72$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
    #Filter for the 168 hour comparison:    
      explore.filter.taxa(DATA=data.melted.168, trt.code.1=NULL, trt.code.2=Tcompare168$trt.code.1[1], trt.refs=Tcompare168$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)
      dev.off()
      data.melted.168 <- filter.taxa(DATA=data.melted.168, trt.code.1=NULL, trt.code.2=Tcompare168$trt.code.1[1], trt.refs=Tcompare168$trt.refs[1], vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt"), min.reps=3)

    #Renumber row.names
      row.names(data.melted.3)<- 1:dim(data.melted.3)[1]
      row.names(data.melted.24)<- 1:dim(data.melted.24)[1]
      row.names(data.melted.72)<- 1:dim(data.melted.72)[1]
      row.names(data.melted.168)<- 1:dim(data.melted.168)[1]
    #Convert the factor columns to factor:
      data.melted.3 <- as.data.frame(lapply(data.melted.3, function(x) if(is.factor(x)) factor(x) else x))
      data.melted.24 <- as.data.frame(lapply(data.melted.24, function(x) if(is.factor(x)) factor(x) else x))
      data.melted.48 <- as.data.frame(lapply(data.melted.48, function(x) if(is.factor(x)) factor(x) else x))
      data.melted.72 <- as.data.frame(lapply(data.melted.72, function(x) if(is.factor(x)) factor(x) else x))
      data.melted.168 <- as.data.frame(lapply(data.melted.168, function(x) if(is.factor(x)) factor(x) else x))



#Filter rare taxa from copies.tube.melted in the same way as was done for data.melted above:
  #Filter for the 3 hour comparison:
  copies.tube.melted.3 <- copies.tube.melted  
  copies.tube.melted.3 <- copies.tube.melted.3[copies.tube.melted.3$taxon %in% levels(data.melted.3$taxon),]
    #Renumber row.names
      row.names(copies.tube.melted.3)<- 1:dim(copies.tube.melted.3)[1]
    #Convert the factor columns to factor:
      copies.tube.melted.3 <- as.data.frame(lapply(copies.tube.melted.3, function(x) if(is.factor(x)) factor(x) else x))
      sum(levels(copies.tube.melted.3$taxon) != levels(data.melted.3$taxon))

  #Filter for the 24 hour comparison:
  copies.tube.melted.24 <- copies.tube.melted  
  copies.tube.melted.24 <- copies.tube.melted.24[copies.tube.melted.24$taxon %in% levels(data.melted.24$taxon),]
    #Renumber row.names
      row.names(copies.tube.melted.24)<- 1:dim(copies.tube.melted.24)[1]
    #Convert the factor columns to factor:
      copies.tube.melted.24 <- as.data.frame(lapply(copies.tube.melted.24, function(x) if(is.factor(x)) factor(x) else x))
      sum(levels(copies.tube.melted.24$taxon) != levels(data.melted.24$taxon))
      
  #Filter for the 48 hour comparison:
    copies.tube.melted.48 <- copies.tube.melted  
    copies.tube.melted.48 <- copies.tube.melted.48[copies.tube.melted.48$taxon %in% levels(data.melted.48$taxon),]
    #Renumber row.names
      row.names(copies.tube.melted.48)<- 1:dim(copies.tube.melted.48)[1]
    #Convert the factor columns to factor:
      copies.tube.melted.48 <- as.data.frame(lapply(copies.tube.melted.48, function(x) if(is.factor(x)) factor(x) else x))
      sum(levels(copies.tube.melted.48$taxon) != levels(data.melted.48$taxon))
      

  #Filter for the 72 hour comparison:
  copies.tube.melted.72 <- copies.tube.melted  
  copies.tube.melted.72 <- copies.tube.melted.72[copies.tube.melted.72$taxon %in% levels(data.melted.72$taxon),]
    #Renumber row.names
      row.names(copies.tube.melted.72)<- 1:dim(copies.tube.melted.72)[1]
    #Convert the factor columns to factor:
      copies.tube.melted.72 <- as.data.frame(lapply(copies.tube.melted.72, function(x) if(is.factor(x)) factor(x) else x))
      sum(levels(copies.tube.melted.72$taxon) != levels(data.melted.72$taxon))

  #Filter for the 168 hour comparison:
  copies.tube.melted.168 <- copies.tube.melted  
  copies.tube.melted.168 <- copies.tube.melted.168[copies.tube.melted.168$taxon %in% levels(data.melted.168$taxon),]
    #Renumber row.names
      row.names(copies.tube.melted.168)<- 1:dim(copies.tube.melted.168)[1]
    #Convert the factor columns to factor:
      copies.tube.melted.168 <- as.data.frame(lapply(copies.tube.melted.168, function(x) if(is.factor(x)) factor(x) else x))
      sum(levels(copies.tube.melted.168$taxon) != levels(data.melted.168$taxon))

 

#Import soil extraction data with mass of soil represented by the DNA in each tube:
  Sdat <- read.table("./data/copies_per_g_soil_50percent.txt", header=TRUE, sep="\t")
  Sdat
  summary(Sdat)


  
  
#GROSS GROWTH CALCULATIONS:

#FOR 3 HOUR (takes ~5 minutes to run on Ben's MacBook Pro):
#Calculate wad.diff, ape, r, & flux for all taxa and all comparisons - GROSS REPLICATION:
#NOTES: using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
####NOTE - when use all.taxa.calcs function r = b = gross growth. There is no death in this output
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       set prop.O.from.water = 60% based on teh Ecospher egrowth paper (for now)  This is more conservative than the 33% from E. coli data
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
# volume of fractions for 4th wedge are not equal, so I changed v.frac=30 (which is stating that all fractions are 30 ul) to 1 and normalized the qPCR data for volume before imputing (calculated absolute abundance of copies per combined SIP fraction)
  set.seed(100)
  system.time(all.comparisons.3 <- all.taxa.calcs(X.all=data.melted.3, comparisons=Tcompare3, M.soil=Sdat, vars=c("taxon", "Density", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", prop.O.from.water=0.60, v.frac=1, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_3hour_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.3)
  dim(all.comparisons.3)
  
  #####Linnea's notes: is this where we could modify and constrain all ape values to be > 0? Because there shouldn't be negative values. I ended up 
  #####changing the ape_calc function to not allow negative values for ape.
  
 

#Write the results (all.comparisons.3) to a text file:
  #dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons.3, paste(Dir.o,"SB_3hour_all_comparisons_100percent_no_0.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, paste(Dir.o,"SB_taxa_ID_100percent.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)




#FOR 24 HOUR:
#Calculate wad.diff, ape, r, & flux for all taxa and all comparisons - GROSS REPLICATION:
#NOTES: using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       set prop.O.from.water = 60% based on teh Ecospher egrowth paper (for now)  This is more conservative than the 33% from E. coli data
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.24 <- all.taxa.calcs(X.all=data.melted.24, comparisons=Tcompare24, M.soil=Sdat, vars=c("taxon", "Density", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", prop.O.from.water=0.60, v.frac=1, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o,"SB_24hour_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.24)
  dim(all.comparisons.24)


#Write the results (all.comparisons.24) to a text file:
  write.table(all.comparisons.24, paste(Dir.o,"SB_24hour_all_comparisons_100percent.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #FOR 48 HOUR:
  #Calculate wad.diff, ape, r, & flux for all taxa and all comparisons - GROSS REPLICATION:
  #NOTES: using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use linear growth model
  #       set prop.O.from.water = 60% based on teh Ecospher egrowth paper (for now)  This is more conservative than the 33% from E. coli data
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.48 <- all.taxa.calcs(X.all=data.melted.48, comparisons=Tcompare48, M.soil=Sdat, vars=c("taxon", "Density", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", prop.O.from.water=0.60, v.frac=1, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o,"SB_48hour_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.48)
  dim(all.comparisons.48)
  
  
  #Write the results (all.comparisons.48) to a text file:
  write.table(all.comparisons.48, paste(Dir.o, "SB_48hour_all_comparisons_100percent.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#FOR 72 HOUR:
#Calculate wad.diff, ape, r, & flux for all taxa and all comparisons - GROSS REPLICATION:
#NOTES: using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       set prop.O.from.water = 60% based on teh Ecospher egrowth paper (for now)  This is more conservative than the 33% from E. coli data
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.72 <- all.taxa.calcs(X.all=data.melted.72, comparisons=Tcompare72, M.soil=Sdat, vars=c("taxon", "Density", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", prop.O.from.water=0.60, v.frac=1, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_72hour_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.72)
  dim(all.comparisons.72)


#Write the results (all.comparisons.72) to a text file:
  write.table(all.comparisons.72, paste(Dir.o, "SB_72hour_all_comparisons_100percent.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)




#FOR 168 HOUR:
#Calculate wad.diff, ape, r, & flux for all taxa and all comparisons - GROSS REPLICATION:
#NOTES: using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       set prop.O.from.water = 60% based on teh Ecospher egrowth paper (for now)  This is more conservative than the 33% from E. coli data
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.168 <- all.taxa.calcs(X.all=data.melted.168, comparisons=Tcompare168, M.soil=Sdat, vars=c("taxon", "Density", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", prop.O.from.water=0.60, v.frac=1, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_168hour_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.168)
  dim(all.comparisons.168)


#Write the results (all.comparisons.168) to a text file:
  write.table(all.comparisons.168, paste(Dir.o, "SB_168hour_all_comparisons_100percent.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)




#Glue together the gross growth output from each time point:
  all.comparisons <- rbind(all.comparisons.3, all.comparisons.24, all.comparisons.48, all.comparisons.72, all.comparisons.168)
  row.names(all.comparisons)<- 1:dim(all.comparisons)[1]     #renumber row.names
  all.comparisons <- as.data.frame(lapply(all.comparisons, function(x) if(is.factor(x)) factor(x) else x))     #convert the factor columns to factor:

#Write the combined results (all.comparisons) to a text file:
  write.table(all.comparisons, paste(Dir.o, "SB_all_comparisons_100percent.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  




#NET GROWTH CALCULATIONS:
    
  #FOR 3 HOUR (takes ~5 minutes to run on Ben's MacBook Pro):
    #Calculate population-based r & flux for all taxa - NET REPLICATION:
    #NOTES: the treatments for time t abundances are done in 3 ways: 16O, 18O, and both 16O & 18O treatment replicates combined ('Time0' is the only treatment used for time 0)
    #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
    #       using system.time calculates how long it took the function to run
    #       use linear growth model
    #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  # all.taxa.calcs.pop output r = net growth
  #### note that v=c(30, 30) means volume of fraction for t1 and t2 fractions with the assumption that fractions volumes are equal within a tube. I changed this to vol=c(1, 1) because my volumes for sip fractions differed, so I used total copy numbers in SIP fractions and volume of 1 instead
set.seed(100)
system.time(all.comparisons.pop.3 <- all.taxa.calcs.pop(X.all=copies.tube.melted.3, comparisons=TcompareTime0[1:3,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_3hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.3)
  dim(all.comparisons.pop.3)


#Write the results (all.comparisons.pop.3) to a text file:
  write.table(all.comparisons.pop.3, paste(Dir.o, "SB_3hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)




#FOR 24 HOUR:
#Calculate population-based r & flux for all taxa - NET REPLICATION:
#NOTES: the treatments for time t abundances are done in 3 ways: 16O, 18O, and both 16O & 18O treatment replicates combined ('Time0' is the only treatment used for time 0)
#       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.24 <- all.taxa.calcs.pop(X.all=copies.tube.melted.24, comparisons=TcompareTime0[4:6,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_24hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.24)
  dim(all.comparisons.pop.24)


#Write the results (all.comparisons.pop.24) to a text file:
  write.table(all.comparisons.pop.24, paste(Dir.o, "SB_24hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #FOR 48 HOUR:
  #Calculate population-based r & flux for all taxa - NET REPLICATION:
  #NOTES: the treatments for time t abundances are done in 3 ways: 16O, 18O, and both 16O & 18O treatment replicates combined ('Time0' is the only treatment used for time 0)
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use linear growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.48 <- all.taxa.calcs.pop(X.all=copies.tube.melted.48, comparisons=TcompareTime0[7:9,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_48hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.48)
  dim(all.comparisons.pop.48)
  
  
  #Write the results (all.comparisons.pop.48) to a text file:
  write.table(all.comparisons.pop.48, paste(Dir.o, "SB_48hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#FOR 72 HOUR:
#Calculate population-based r & flux for all taxa - NET REPLICATION:
#NOTES: the treatments for time t abundances are done in 3 ways: 16O, 18O, and both 16O & 18O treatment replicates combined ('Time0' is the only treatment used for time 0)
#       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.72 <- all.taxa.calcs.pop(X.all=copies.tube.melted.72, comparisons=TcompareTime0[10:12,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_72hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.72)
  dim(all.comparisons.pop.72)


#Write the results (all.comparisons.pop.72) to a text file:
  write.table(all.comparisons.pop.72, paste(Dir.o, "SB_72hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)




#FOR 168 HOUR:
#Calculate population-based r & flux for all taxa - NET REPLICATION:
#NOTES: the treatments for time t abundances are done in 3 ways: 16O, 18O, and both 16O & 18O treatment replicates combined ('Time0' is the only treatment used for time 0)
#       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
#       using system.time calculates how long it took the function to run
#       use linear growth model
#       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.168 <- all.taxa.calcs.pop(X.all=copies.tube.melted.168, comparisons=TcompareTime0[13:15,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_168hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.168)
  dim(all.comparisons.pop.168)


#Write the results (all.comparisons.pop.168) to a text file:
  write.table(all.comparisons.pop.168, paste(Dir.o, "SB_168hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


  #FOR 3-24 HOUR INTERVAL:
  #Calculate population-based r & flux for all taxa - NET REPLICATION:
  #NOTES: the treatments for time 0 and time t abundances are done in only one way: 3_18O and 24_18O
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use linear growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.3.24 <- all.taxa.calcs.pop(X.all=copies.tube.melted.3.24, comparisons=TcompareTime0[16,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_3-24hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.3.24)
  dim(all.comparisons.pop.3.24)
  
  
  #Write the results (all.comparisons.pop.3.24) to a text file:
  write.table(all.comparisons.pop.3.24, paste(Dir.o, "SB_3-24hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  
  
  
  
  #FOR 24-48 HOUR INTERVAL:
  #Calculate population-based r & flux for all taxa - NET REPLICATION:
  #NOTES: the treatments for time 0 and time t abundances are done in only one way: 24_18O and 48_18O
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use linear growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.24.48 <- all.taxa.calcs.pop(X.all=copies.tube.melted.24.48, comparisons=TcompareTime0[17,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_24-48hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.24.48)
  dim(all.comparisons.pop.24.48)
  
  
  #Write the results (all.comparisons.pop.24.48) to a text file:
  write.table(all.comparisons.pop.24.48, paste(Dir.o, "SB_24-48hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  
  #FOR 48-72 HOUR INTERVAL:
  #Calculate population-based r & flux for all taxa - NET REPLICATION:
  #NOTES: the treatments for time 0 and time t abundances are done in only one way: 48_18O and 72_18O
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use linear growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.48.72 <- all.taxa.calcs.pop(X.all=copies.tube.melted.48.72, comparisons=TcompareTime0[18,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_48-72hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.48.72)
  dim(all.comparisons.pop.48.72)
  
  
  #Write the results (all.comparisons.pop.48.72) to a text file:
  write.table(all.comparisons.pop.48.72, paste(Dir.o, "SB_48-72hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  
  
  
  #FOR 72-168 HOUR INTERVAL:
  #Calculate population-based r & flux for all taxa - NET REPLICATION:
  #NOTES: the treatments for time 0 and time t abundances are done in only one way: 72_18O and 168_18O
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use linear growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.pop.72.168 <- all.taxa.calcs.pop(X.all=copies.tube.melted.72.168, comparisons=TcompareTime0[19,], M.soil=Sdat, vars=c("taxon", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"), growth.model="linear", vol=c(1, 1), copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000))
  bootstrapped.filenames <- paste(Dir.o, c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  new.bootstrapped.filenames <- paste(Dir.o, "SB_72-168hour_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.pop.72.168)
  dim(all.comparisons.pop.72.168)
  
  
  #Write the results (all.comparisons.pop.72.168) to a text file:
  write.table(all.comparisons.pop.72.168, paste(Dir.o, "SB_72-168hour_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  
  
  
  
  #Glue together the net growth output from each time interval:
  all.comparisons.pop <- rbind(all.comparisons.pop.3, all.comparisons.pop.24, all.comparisons.pop.72, all.comparisons.pop.168, all.comparisons.pop.3.24, all.comparisons.pop.24.48, all.comparisons.pop.48.72, all.comparisons.pop.72.168)
  row.names(all.comparisons.pop)<- 1:dim(all.comparisons.pop)[1]     #renumber row.names
  all.comparisons.pop <- as.data.frame(lapply(all.comparisons.pop, function(x) if(is.factor(x)) factor(x) else x))     #convert the factor columns to factor:
  
  #Write the combined results (all.comparisons.pop) to a text file:
  write.table(all.comparisons.pop, paste(Dir.o, "SB_all_comparisons_pop.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  
  #Save workspace:
  #save(list=ls(all.names=TRUE), file="./qSIP_workspaces/100_SB_01_.RData", envir=.GlobalEnv)
  
  
  
  

