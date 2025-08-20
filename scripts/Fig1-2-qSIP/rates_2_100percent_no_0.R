#Reload the saved workspace resulting from the previous script (above):
Dir.o = "./output/100-growth-rates/"
#load("./qSIP_workspaces/100_SB_01.RData")
      
    
      ##########  Calculate & verify birth, death, net growth, and flux rates & ape & N:
      
      #Read in the necessary bootstrapped growth, C flux, and abundance estimates from the all.taxa.calcs & all.taxa.calcs.pop output:
      r.gross.boots.3 <- read.table(paste(Dir.o, "SB_3hour_bootstrapped_r.txt", sep=""), header=TRUE, sep="")
      f.gross.boots.3 <- read.table(paste(Dir.o, "SB_3hour_bootstrapped_C_fluxes.txt", sep=""), header=TRUE, sep="")
      r.net.boots.3 <- read.table(paste(Dir.o, "SB_3hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.3 <- read.table(paste(Dir.o, "SB_3hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.3 <- read.table(paste(Dir.o, "SB_3hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.3 <- read.table(paste(Dir.o, "SB_3hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.gross.boots.24 <- read.table(paste(Dir.o, "SB_24hour_bootstrapped_r.txt", sep=""), header=TRUE, sep="")
      f.gross.boots.24 <- read.table(paste(Dir.o, "SB_24hour_bootstrapped_C_fluxes.txt", sep=""), header=TRUE, sep="")
      r.net.boots.24 <- read.table(paste(Dir.o, "SB_24hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.24 <- read.table(paste(Dir.o, "SB_24hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.24 <- read.table(paste(Dir.o, "SB_24hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.24 <- read.table(paste(Dir.o, "SB_24hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.gross.boots.48 <- read.table(paste(Dir.o, "SB_48hour_bootstrapped_r.txt", sep=""), header=TRUE, sep="")
      f.gross.boots.48 <- read.table(paste(Dir.o, "SB_48hour_bootstrapped_C_fluxes.txt", sep=""), header=TRUE, sep="")
      r.net.boots.48 <- read.table(paste(Dir.o, "SB_48hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.48 <- read.table(paste(Dir.o, "SB_48hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.48 <- read.table(paste(Dir.o, "SB_48hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.48 <- read.table(paste(Dir.o, "SB_48hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.gross.boots.72 <- read.table(paste(Dir.o, "SB_72hour_bootstrapped_r.txt", sep=""), header=TRUE, sep="")
      f.gross.boots.72 <- read.table(paste(Dir.o, "SB_72hour_bootstrapped_C_fluxes.txt", sep=""), header=TRUE, sep="")
      r.net.boots.72 <- read.table(paste(Dir.o, "SB_72hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.72 <- read.table(paste(Dir.o, "SB_72hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.72 <- read.table(paste(Dir.o, "SB_72hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.72 <- read.table(paste(Dir.o, "SB_72hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.gross.boots.168 <- read.table(paste(Dir.o, "SB_168hour_bootstrapped_r.txt", sep=""), header=TRUE, sep="")
      f.gross.boots.168 <- read.table(paste(Dir.o, "SB_168hour_bootstrapped_C_fluxes.txt", sep=""), header=TRUE, sep="")
      r.net.boots.168 <- read.table(paste(Dir.o, "SB_168hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.168 <- read.table(paste(Dir.o, "SB_168hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.168 <- read.table(paste(Dir.o, "SB_168hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.168 <- read.table(paste(Dir.o, "SB_168hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.net.boots.3.24 <- read.table(paste(Dir.o, "SB_3-24hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.3.24 <- read.table(paste(Dir.o, "SB_3-24hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.3.24 <- read.table(paste(Dir.o, "SB_3-24hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.3.24 <- read.table(paste(Dir.o, "SB_3-24hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.net.boots.24.48 <- read.table(paste(Dir.o, "SB_24-48hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.24.48 <- read.table(paste(Dir.o, "SB_24-48hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.24.48 <- read.table(paste(Dir.o, "SB_24-48hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.24.48 <- read.table(paste(Dir.o, "SB_24-48hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.net.boots.48.72 <- read.table(paste(Dir.o, "SB_48-72hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.48.72 <- read.table(paste(Dir.o, "SB_48-72hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.48.72 <- read.table(paste(Dir.o, "SB_48-72hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.48.72 <- read.table(paste(Dir.o, "SB_48-72hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      r.net.boots.72.168 <- read.table(paste(Dir.o, "SB_72-168hour_bootstrapped_r_pop.txt", sep=""), header=TRUE, sep="")
      f.net.boots.72.168 <- read.table(paste(Dir.o, "SB_72-168hour_bootstrapped_C_fluxes_pop.txt", sep=""), header=TRUE, sep="")
      N.Tt.boots.72.168 <- read.table(paste(Dir.o, "SB_72-168hour_bootstrapped_N_Tt_pop.txt", sep=""), header=TRUE, sep="")
      N.T0.boots.72.168 <- read.table(paste(Dir.o, "SB_72-168hour_bootstrapped_N_T0_pop.txt", sep=""), header=TRUE, sep="")
      
      #Note that they are not all of the same dimension:
      dim(r.gross.boots.3)
      dim(f.gross.boots.3)
      dim(r.net.boots.3)
      dim(f.net.boots.3)
      dim(N.Tt.boots.3)
      dim(N.T0.boots.3)
      
      dim(r.gross.boots.24)
      dim(f.gross.boots.24)
      dim(r.net.boots.24)
      dim(f.net.boots.24)
      dim(N.Tt.boots.24)
      dim(N.T0.boots.24)
      
      dim(r.gross.boots.48)
      dim(f.gross.boots.48)
      dim(r.net.boots.48)
      dim(f.net.boots.48)
      dim(N.Tt.boots.48)
      dim(N.T0.boots.48)
      
      dim(r.gross.boots.72)
      dim(f.gross.boots.72)
      dim(r.net.boots.72)
      dim(f.net.boots.72)
      dim(N.Tt.boots.72)
      dim(N.T0.boots.72)
      
      dim(r.gross.boots.168)
      dim(f.gross.boots.168)
      dim(r.net.boots.168)
      dim(f.net.boots.168)
      dim(N.Tt.boots.168)
      dim(N.T0.boots.168)
      
      dim(r.net.boots.3.24)
      dim(f.net.boots.3.24)
      dim(N.Tt.boots.3.24)
      dim(N.T0.boots.3.24)
      
      dim(r.net.boots.24.48)
      dim(f.net.boots.24.48)
      dim(N.Tt.boots.24.48)
      dim(N.T0.boots.24.48)
      
      dim(r.net.boots.48.72)
      dim(f.net.boots.48.72)
      dim(N.Tt.boots.48.72)
      dim(N.T0.boots.48.72)
      
      dim(r.net.boots.72.168)
      dim(f.net.boots.72.168)
      dim(N.Tt.boots.72.168)
      dim(N.T0.boots.72.168)
      
      #Get the complete set of taxa for which to build the master 'all.rates' data.frames:
      taxa.set.3 <- union(levels(factor(as.character(r.gross.boots.3$taxonID))), levels(factor(as.character(r.net.boots.3$taxonID))))
      taxa.set.24 <- union(levels(factor(as.character(r.gross.boots.24$taxonID))), levels(factor(as.character(r.net.boots.24$taxonID))))
      taxa.set.48 <- union(levels(factor(as.character(r.gross.boots.48$taxonID))), levels(factor(as.character(r.net.boots.48$taxonID))))
      taxa.set.72 <- union(levels(factor(as.character(r.gross.boots.72$taxonID))), levels(factor(as.character(r.net.boots.72$taxonID))))
      taxa.set.168 <- union(levels(factor(as.character(r.gross.boots.168$taxonID))), levels(factor(as.character(r.net.boots.168$taxonID))))
      taxa.set.serials <- union(union(taxa.set.24, union(taxa.set.72, taxa.set.168)), taxa.set.48)
      taxa.set.intervals <- union(union(levels(factor(as.character(r.net.boots.24.48$taxonID))), levels(factor(as.character(r.net.boots.48.72$taxonID)))), levels(factor(as.character(r.net.boots.72.168$taxonID))))
      taxa.set.all <- union(taxa.set.serials, taxa.set.intervals)
      taxa.set.all.order <- order(as.numeric(gsub(pattern="^sip(\\d+)$", replacement="\\1", x=taxa.set.all, perl=TRUE)))    #sort by taxonID
      taxa.set.all <- taxa.set.all[taxa.set.all.order]
      taxa.set.all
      
      
      #Create 15 data.frames to store all combined results (Note: each data.frame corresponds to one row (comparison) in TcompareTime0):
      
      #Create a quick & dirty function to do this quickly for all 12 'serially-based' data.frames (i.e., the ones that go from Time 0 to Time t):
      create.all.rates.func <- function(qSIP.data, pop.data, Time.t.group.code, set.of.all.taxa){
        #First write all the taxa:
        all.rates <- data.frame(taxonID=set.of.all.taxa)
        #Second, include the qSIP-derived 'gross' rates:
        qSIP.data.subset <- qSIP.data[, names(qSIP.data) %in% c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", "ape.obs", "ape.boot.mean", "ape.boot.median", "ape.boot.CI.L", "ape.boot.CI.U", "r.obs", "r.boot.mean", "r.boot.median", "r.boot.CI.L", "r.boot.CI.U", "f.obs", "f.boot.mean", "f.boot.median", "f.boot.CI.L", "f.boot.CI.U", "message")]
        all.rates <- merge(all.rates, qSIP.data.subset, by="taxonID", all.x=TRUE, all.y=FALSE)
        names(all.rates)[10:14] <- c("b.obs", "b.boot.mean", "b.boot.median", "b.boot.CI.L", "b.boot.CI.U")
        #Third, include the population-derived 'net' rates from the relevant rows in the 'all.comparisons.pop' data.frame:
        pop.data.subset <- pop.data[pop.data$group.code.t == Time.t.group.code,]
        row.names(pop.data.subset) <- 1:dim(pop.data.subset)[1]     #renumber row.names
        pop.data.subset <- as.data.frame(lapply(pop.data.subset, function(x) if(is.factor(x)) factor(x) else x))     #convert the factor columns to factor
        names(pop.data.subset)[c(2,25)] <- c("comparisonID.pop", "message.pop")
        all.rates <- merge(all.rates, pop.data.subset, by="taxonID", all.x=TRUE, all.y=FALSE)
        names(all.rates)[15:19] <- c("f.gross.obs", "f.gross.boot.mean", "f.gross.boot.median", "f.gross.boot.CI.L", "f.gross.boot.CI.U")
        names(all.rates)[24:33] <- c("r.net.obs", "r.net.boot.mean", "r.net.boot.median", "r.net.boot.CI.L", "r.net.boot.CI.U", "f.net.obs", "f.net.boot.mean", "f.net.boot.median", "f.net.boot.CI.L", "f.net.boot.CI.U")
        names(all.rates)[34:43] <- c("N.tot.T0.obs.mean", "N.tot.T0.boot.mean", "N.tot.T0.boot.median", "N.tot.T0.boot.CI.L", "N.tot.T0.boot.CI.U", "N.tot.Tt.obs.mean", "N.tot.Tt.boot.mean", "N.tot.Tt.boot.median", "N.tot.Tt.boot.CI.L", "N.tot.Tt.boot.CI.U")
        all.rates
      }
      
      #Create an additional function to add in the death rates, C fluxes to dead cells (turnover), and abundance of light copies at time T:
      add.death.rates.func <- function(all.rates, r.gross.boots, r.net.boots, f.gross.boots, f.net.boots, N.T0.boots, N.Tt.boots, growth.model, days, CI=0.95){
        Time.t.group.code <- levels(factor(as.character(all.rates$group.code.t)))
        #Calculate death rate as d = r - b (convention is that death rates are negative):
        #Observed:
        all.rates$d.obs <- all.rates$r.net.obs - all.rates$b.obs
        #Bootstrapped:
        taxa.only <- data.frame(taxonID=all.rates$taxonID)
        r.net.boots.all <- merge(taxa.only, r.net.boots[r.net.boots$group.code.t == Time.t.group.code, ], by="taxonID", all.x=TRUE, all.y=FALSE)
        r.net.boots.all <- r.net.boots.all[order(r.net.boots.all$taxonID), ]
        r.net.boots.only <- r.net.boots.all[, 5:dim(r.net.boots.all)[2]]
        r.gross.boots.all <- merge(taxa.only, r.gross.boots, by="taxonID", all.x=TRUE, all.y=FALSE)
        r.gross.boots.all <- r.gross.boots.all[order(r.gross.boots.all$taxonID), ]
        b.boots.only <- r.gross.boots.all[, 5:dim(r.gross.boots.all)[2]]
        d.boots.only <- r.net.boots.only - b.boots.only
        d.boot.mean <- as.numeric(apply(d.boots.only, 1, mean, na.rm=TRUE))
        d.boot.median <- as.numeric(apply(d.boots.only, 1, median, na.rm=TRUE))
        d.boot.CI.L <- as.numeric(apply(d.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        d.boot.CI.U <- as.numeric(apply(d.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
        death.rates.to.merge <- data.frame(taxonID=r.gross.boots.all$taxonID, d.boot.mean, d.boot.median, d.boot.CI.L, d.boot.CI.U)
        all.rates <- merge(all.rates, death.rates.to.merge, by="taxonID", all.x=TRUE, all.y=FALSE)
        #Calculate C flux into dead cells (turnover) as f.death = f.net - f.gross (convention is that the flux into dead cells is negative):
        #Observed:
        all.rates$f.death.obs <- all.rates$f.net.obs - all.rates$f.gross.obs
        #Bootstrapped:
        f.net.boots.all <- merge(taxa.only, f.net.boots[f.net.boots$group.code.t == Time.t.group.code, ], by="taxonID", all.x=TRUE, all.y=FALSE)
        f.net.boots.all <- f.net.boots.all[order(f.net.boots.all$taxonID), ]
        f.net.boots.only <- f.net.boots.all[, 5:dim(f.net.boots.all)[2]]
        f.gross.boots.all <- merge(taxa.only, f.gross.boots, by="taxonID", all.x=TRUE, all.y=FALSE)
        f.gross.boots.all <- f.gross.boots.all[order(f.gross.boots.all$taxonID), ]
        f.gross.boots.only <- f.gross.boots.all[, 5:dim(f.gross.boots.all)[2]]
        f.death.boots.only <- f.net.boots.only - f.gross.boots.only
        f.death.boot.mean <- as.numeric(apply(f.death.boots.only, 1, mean, na.rm=TRUE))
        f.death.boot.median <- as.numeric(apply(f.death.boots.only, 1, median, na.rm=TRUE))
        f.death.boot.CI.L <- as.numeric(apply(f.death.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        f.death.boot.CI.U <- as.numeric(apply(f.death.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
        death.fluxes.to.merge <- data.frame(taxonID=f.gross.boots.all$taxonID, f.death.boot.mean, f.death.boot.median, f.death.boot.CI.L, f.death.boot.CI.U)
        all.rates <- merge(all.rates, death.fluxes.to.merge, by="taxonID", all.x=TRUE, all.y=FALSE)
        #Calculate the abundance of light copies at time T (light copies/g soil) according to the growth model used:
        N.Tt.boots.all <- merge(taxa.only, N.Tt.boots[N.Tt.boots$group.code.t == Time.t.group.code, ], by="taxonID", all.x=TRUE, all.y=FALSE)
        N.Tt.boots.all <- N.Tt.boots.all[order(N.Tt.boots.all$taxonID), ]
        N.Tt.boots.only <- N.Tt.boots.all[, 5:dim(N.Tt.boots.all)[2]]
        N.T0.boots.all <- merge(taxa.only, N.T0.boots[N.T0.boots$group.code.t == Time.t.group.code, ], by="taxonID", all.x=TRUE, all.y=FALSE)
        N.T0.boots.all <- N.T0.boots.all[order(N.T0.boots.all$taxonID), ]
        N.T0.boots.only <- N.T0.boots.all[, 5:dim(N.T0.boots.all)[2]]
        if (growth.model == "exponential"){
          all.rates$N.light.Tt.obs <- all.rates$N.tot.Tt.obs.mean / (exp(all.rates$b.obs*days))
          N.light.Tt.boots.only <- N.Tt.boots.only / (exp(b.boots.only*days))
        }
        if (growth.model == "linear"){
          all.rates$N.light.Tt.obs <- all.rates$N.tot.Tt.obs.mean - (all.rates$b.obs*days)
          N.light.Tt.boots.only <- N.Tt.boots.only - (b.boots.only*days)
        }
        N.light.Tt.boot.mean <- as.numeric(apply(N.light.Tt.boots.only, 1, mean, na.rm=TRUE))
        N.light.Tt.boot.median <- as.numeric(apply(N.light.Tt.boots.only, 1, median, na.rm=TRUE))
        N.light.Tt.boot.CI.L <- as.numeric(apply(N.light.Tt.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        N.light.Tt.boot.CI.U <- as.numeric(apply(N.light.Tt.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
        N.light.Tt.to.merge <- data.frame(taxonID=N.Tt.boots.all$taxonID, N.light.Tt.boot.mean, N.light.Tt.boot.median, N.light.Tt.boot.CI.L, N.light.Tt.boot.CI.U)
        all.rates <- merge(all.rates, N.light.Tt.to.merge, by="taxonID", all.x=TRUE, all.y=FALSE)
        taxa.order <- order(as.numeric(gsub(pattern="^sip(\\d+)$", replacement="\\1", x=all.rates$taxonID, perl=TRUE)))    #sort by taxonID
        all.rates <- all.rates[taxa.order, ]
        all.rates      
      }
      
      #Create another quick & dirty function to create a data.frame for storing all combined results quickly for the 3 'interval-based' data.frames:
      create.all.rates.interval.func <- function(pop.data, set.of.all.taxa){
        #First write all the taxa:
        all.rates <- data.frame(taxonID=set.of.all.taxa)
        #Second, include columns for the qSIP-derived 'gross' rates:
        names.to.add <- c("comparisonID", "trt.code.1", "trt.code.2", "ape.obs", "ape.boot.mean", "ape.boot.median", "ape.boot.CI.L", "ape.boot.CI.U", "b.obs", "b.boot.mean", "b.boot.median", "b.boot.CI.L", "b.boot.CI.U", "f.gross.obs", "f.gross.boot.mean", "f.gross.boot.median", "f.gross.boot.CI.L", "f.gross.boot.CI.U", "message")
        all.rates <- data.frame(all.rates, data.frame(matrix(NA, nrow=dim(all.rates)[1], ncol=length(names.to.add))))
        names(all.rates)[2:dim(all.rates)[2]] <- names.to.add
        #Third, include the population-derived 'net' rates from the 'all.comparisons.pop' data.frame:
        pop.data.subset <- pop.data
        names(pop.data.subset)[c(2,25)] <- c("comparisonID.pop", "message.pop")
        all.rates <- merge(all.rates, pop.data.subset, by="taxonID", all.x=TRUE, all.y=FALSE)
        names(all.rates)[15:19] <- c("f.gross.obs", "f.gross.boot.mean", "f.gross.boot.median", "f.gross.boot.CI.L", "f.gross.boot.CI.U")
        names(all.rates)[24:33] <- c("r.net.obs", "r.net.boot.mean", "r.net.boot.median", "r.net.boot.CI.L", "r.net.boot.CI.U", "f.net.obs", "f.net.boot.mean", "f.net.boot.median", "f.net.boot.CI.L", "f.net.boot.CI.U")
        names(all.rates)[34:43] <- c("N.tot.T0.obs.mean", "N.tot.T0.boot.mean", "N.tot.T0.boot.median", "N.tot.T0.boot.CI.L", "N.tot.T0.boot.CI.U", "N.tot.Tt.obs.mean", "N.tot.Tt.boot.mean", "N.tot.Tt.boot.median", "N.tot.Tt.boot.CI.L", "N.tot.Tt.boot.CI.U")
        all.rates
      }
      
      #Create an additional function to add in the birth rates & C fluxes, and also create r.gross.boots & f.gross.boots for the 'interval-based' data.frames:
      add.birth.flux.interval.func <- function(all.rates, T1.all.rates, T2.all.rates, T1.r.gross.boots, T2.r.gross.boots, T1.N.Tt.boots, T2.N.Tt.boots, T2.f.gross.boots, growth.model, interval.days, T1.days, T2.days, copies.cell=6, pgC.cell=0.1, CI=0.95){
        T1.all.rates <- T1.all.rates[order(T1.all.rates$taxonID), ]
        T2.all.rates <- T2.all.rates[order(T2.all.rates$taxonID), ]
        #Important: values for 'copies.cell' & 'pgC.cell' must be the same as those used in all.taxa.calcs!
        all.rates <- all.rates[order(all.rates$taxonID), ]
        #Calculate birth rate for the interval (T1 to T2):
        T1.group.code <- levels(factor(as.character(all.rates$group.code.0)))
        T2.group.code <- levels(factor(as.character(all.rates$group.code.t)))
        taxa.only <- data.frame(taxonID=all.rates$taxonID)
        T1.b.boots.all <- merge(taxa.only, T1.r.gross.boots, by="taxonID", all.x=TRUE, all.y=FALSE)
        T1.b.boots.all <- T1.b.boots.all[order(T1.b.boots.all$taxonID), ]
        T2.b.boots.all <- merge(taxa.only, T2.r.gross.boots, by="taxonID", all.x=TRUE, all.y=FALSE)
        T2.b.boots.all <- T2.b.boots.all[order(T2.b.boots.all$taxonID), ]
        T1.N.Tt.boots.all <- merge(taxa.only, T1.N.Tt.boots[T1.N.Tt.boots$group.code.t == T1.group.code, ], by="taxonID", all.x=TRUE, all.y=FALSE)
        T1.N.Tt.boots.all <- T1.N.Tt.boots.all[order(T1.N.Tt.boots.all$taxonID), ]
        T2.N.Tt.boots.all <- merge(taxa.only, T2.N.Tt.boots[T2.N.Tt.boots$group.code.t == T2.group.code, ], by="taxonID", all.x=TRUE, all.y=FALSE)
        T2.N.Tt.boots.all <- T2.N.Tt.boots.all[order(T2.N.Tt.boots.all$taxonID), ]
        T1.b.boots.only <- T1.b.boots.all[, 5:dim(T1.b.boots.all)[2]]
        T2.b.boots.only <- T2.b.boots.all[, 5:dim(T2.b.boots.all)[2]]
        T1.N.Tt.boots.only <- T1.N.Tt.boots.all[, 5:dim(T1.N.Tt.boots.all)[2]]
        T2.N.Tt.boots.only <- T2.N.Tt.boots.all[, 5:dim(T2.N.Tt.boots.all)[2]]
        if (growth.model == "exponential"){
          all.rates$b.obs <- (1/interval.days)*log((all.rates$N.tot.Tt.obs.mean*T1.all.rates$N.light.Tt.obs)/(all.rates$N.tot.T0.obs.mean*T2.all.rates$N.light.Tt.obs))
          T1.N.light.Tt.boots.only <- T1.N.Tt.boots.only / (exp(T1.b.boots.only*T1.days))
          T2.N.light.Tt.boots.only <- T2.N.Tt.boots.only / (exp(T2.b.boots.only*T2.days))
          interval.b.boots <- (1/interval.days)*log((T2.N.Tt.boots.only*T1.N.light.Tt.boots.only)/(T1.N.Tt.boots.only*T2.N.light.Tt.boots.only))
        }
        if (growth.model == "linear"){
          all.rates$b.obs <- (1/interval.days)*(all.rates$N.tot.Tt.obs.mean - all.rates$N.tot.T0.obs.mean - T2.all.rates$N.light.Tt.obs + T1.all.rates$N.light.Tt.obs)
          T1.N.light.Tt.boots.only <- T1.N.Tt.boots.only - (T1.b.boots.only*T1.days)
          T2.N.light.Tt.boots.only <- T2.N.Tt.boots.only - (T2.b.boots.only*T2.days)
          interval.b.boots <- (1/interval.days)*(T2.N.Tt.boots.only - T1.N.Tt.boots.only - T2.N.light.Tt.boots.only + T1.N.light.Tt.boots.only)
        }
        names(interval.b.boots) <- names(T1.b.boots.only)
        all.rates$b.boot.mean <- as.numeric(apply(interval.b.boots, 1, mean, na.rm=TRUE))
        all.rates$b.boot.median <- as.numeric(apply(interval.b.boots, 1, median, na.rm=TRUE))
        all.rates$b.boot.CI.L <- as.numeric(apply(interval.b.boots, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        all.rates$b.boot.CI.U <- as.numeric(apply(interval.b.boots, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
        #Calculate C flux into all cells (production) for the interval (T1 to T2):
        copies.soil.observed <- T2.all.rates$f.gross.obs * (1/T2.all.rates$b.obs) * copies.cell * (1/pgC.cell)
        T2.f.gross.boots.all <- merge(taxa.only, T2.f.gross.boots, by="taxonID", all.x=TRUE, all.y=FALSE)
        T2.f.gross.boots.all <- T2.f.gross.boots.all[order(T2.f.gross.boots.all$taxonID), ]
        T2.f.gross.boots.only <- T2.f.gross.boots.all[, 5:dim(T2.f.gross.boots.all)[2]]
        copies.soil.boots <- T2.f.gross.boots.only * (1/T2.b.boots.only) * copies.cell * (1/pgC.cell)
        if (growth.model == "exponential"){
          #(uses same calculation approach as for all.taxa.calcs, where abundance from T2 is used to convert the per-day growth rate to a flux)
          all.rates$f.gross.obs <- all.rates$b.obs * copies.soil.observed * (1/copies.cell) * pgC.cell
          interval.f.gross.boots <- interval.b.boots * copies.soil.boots * (1/copies.cell) * pgC.cell
        }
        if (growth.model == "linear"){
          all.rates$f.gross.obs <- all.rates$b.obs * (1/copies.cell) * pgC.cell
          interval.f.gross.boots <- interval.b.boots * (1/copies.cell) * pgC.cell
        }
        names(interval.f.gross.boots) <- paste("f", 1:1000, sep="")
        all.rates$f.gross.boot.mean <- as.numeric(apply(interval.f.gross.boots, 1, mean, na.rm=TRUE))
        all.rates$f.gross.boot.median <- as.numeric(apply(interval.f.gross.boots, 1, median, na.rm=TRUE))
        all.rates$f.gross.boot.CI.L <- as.numeric(apply(interval.f.gross.boots, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        all.rates$f.gross.boot.CI.U <- as.numeric(apply(interval.f.gross.boots, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
        prefix <- data.frame(matrix(NA, nrow=dim(interval.b.boots)[1], ncol=3))
        names(prefix) <- names(T1.r.gross.boots)[2:4]
        r.gross.boots.new <- data.frame(taxa.only, prefix, interval.b.boots)
        r.gross.boots.new <- r.gross.boots.new[order(r.gross.boots.new$taxonID), ]
        f.gross.boots.new <- data.frame(taxa.only, prefix, interval.f.gross.boots)
        f.gross.boots.new <- f.gross.boots.new[order(f.gross.boots.new$taxonID), ]
        list(all.rates=all.rates, r.gross.boots=r.gross.boots.new, f.gross.boots=f.gross.boots.new)
      }
      
      
      
#################################################################### end other code     
    #For 0-3 hours, using only 16O reps to calculate T0 & Tt abundances & net rates:
      all.rates.3h.16 <- create.all.rates.func(qSIP.data=all.comparisons.3, pop.data=all.comparisons.pop.3, Time.t.group.code="t3_16O", set.of.all.taxa=taxa.set.all)
      all.rates.3h.16 <- add.death.rates.func(all.rates=all.rates.3h.16, r.gross.boots=r.gross.boots.3, r.net.boots=r.net.boots.3, f.gross.boots=f.gross.boots.3, f.net.boots=f.net.boots.3, N.T0.boots=N.T0.boots.3, N.Tt.boots=N.Tt.boots.3, growth.model="linear", days=TcompareTime0$days[1], CI=0.95)

    #For 0-3 hours, using 16O reps to calculate T0 abundance & 18O reps to calculate Tt abundance & net rates:
      all.rates.3h.18 <- create.all.rates.func(qSIP.data=all.comparisons.3, pop.data=all.comparisons.pop.3, Time.t.group.code="t3_18O", set.of.all.taxa=taxa.set.all)
      all.rates.3h.18 <- add.death.rates.func(all.rates=all.rates.3h.18, r.gross.boots=r.gross.boots.3, r.net.boots=r.net.boots.3, f.gross.boots=f.gross.boots.3, f.net.boots=f.net.boots.3, N.T0.boots=N.T0.boots.3, N.Tt.boots=N.Tt.boots.3, growth.model="linear", days=TcompareTime0$days[2], CI=0.95)

    #For 0-3 hours, using 16O reps to calculate T0 abundance & 16O+18O reps to calculate Tt abundance & net rates:
      all.rates.3h.16.18 <- create.all.rates.func(qSIP.data=all.comparisons.3, pop.data=all.comparisons.pop.3, Time.t.group.code="t3_16O.t3_18O", set.of.all.taxa=taxa.set.all)
      all.rates.3h.16.18 <- add.death.rates.func(all.rates=all.rates.3h.16.18, r.gross.boots=r.gross.boots.3, r.net.boots=r.net.boots.3, f.gross.boots=f.gross.boots.3, f.net.boots=f.net.boots.3, N.T0.boots=N.T0.boots.3, N.Tt.boots=N.Tt.boots.3, growth.model="linear", days=TcompareTime0$days[3], CI=0.95)

    #For 0-24 hours, using only 16O reps to calculate T0 & Tt abundances & net rates:
      all.rates.24h.16 <- create.all.rates.func(qSIP.data=all.comparisons.24, pop.data=all.comparisons.pop.24, Time.t.group.code="t24_16O", set.of.all.taxa=taxa.set.all)
      all.rates.24h.16 <- add.death.rates.func(all.rates=all.rates.24h.16, r.gross.boots=r.gross.boots.24, r.net.boots=r.net.boots.24, f.gross.boots=f.gross.boots.24, f.net.boots=f.net.boots.24, N.T0.boots=N.T0.boots.24, N.Tt.boots=N.Tt.boots.24, growth.model="linear", days=TcompareTime0$days[4], CI=0.95)

    #For 0-24 hours, using 16O reps to calculate T0 abundance & 18O reps to calculate Tt abundance & net rates:
      all.rates.24h.18 <- create.all.rates.func(qSIP.data=all.comparisons.24, pop.data=all.comparisons.pop.24, Time.t.group.code="t24_18O", set.of.all.taxa=taxa.set.all)
      all.rates.24h.18 <- add.death.rates.func(all.rates=all.rates.24h.18, r.gross.boots=r.gross.boots.24, r.net.boots=r.net.boots.24, f.gross.boots=f.gross.boots.24, f.net.boots=f.net.boots.24, N.T0.boots=N.T0.boots.24, N.Tt.boots=N.Tt.boots.24, growth.model="linear", days=TcompareTime0$days[5], CI=0.95)

    #For 0-24 hours, using 16O reps to calculate T0 abundance & 16O+18O reps to calculate Tt abundance & net rates:
      all.rates.24h.16.18 <- create.all.rates.func(qSIP.data=all.comparisons.24, pop.data=all.comparisons.pop.24, Time.t.group.code="t24_16O.t24_18O", set.of.all.taxa=taxa.set.all)
      all.rates.24h.16.18 <- add.death.rates.func(all.rates=all.rates.24h.16.18, r.gross.boots=r.gross.boots.24, r.net.boots=r.net.boots.24, f.gross.boots=f.gross.boots.24, f.net.boots=f.net.boots.24, N.T0.boots=N.T0.boots.24, N.Tt.boots=N.Tt.boots.24, growth.model="linear", days=TcompareTime0$days[6], CI=0.95)

      #For 0-48 hours, using only 16O reps to calculate T0 & Tt abundances & net rates:
      all.rates.48h.16 <- create.all.rates.func(qSIP.data=all.comparisons.48, pop.data=all.comparisons.pop.48, Time.t.group.code="t48_16O", set.of.all.taxa=taxa.set.all)
      all.rates.48h.16 <- add.death.rates.func(all.rates=all.rates.48h.16, r.gross.boots=r.gross.boots.48, r.net.boots=r.net.boots.48, f.gross.boots=f.gross.boots.48, f.net.boots=f.net.boots.48, N.T0.boots=N.T0.boots.48, N.Tt.boots=N.Tt.boots.48, growth.model="linear", days=TcompareTime0$days[4], CI=0.95)
      
      #For 0-48 hours, using 16O reps to calculate T0 abundance & 18O reps to calculate Tt abundance & net rates:
      all.rates.48h.18 <- create.all.rates.func(qSIP.data=all.comparisons.48, pop.data=all.comparisons.pop.48, Time.t.group.code="t48_18O", set.of.all.taxa=taxa.set.all)
      all.rates.48h.18 <- add.death.rates.func(all.rates=all.rates.48h.18, r.gross.boots=r.gross.boots.48, r.net.boots=r.net.boots.48, f.gross.boots=f.gross.boots.48, f.net.boots=f.net.boots.48, N.T0.boots=N.T0.boots.48, N.Tt.boots=N.Tt.boots.48, growth.model="linear", days=TcompareTime0$days[5], CI=0.95)
      
      #For 0-48 hours, using 16O reps to calculate T0 abundance & 16O+18O reps to calculate Tt abundance & net rates:
      all.rates.48h.16.18 <- create.all.rates.func(qSIP.data=all.comparisons.48, pop.data=all.comparisons.pop.48, Time.t.group.code="t48_16O.t48_18O", set.of.all.taxa=taxa.set.all)
      all.rates.48h.16.18 <- add.death.rates.func(all.rates=all.rates.48h.16.18, r.gross.boots=r.gross.boots.48, r.net.boots=r.net.boots.48, f.gross.boots=f.gross.boots.48, f.net.boots=f.net.boots.48, N.T0.boots=N.T0.boots.48, N.Tt.boots=N.Tt.boots.48, growth.model="linear", days=TcompareTime0$days[6], CI=0.95)
      
      
    #For 0-72 hours, using only 16O reps to calculate T0 & Tt abundances & net rates:
      all.rates.72h.16 <- create.all.rates.func(qSIP.data=all.comparisons.72, pop.data=all.comparisons.pop.72, Time.t.group.code="t72_16O", set.of.all.taxa=taxa.set.all)
      all.rates.72h.16 <- add.death.rates.func(all.rates=all.rates.72h.16, r.gross.boots=r.gross.boots.72, r.net.boots=r.net.boots.72, f.gross.boots=f.gross.boots.72, f.net.boots=f.net.boots.72, N.T0.boots=N.T0.boots.72, N.Tt.boots=N.Tt.boots.72, growth.model="linear", days=TcompareTime0$days[7], CI=0.95)

    #For 0-72 hours, using 16O reps to calculate T0 abundance & 18O reps to calculate Tt abundance & net rates:
      all.rates.72h.18 <- create.all.rates.func(qSIP.data=all.comparisons.72, pop.data=all.comparisons.pop.72, Time.t.group.code="t72_18O", set.of.all.taxa=taxa.set.all)
      all.rates.72h.18 <- add.death.rates.func(all.rates=all.rates.72h.18, r.gross.boots=r.gross.boots.72, r.net.boots=r.net.boots.72, f.gross.boots=f.gross.boots.72, f.net.boots=f.net.boots.72, N.T0.boots=N.T0.boots.72, N.Tt.boots=N.Tt.boots.72, growth.model="linear", days=TcompareTime0$days[8], CI=0.95)

    #For 0-72 hours, using 16O reps to calculate T0 abundance & 16O+18O reps to calculate Tt abundance & net rates:
      all.rates.72h.16.18 <- create.all.rates.func(qSIP.data=all.comparisons.72, pop.data=all.comparisons.pop.72, Time.t.group.code="t72_16O.t72_18O", set.of.all.taxa=taxa.set.all)
      all.rates.72h.16.18 <- add.death.rates.func(all.rates=all.rates.72h.16.18, r.gross.boots=r.gross.boots.72, r.net.boots=r.net.boots.72, f.gross.boots=f.gross.boots.72, f.net.boots=f.net.boots.72, N.T0.boots=N.T0.boots.72, N.Tt.boots=N.Tt.boots.72, growth.model="linear", days=TcompareTime0$days[9], CI=0.95)

    #For 0-168 hours, using only 16O reps to calculate T0 & Tt abundances & net rates:
      all.rates.168h.16 <- create.all.rates.func(qSIP.data=all.comparisons.168, pop.data=all.comparisons.pop.168, Time.t.group.code="t168_16O", set.of.all.taxa=taxa.set.all)
      all.rates.168h.16 <- add.death.rates.func(all.rates=all.rates.168h.16, r.gross.boots=r.gross.boots.168, r.net.boots=r.net.boots.168, f.gross.boots=f.gross.boots.168, f.net.boots=f.net.boots.168, N.T0.boots=N.T0.boots.168, N.Tt.boots=N.Tt.boots.168, growth.model="linear", days=TcompareTime0$days[10], CI=0.95)

    #For 0-168 hours, using 16O reps to calculate T0 abundance & 18O reps to calculate Tt abundance & net rates:
      all.rates.168h.18 <- create.all.rates.func(qSIP.data=all.comparisons.168, pop.data=all.comparisons.pop.168, Time.t.group.code="t168_18O", set.of.all.taxa=taxa.set.all)
      all.rates.168h.18 <- add.death.rates.func(all.rates=all.rates.168h.18, r.gross.boots=r.gross.boots.168, r.net.boots=r.net.boots.168, f.gross.boots=f.gross.boots.168, f.net.boots=f.net.boots.168, N.T0.boots=N.T0.boots.168, N.Tt.boots=N.Tt.boots.168, growth.model="linear", days=TcompareTime0$days[11], CI=0.95)

    #For 0-168 hours, using 16O reps to calculate T0 abundance & 16O+18O reps to calculate Tt abundance & net rates:
      all.rates.168h.16.18 <- create.all.rates.func(qSIP.data=all.comparisons.168, pop.data=all.comparisons.pop.168, Time.t.group.code="t168_16O.t168_18O", set.of.all.taxa=taxa.set.all)
      all.rates.168h.16.18 <- add.death.rates.func(all.rates=all.rates.168h.16.18, r.gross.boots=r.gross.boots.168, r.net.boots=r.net.boots.168, f.gross.boots=f.gross.boots.168, f.net.boots=f.net.boots.168, N.T0.boots=N.T0.boots.168, N.Tt.boots=N.Tt.boots.168, growth.model="linear", days=TcompareTime0$days[12], CI=0.95)

    #For the 3-24 hour interval, using only 18O reps to calculate T0 & Tt abundances & net rates:
      all.rates.3h.24h <- create.all.rates.interval.func(pop.data=all.comparisons.pop.3.24, set.of.all.taxa=taxa.set.all)
      interval.3h.24h.list <- add.birth.flux.interval.func(all.rates=all.rates.3h.24h, T1.all.rates=all.rates.3h.18, T2.all.rates=all.rates.24h.18, T1.r.gross.boots=r.gross.boots.3, T2.r.gross.boots=r.gross.boots.24, T1.N.Tt.boots=N.Tt.boots.3, T2.N.Tt.boots=N.Tt.boots.24, T2.f.gross.boots=f.gross.boots.24, growth.model="linear", interval.days=TcompareTime0$days[13], T1.days=TcompareTime0$days[2], T2.days=TcompareTime0$days[5], copies.cell=6, pgC.cell=0.1, CI=0.95)
      names(interval.3h.24h.list)
      all.rates.3h.24h <- interval.3h.24h.list$all.rates
      all.rates.3h.24h <- add.death.rates.func(all.rates=all.rates.3h.24h, r.gross.boots=interval.3h.24h.list$r.gross.boots, r.net.boots=r.net.boots.3.24, f.gross.boots=interval.3h.24h.list$f.gross.boots, f.net.boots=f.net.boots.3.24, N.T0.boots=N.T0.boots.3.24, N.Tt.boots=N.Tt.boots.3.24, growth.model="linear", days=TcompareTime0$days[13], CI=0.95)
 
    #For the 24-48 hour interval, using only 18O reps to calculate T0 & Tt abundances & net rates:
      all.rates.24h.48h <- create.all.rates.interval.func(pop.data=all.comparisons.pop.24.48, set.of.all.taxa=taxa.set.all)
      interval.24h.48h.list <- add.birth.flux.interval.func(all.rates=all.rates.24h.48h, T1.all.rates=all.rates.24h.18, T2.all.rates=all.rates.48h.18, T1.r.gross.boots=r.gross.boots.24, T2.r.gross.boots=r.gross.boots.48, T1.N.Tt.boots=N.Tt.boots.24, T2.N.Tt.boots=N.Tt.boots.48, T2.f.gross.boots=f.gross.boots.48, growth.model="linear", interval.days=TcompareTime0$days[14], T1.days=TcompareTime0$days[5], T2.days=TcompareTime0$days[8], copies.cell=6, pgC.cell=0.1, CI=0.95)
      names(interval.24h.48h.list)
      all.rates.24h.48h <- interval.24h.48h.list$all.rates
      all.rates.24h.48h <- add.death.rates.func(all.rates=all.rates.24h.48h, r.gross.boots=interval.24h.48h.list$r.gross.boots, r.net.boots=r.net.boots.24.48, f.gross.boots=interval.24h.48h.list$f.gross.boots, f.net.boots=f.net.boots.24.48, N.T0.boots=N.T0.boots.24.48, N.Tt.boots=N.Tt.boots.24.48, growth.model="linear", days=TcompareTime0$days[14], CI=0.95)

      #For the 48-72 hour interval, using only 18O reps to calculate T0 & Tt abundances & net rates:
      all.rates.48h.72h <- create.all.rates.interval.func(pop.data=all.comparisons.pop.48.72, set.of.all.taxa=taxa.set.all)
      interval.48h.72h.list <- add.birth.flux.interval.func(all.rates=all.rates.48h.72h, T1.all.rates=all.rates.48h.18, T2.all.rates=all.rates.72h.18, T1.r.gross.boots=r.gross.boots.48, T2.r.gross.boots=r.gross.boots.72, T1.N.Tt.boots=N.Tt.boots.48, T2.N.Tt.boots=N.Tt.boots.72, T2.f.gross.boots=f.gross.boots.72, growth.model="linear", interval.days=TcompareTime0$days[14], T1.days=TcompareTime0$days[5], T2.days=TcompareTime0$days[8], copies.cell=6, pgC.cell=0.1, CI=0.95)
      names(interval.48h.72h.list)
      all.rates.48h.72h <- interval.48h.72h.list$all.rates
      all.rates.48h.72h <- add.death.rates.func(all.rates=all.rates.48h.72h, r.gross.boots=interval.48h.72h.list$r.gross.boots, r.net.boots=r.net.boots.48.72, f.gross.boots=interval.48h.72h.list$f.gross.boots, f.net.boots=f.net.boots.48.72, N.T0.boots=N.T0.boots.48.72, N.Tt.boots=N.Tt.boots.48.72, growth.model="linear", days=TcompareTime0$days[14], CI=0.95)
      
    #For the 72-168 hour interval, using only 18O reps to calculate T0 & Tt abundances & net rates:
      all.rates.72h.168h <- create.all.rates.interval.func(pop.data=all.comparisons.pop.72.168, set.of.all.taxa=taxa.set.all)
      interval.72h.168h.list <- add.birth.flux.interval.func(all.rates=all.rates.72h.168h, T1.all.rates=all.rates.72h.18, T2.all.rates=all.rates.168h.18, T1.r.gross.boots=r.gross.boots.72, T2.r.gross.boots=r.gross.boots.168, T1.N.Tt.boots=N.Tt.boots.72, T2.N.Tt.boots=N.Tt.boots.168, T2.f.gross.boots=f.gross.boots.168, growth.model="linear", interval.days=TcompareTime0$days[15], T1.days=TcompareTime0$days[8], T2.days=TcompareTime0$days[11], copies.cell=6, pgC.cell=0.1, CI=0.95)
      names(interval.72h.168h.list)
      all.rates.72h.168h <- interval.72h.168h.list$all.rates
      all.rates.72h.168h <- add.death.rates.func(all.rates=all.rates.72h.168h, r.gross.boots=interval.72h.168h.list$r.gross.boots, r.net.boots=r.net.boots.72.168, f.gross.boots=interval.72h.168h.list$f.gross.boots, f.net.boots=f.net.boots.72.168, N.T0.boots=N.T0.boots.72.168, N.Tt.boots=N.Tt.boots.72.168, growth.model="linear", days=TcompareTime0$days[15], CI=0.95)

    #For the three 'interval-based' sets of results, the estimates of 'N.light.Tt' differ from those in the corresponding 'serial' sets of results
    #As yet, I'm not 100% sure why this is.  It could be due to differences in random bootstrap resampling in all.taxa.calcs, or (more likely) it may be due to the fact that N.light.Tt cannot be calculated the same way for interval that do not start at zero, or it could be due to the 0.5*detection limit assignment done above for these 'interval-based' calculations.
    #Therefore, get rid of the 'N.light.Tt' estimates from the 'interval-based' results to avoid confusion (if needed, use the 'N.light.Tt' estimates from the corresponding 'serial' sets fo results instead:
      #Take a quick look at the differences:
        #24 hours:
          par(mfrow=c(2,2))
          plot(x=all.rates.24h.18$N.light.Tt.obs, y=all.rates.3h.24h$N.light.Tt.obs, xlab="all.rates.24h.18", ylab="all.rates.3h.24h", main='N.light.Tt.obs (24h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.24h.18$N.light.Tt.boot.median, y=all.rates.3h.24h$N.light.Tt.boot.median, xlab="all.rates.24h.18", ylab="all.rates.3h.24h", main='N.light.Tt.boot.median (24h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.24h.18$N.light.Tt.boot.CI.L, y=all.rates.3h.24h$N.light.Tt.boot.CI.L, xlab="all.rates.24h.18", ylab="all.rates.3h.24h", main='N.light.Tt.boot.CI.L (24h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.24h.18$N.light.Tt.boot.CI.U, y=all.rates.3h.24h$N.light.Tt.boot.CI.U, xlab="all.rates.24h.18", ylab="all.rates.3h.24h", main='N.light.Tt.boot.CI.U (24h)')
          abline(a=0, b=1, col="blue")      
          par(mfrow=c(1,1))
        #48 hours:
          par(mfrow=c(2,2))
          plot(x=all.rates.48h.18$N.light.Tt.obs, y=all.rates.24h.48h$N.light.Tt.obs, xlab="all.rates.48h.18", ylab="all.rates.24h.48h", main='N.light.Tt.obs (48h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.48h.18$N.light.Tt.boot.median, y=all.rates.24h.48h$N.light.Tt.boot.median, xlab="all.rates.48h.18", ylab="all.rates.24h.48h", main='N.light.Tt.boot.median (48h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.48h.18$N.light.Tt.boot.CI.L, y=all.rates.24h.48h$N.light.Tt.boot.CI.L, xlab="all.rates.48h.18", ylab="all.rates.24h.48h", main='N.light.Tt.boot.CI.L (48h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.48h.18$N.light.Tt.boot.CI.U, y=all.rates.24h.48h$N.light.Tt.boot.CI.U, xlab="all.rates.48h.18", ylab="all.rates.24h.48h", main='N.light.Tt.boot.CI.U (48h)')
          abline(a=0, b=1, col="blue")      
          par(mfrow=c(1,1))
        #168 hours:
          par(mfrow=c(2,2))
          plot(x=all.rates.168h.18$N.light.Tt.obs, y=all.rates.72h.168h$N.light.Tt.obs, xlab="all.rates.168h.18", ylab="all.rates.72h.168h", main='N.light.Tt.obs (168h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.168h.18$N.light.Tt.boot.median, y=all.rates.72h.168h$N.light.Tt.boot.median, xlab="all.rates.168h.18", ylab="all.rates.72h.168h", main='N.light.Tt.boot.median (168h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.168h.18$N.light.Tt.boot.CI.L, y=all.rates.72h.168h$N.light.Tt.boot.CI.L, xlab="all.rates.168h.18", ylab="all.rates.72h.168h", main='N.light.Tt.boot.CI.L (168h)')
          abline(a=0, b=1, col="blue")
          plot(x=all.rates.168h.18$N.light.Tt.boot.CI.U, y=all.rates.72h.168h$N.light.Tt.boot.CI.U, xlab="all.rates.168h.18", ylab="all.rates.72h.168h", main='N.light.Tt.boot.CI.U (168h)')
          abline(a=0, b=1, col="blue")      
          par(mfrow=c(1,1))
      #Replace all the 'N.light.Tt' entries in the interval-based sets of results with NA's:
        #3-24hours:
          all.rates.3h.24h$N.light.Tt.obs <- NA
          all.rates.3h.24h$N.light.Tt.boot.mean <- NA
          all.rates.3h.24h$N.light.Tt.boot.median <- NA
          all.rates.3h.24h$N.light.Tt.boot.CI.L <- NA
          all.rates.3h.24h$N.light.Tt.boot.CI.U <- NA
        #24-48hours:
          all.rates.24h.48h$N.light.Tt.obs <- NA
          all.rates.24h.48h$N.light.Tt.boot.mean <- NA
          all.rates.24h.48h$N.light.Tt.boot.median <- NA
          all.rates.24h.48h$N.light.Tt.boot.CI.L <- NA
          all.rates.24h.48h$N.light.Tt.boot.CI.U <- NA
        #72-168hours:
          all.rates.72h.168h$N.light.Tt.obs <- NA
          all.rates.72h.168h$N.light.Tt.boot.mean <- NA
          all.rates.72h.168h$N.light.Tt.boot.median <- NA
          all.rates.72h.168h$N.light.Tt.boot.CI.L <- NA
          all.rates.72h.168h$N.light.Tt.boot.CI.U <- NA

  #Distributions of the 'serial' growth rates:
    dev.off()
    dev.new(width=11, height=7)
    par(mfcol=c(3,4))
    hist(all.rates.3h.16.18[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.3h.16.18[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.3h.16.18[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    hist(all.rates.24h.16.18[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.24h.16.18[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.24h.16.18[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    hist(all.rates.72h.16.18[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.72h.16.18[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.72h.16.18[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    hist(all.rates.168h.16.18[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.168h.16.18[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.168h.16.18[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    par(mfcol=c(1,1))

  #Distributions of the 'interval' growth rates:
    dev.off()
    dev.new(width=11, height=7)
    par(mfcol=c(3,3))
    hist(all.rates.3h.24h[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.3h.24h[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.3h.24h[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    hist(all.rates.24h.48h[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.24h.48h[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.24h.48h[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    hist(all.rates.72h.168h[,c(10)], xlim=c(-1.5e+10, 2.5e+10), main="birth")
    hist(all.rates.72h.168h[,c(24)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="net")
    hist(all.rates.72h.168h[,c(45)], xlim=c(-1.5e+10, 2.5e+10), breaks=20, main="death")
    par(mfcol=c(1,1))

    dev.off()


#Write all 15 sets of results to text files:
  write.table(all.rates.3h.16, paste(Dir.o, "all_rates_3h_16.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.3h.18, paste(Dir.o, "all_rates_3h_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.3h.16.18, paste(Dir.o, "all_rates_3h_16_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  write.table(all.rates.24h.16, paste(Dir.o, "all_rates_24h_16.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.24h.18, paste(Dir.o, "all_rates_24h_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.24h.16.18, paste(Dir.o, "all_rates_24h_16_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  write.table(all.rates.48h.16, paste(Dir.o, "all_rates_48h_16.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.48h.18, paste(Dir.o, "all_rates_48h_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.48h.16.18, paste(Dir.o, "all_rates_48h_16_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  
  write.table(all.rates.72h.16, paste(Dir.o, "all_rates_72h_16.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.72h.18, paste(Dir.o, "all_rates_72h_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.72h.16.18, paste(Dir.o, "all_rates_72h_16_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  write.table(all.rates.168h.16, paste(Dir.o, "all_rates_168h_16.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.168h.18, paste(Dir.o, "all_rates_168h_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.168h.16.18, paste(Dir.o, "all_rates_168h_16_18.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  write.table(all.rates.3h.24h, paste(Dir.o, "all_rates_3h_24h.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.24h.48h, paste(Dir.o, "all_rates_24h_48h.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.48h.72h, paste(Dir.o, "all_rates_48h_72h.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(all.rates.72h.168h, paste(Dir.o, "all_rates_72h_168h.txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/100-SB_01-02.RData", envir=.GlobalEnv)
    


