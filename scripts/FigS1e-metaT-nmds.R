##Linnea Honeker
##NMDS and PCA on metaT data from IMG
##2/26/24

#import libraries
library(ggfortify)
library(tidyverse)
library(factoextra)
library(ggnewscale)
library(ggrepel)
library(viridis)
library(vegan)
library(nlme)
library(ecodist)



#set aesthetics
col_list = c("red", "orange", "green",  "blue",  "purple", "pink")
col_list_subject = c("turquoise3", "magenta3", "red3", "darkgreen", "yellow3", "blue3", "orange3", "grey")
col_list_moisture <- c("darkgreen", "brown")
list_of_shapes = c(17, 16,15,18,19, 20, 21)

# parameters for plots
h = 6
w = 7
res = 300
size = 10

# set figure output directories
out_dir <- "./figures/FigS1-NMDS/"

#import metaT (ko)
metaT.0 <- read.csv("./output/metaT/vst-ko.csv", header = TRUE)
metaT <- metaT.0 %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()

#import metaT (ko) no time 0
metaT.no.t0.0 <- read.csv("./output/metaT/vst-ko-no-t0.csv", header = TRUE)
metaT.no.t0 <- metaT.no.t0.0 %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()

#import metaT (cazy) 
metaT.cazy.0 <- read.csv("./output/metaT/vst-cazy.csv", header = TRUE)
metaT.cazy <- metaT.cazy.0 %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()

#import metaT (cazy) 
metaT.cazy.no.t0.0 <- read.csv("./output/metaT/vst-cazy-no-t0.csv", header = TRUE)
metaT.cazy.no.t0 <- metaT.cazy.no.t0.0 %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()

#import vip
kegg.vip <- read.csv("./output/PLSR/TableS1_10_v2.csv", header = TRUE) 
cazy.vip <- read.csv("./output/PLSR/VIP_cazy_top_10.csv", header = TRUE)

#filter dfs to only vip cutoff
metaT.vip <- metaT.no.t0 %>%
  select(kegg.vip$KO)

metaT.cazy.vip <- metaT.cazy.no.t0 %>%
  select(cazy.vip$cazy)


########################################## nmds on metaT-mapped to gene catalog#####################
# Calculate nmds
set.seed(123)
nmds = metaMDS(metaT, distance = "bray")
plot(nmds)


scores <- vegan::scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  split_sample_id("SampleID") 

### nmds plot (timepoint and moisture)
nmds_plot <-  make_nmds_plot_time(data.scores)

filename <- paste0(out_dir, "metaT_nmds_ko.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(out_dir,"metaT_nmds_ko.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

### nmds plot (subject_id and moisture)
nmds_plot <-  make_nmds_plot_moisture(data.scores)

filename <- paste0(out_dir, "metaT_nmds_moisture_ko.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(out_dir, "metaT_nmds_moisture_ko.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

nmds_plot <- make_nmds_plot_plot(data.scores)

####statistical tests#####
##metaT###
df.df <- as.data.frame(metaT)

df.df <- df.df[sapply(df.df, is.numeric)]


df.stat <- df.df %>% 
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c("subject_id", "moisture", "timepoint"), sep = "_", remove = T) %>%
  mutate(CGE = case_when(
    timepoint == "t24" ~ "high",
    moisture == "x50" ~ "low",
    moisture == "x100" ~ "high")
)

####permanova on bray-curtis distance matrix####
##all time points##
dist<- ecodist::distance(df.df, "bray-curtis")

adonis2(dist ~ timepoint*moisture,  df.stat)
# adonis2(formula = dist ~ timepoint * moisture, data = df.stat)
#         Df SumOfSqs      R2      F Pr(>F)    
# Model    11 0.103158 0.68799 6.8154  0.001 ***
#   Residual 34 0.046784 0.31201                  
# Total    45 0.149942 1.00000     
adonis2(dist ~ CGE,  df.stat)

adonis2(dist ~ moisture,  df.stat)
# adonis2(formula = dist ~ moisture, data = df.stat)
#         Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.00149 0.00994 0.4416  0.941
# Residual 44  0.14845 0.99006              
# Total    45  0.14994 1.00000

adonis2(dist ~ timepoint,  df.stat)
# adonis2(formula = dist ~ timepoint, data = df.stat)
# D.      f SumOfSqs      R2      F Pr(>F)    
# Model     5 0.095261 0.63532 13.937  0.001 ***
#   Residual 40 0.054681 0.36468                  
# Total    45 0.149942 1.00000   

adonis2(dist ~ subject_id, df.stat)
# adonis2(formula = dist ~ subject_id, data = df.stat)
#        Df SumOfSqs      R2      F Pr(>F)
# Model     7  0.01529 0.10198 0.6164  0.986
# Residual 38  0.13465 0.89802              
# Total    45  0.14994 1.00000    


########################################## nmds on metaT mapped to gc cazy#####################
# Calculate nmds
set.seed(123)
nmds = metaMDS(metaT.cazy, distance = "bray", autotransform = FALSE)
plot(nmds)

scores <- vegan::scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  separate(., SampleID, into = c("plot", "moisture", "timepoint"), sep = "_", remove = T) 

data.scores$timepoint <- factor(data.scores$timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168"))
data.scores$moisture <- factor(data.scores$moisture, levels = c("x50", "x100"))



### nmds plot (timepoint and moisture)
nmds_plot <-  make_nmds_plot_time(data.scores)

filename <- paste0(out_dir, "metaT_nmds_cazy.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

filename <- paste0(out_dir, "metaT_nmds_cazy.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

stressplot(nmds)
nmds$stress

### nmds plot (moisture)
nmds_plot <-  make_nmds_plot_moisture(data.scores)

filename <- paste0(out_dir, "metaT_nmds_cazy_moisture.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(out_dir, "metaT_nmds_cazy_moisture.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)


####statistical tests#####
##metaT###
#df.df <- as.data.frame(df)
df.df <- as.data.frame(metaT.cazy)

df.df <- df.df[sapply(df.df, is.numeric)]


df.stat <- df.df %>% 
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c("subject_id", "moisture", "timepoint"), sep = "_", remove = T) 


####permanova on bray-curtis distance matrix####
##all time points##
dist<- ecodist::distance(df.df, "bray-curtis")

adonis2(dist ~ timepoint*moisture,  df.stat)
# adonis2(formula = dist ~ timepoint * moisture, data = df.stat)
# Df SumOfSqs      R2     F Pr(>F)    
# Model    11 0.056370 0.67316 6.366  0.001 ***
#   Residual 34 0.027369 0.32684                 
# Total    45 0.083739 1.00000                 
# ---

adonis2(dist ~ timepoint,  df.stat)
# adonis2(formula = dist ~ timepoint, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     5 0.051789 0.61845 12.967  0.001 ***
#   Residual 40 0.031951 0.38155                  
# Total    45 0.083739 1.00000   

adonis2(dist ~ moisture,  df.stat)
# adonis2(formula = dist ~ moisture, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)
# Model     1 0.001351 0.01613 0.7215    0.6
# Residual 44 0.082388 0.98387              
# Total    45 0.083739 1.00000        

adonis2(dist ~ subject_id, df.stat)
# adonis2(formula = dist ~ subject_id, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)
# Model     7 0.009031 0.10784 0.6562  0.967
# Residual 38 0.074709 0.89216              
# Total    45 0.083739 1.00000    

###########################################only genes associated with CGE############################
########################################## nmds on metaT-mapped to gene catalog#####################
# Calculate nmds
set.seed(123)
nmds = metaMDS(metaT.vip, distance = "bray")
plot(nmds)


scores <- vegan::scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  split_sample_id("SampleID") 

### nmds plot (timepoint and moisture)
nmds_plot <-  make_nmds_plot_time_no_t0(data.scores)

filename <- paste0(out_dir, "metaT_nmds_ko_vip.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(out_dir,"metaT_nmds_ko_vip.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

### nmds plot (subject_id and moisture)
nmds_plot <-  make_nmds_plot_moisture(data.scores)

filename <- paste0(out_dir, "metaT_nmds_moisture_ko_vip.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(out_dir, "metaT_nmds_moisture_ko_vip.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

nmds_plot <- make_nmds_plot_plot(data.scores)

####statistical tests#####
##metaT###
df.df <- as.data.frame(metaT.vip)

df.df <- df.df[sapply(df.df, is.numeric)]


df.stat <- df.df %>% 
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c("subject_id", "moisture", "timepoint"), sep = "_", remove = T) 


####permanova on bray-curtis distance matrix####
##all time points##
dist<- ecodist::distance(df.df, "bray-curtis")

adonis2(dist ~ timepoint*moisture,  df.stat)
# adonis2(formula = dist ~ timepoint * moisture, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     9 0.035536 0.56119 4.1209  0.001 ***
#   Residual 29 0.027786 0.43881                  
# Total    38 0.063323 1.00000  


adonis2(dist ~ moisture,  df.stat)
# adonis2(formula = dist ~ moisture, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1 0.003585 0.05662 2.2205  0.033 *
#   Residual 37 0.059737 0.94338                
# Total    38 0.063323 1.00000 

adonis2(dist ~ timepoint,  df.stat)
# adonis2(formula = dist ~ timepoint, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     4 0.028892 0.45627 7.1328  0.001 ***
#   Residual 34 0.034430 0.54373                  
# Total    38 0.063323 1.00000      

adonis2(dist ~ subject_id, df.stat)
# adonis2(formula = dist ~ subject_id, data = df.stat)
# Df SumOfSqs      R2    F Pr(>F)
# Model     7 0.012322 0.19459 1.07  0.333
# Residual 31 0.051000 0.80541            
# Total    38 0.063323 1.00000 


########################################## nmds on metaT mapped to gc cazy#####################
# Calculate nmds
set.seed(123)
nmds = metaMDS(metaT.cazy.vip, distance = "bray", autotransform = FALSE)
plot(nmds)

scores <- vegan::scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  separate(., SampleID, into = c("plot", "moisture", "timepoint"), sep = "_", remove = T) 

data.scores$timepoint <- factor(data.scores$timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168"))
data.scores$moisture <- factor(data.scores$moisture, levels = c("x50", "x100"))



### nmds plot (timepoint and moisture)
nmds_plot <-  make_nmds_plot_time_no_t0(data.scores)

filename <- paste0(out_dir, "metaT_nmds_cazy_vip.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

filename <- paste0(out_dir, "metaT_nmds_cazy-vip.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

stressplot(nmds)
nmds$stress

### nmds plot (moisture)
nmds_plot <-  make_nmds_plot_moisture(data.scores)

filename <- paste0(out_dir, "metaT_nmds_cazy_moisture_vip.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(out_dir, "metaT_nmds_cazy_moisture_vip.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)


####statistical tests#####
##metaT###
#df.df <- as.data.frame(df)
df.df <- as.data.frame(metaT.cazy.vip)

df.df <- df.df[sapply(df.df, is.numeric)]


df.stat <- df.df %>% 
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c("subject_id", "moisture", "timepoint"), sep = "_", remove = T) 


####permanova on bray-curtis distance matrix####
##all time points##
dist<- ecodist::distance(df.df, "bray-curtis")

adonis2(dist ~ timepoint*moisture,  df.stat)
# adonis2(formula = dist ~ timepoint * moisture, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     9 0.027453 0.55556 4.0278  0.001 ***
#   Residual 29 0.021962 0.44444                  
# Total    38 0.049415 1.00000 

adonis2(dist ~ timepoint,  df.stat)
# adonis2(formula = dist ~ timepoint, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     4 0.021403 0.43313 6.4945  0.001 ***
#   Residual 34 0.028012 0.56687                  
# Total    38 0.049415 1.00000  

adonis2(dist ~ moisture,  df.stat)
# adonis2(formula = dist ~ moisture, data = df.stat)
# Df SumOfSqs      R2     F Pr(>F)  
# Model     1 0.003423 0.06928 2.754  0.025 *
#   Residual 37 0.045992 0.93072               
# Total    38 0.049415 1.00000        

adonis2(dist ~ subject_id, df.stat)
# adonis2(formula = dist ~ subject_id, data = df.stat)
# Df SumOfSqs      R2      F Pr(>F)
# Model     7 0.009188 0.18594 1.0115  0.443
# Residual 31 0.040227 0.81406              
# Total    38 0.049415 1.00000   



