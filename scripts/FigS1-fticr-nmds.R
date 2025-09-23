##Linnea Hernandez
##NMDS on FTICR data
##2/26/24


library(tidyverse)
#library(dplyr)
library(vegan)
library(factoextra)
#library(ggrepel)
#library(ggpubr)
#library(ggnewscale)
library('viridis')


#define colors and shapes
col_list = c("red", "orange", "green",  "blue",  "purple", "pink")
col_list_subject = c("red3", "darkgreen", "yellow3", "blue3", "orange3", "grey")
col_list_moisture <- c("darkgreen", "brown")
list_of_shapes = c(17, 16,15,18,19, 20, 21)

# parameters for plots
h = 6
w = 7
res = 300
size = 10


# load source code
source("./scripts/functions/NMDS-functions.R")

# Import data ----

## Defining paths
Dir.f <- "./figures/FigS1-NMDS/"


## Loading tables, transpose, convert to matrix, replace fticr sample names with new names
fticr.h2o <- read.csv("./output/H2O-fticr/metabodirect/Report_processed_noNorm_MolecFormulas.csv", header = TRUE) %>%
  select(-c(C, H, O, N, S, P, NeutralMass, Error_ppm, MolecularFormula, El_comp, Class, OC, HC, NOSC, GFE, DBE, DBE_O, AI, AI_mod, DBE_AI)) %>%
  column_to_rownames(var = "Mass") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  filter(!grepl("DI", SampleID)) %>% #remove any DI controls
  separate(SampleID, c("Sample", "number", "othernumber", "otherothernumber")) %>%
  select(-c(number, othernumber, otherothernumber)) %>%
  separate(Sample, c("time", "other"), sep = "h") %>%
  separate(other, c("moisture", "subjectID"), sep = "pct") %>%
  mutate(time = str_replace(time, "X0", "t0")) %>%
  mutate(time = str_replace(time, "X3", "t3")) %>%
  mutate(time = str_replace(time, "X24", "t24")) %>%
  mutate(time = str_replace(time, "X48", "t48")) %>%
  mutate(time = str_replace(time, "X72", "t72")) %>%
  mutate(time = str_replace(time, "X168", "t168")) %>%
  mutate(moisture = str_replace(moisture, "100", "x100")) %>%
  mutate(moisture = str_replace(moisture, "50", "x50")) %>%
  mutate(subjectID = str_remove(subjectID, "a")) %>%
  mutate(subjectID = str_remove(subjectID, "b")) %>%
  unite(SampleID, c(subjectID, moisture, time)) %>%
  column_to_rownames(var = "SampleID")

fticr.meoh <- read.csv("./output/MeOH-fticr/metabodirect/Report_processed_noNorm_MolecFormulas.csv", header = TRUE) %>%
  select(-c(C, H, O, N, S, P, NeutralMass, Error_ppm, MolecularFormula, El_comp, Class, OC, HC, NOSC, GFE, DBE, DBE_O, AI, AI_mod, DBE_AI)) %>%
  column_to_rownames(var = "Mass") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  filter(!grepl("DI", SampleID)) %>% #remove any DI controls
  separate(SampleID, c("Sample", "number", "othernumber", "otherothernumber")) %>%
  select(-c(number, othernumber, otherothernumber)) %>%
  separate(Sample, c("time", "other"), sep = "h") %>%
  separate(other, c("moisture", "subjectID"), sep = "pct") %>%
  mutate(time = str_replace(time, "X0", "t0")) %>%
  mutate(time = str_replace(time, "X3", "t3")) %>%
  mutate(time = str_replace(time, "X24", "t24")) %>%
  mutate(time = str_replace(time, "X48", "t48")) %>%
  mutate(time = str_replace(time, "X72", "t72")) %>%
  mutate(time = str_replace(time, "X168", "t168")) %>%
  mutate(moisture = str_replace(moisture, "100", "x100")) %>%
  mutate(moisture = str_replace(moisture, "50", "x50")) %>%
  mutate(subjectID = str_remove(subjectID, "a")) %>%
  mutate(subjectID = str_remove(subjectID, "b")) %>%
  unite(SampleID, c(subjectID, moisture, time)) %>%
  column_to_rownames(var = "SampleID")


#convert to matrices
#m.h2o <- fticr.h2o %>%
#  select(-X) %>%
#  column_to_rownames(var = "SampleID")

m.meoh <- fticr.meoh %>%
  filter(rownames(.) != "P13_x50_t168")  #remove out

m.h2o <- fticr.h2o



# NMDS on m.h20

# Calculate nmds
#first, make non-negative by adding lowest number plus a few extra to all values to make positive
m.h2o_bin <- m.h2o > 0  # Converts abundance to presence/absence

set.seed(123)
nmds = metaMDS(m.h2o_bin, distance = "jaccard")
plot(nmds)

scores <- vegan::scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  separate(., SampleID, into = c("subject_id", "moisture", "timepoint"), sep = "_", remove = T) 

data.scores$timepoint <- factor(data.scores$timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168"))
data.scores$moisture <- factor(data.scores$moisture, levels = c("x50", "x100"))

### nmds plot (timepoint and moisture)
nmds_plot <- make_nmds_plot_time(data.scores)

filename <- paste0(Dir.f, "fticr-h2o-nmds-timepoint.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(Dir.f, "fticr-h2o-nmds-timepoint.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)


stressplot(nmds)
nmds$stress
#0.1160701

## NMDS plot ----moisture

nmds_plot <-  make_nmds_plot_moisture(data.scores)

filename <- paste0(Dir.f, "fticr-h2o-nmds-moisture.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(Dir.f, "fticr-h2o-nmds-moisture.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)



# PERMANOVA ----

#permanova
dist<- ecodist::distance(m.h2o_bin, "jaccard")

adonis2(dist ~ timepoint*moisture,  data.scores)
# adonis2(formula = dist ~ timepoint * moisture, data = data.scores)
# Df SumOfSqs     R2     F Pr(>F)   
# Model    11  0.30361 0.4144 1.544  0.003 **
#   Residual 24  0.42903 0.5856                
# Total    35  0.73264 1.0000               

adonis2(dist ~ moisture,  data.scores)
# adonis2(formula = dist ~ moisture, data = data.scores)
# Df SumOfSqs      R2      F Pr(>F)   
# Model     1  0.05646 0.07706 2.8388  0.006 **
#   Residual 34  0.67618 0.92294                 
# Total    35  0.73264 1.00000                   

adonis2(dist ~ timepoint,  data.scores)
# adonis2(formula = dist ~ timepoint, data = data.scores)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     5  0.18360 0.25059 2.0064  0.001 ***
#   Residual 30  0.54904 0.74941                  
# Total    35  0.73264 1.00000                  
# ---

adonis2(dist ~ subject_id, data.scores)
# adonis2(formula = dist ~ subject_id, data = data.scores)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     5  0.21466 0.29299 2.4865  0.001 ***
#   Residual 30  0.51798 0.70701                  
# Total    35  0.73264 1.00000     


###### NMDS on m.meoh
# Calculate nmds
m.meoh_bin <- m.meoh > 0  # Converts abundance to presence/absence

set.seed(123)
nmds = metaMDS(m.meoh_bin, distance = "jaccard")
plot(nmds)

scores <- vegan::scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples) %>%
  rownames_to_column(var = "SampleID") %>%
  #filter(SampleID != "P13_x50_t168") %>% #remove outlier
  column_to_rownames(var = "SampleID")

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  separate(., SampleID, into = c("plot", "moisture", "timepoint"), sep = "_", remove = T) 

data.scores$timepoint <- factor(data.scores$timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168"))
data.scores$moisture <- factor(data.scores$moisture, levels = c("x50", "x100"))

### nmds plot (timepoint)
nmds_plot <- make_nmds_plot_time(data.scores)

filename <- paste0(Dir.f, "fticr-meoh-nmds-timepoint.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(Dir.f, "fticr-meoh-nmds-timepoint.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)


stressplot(nmds)
nmds$stress
# [1] 0.1259875

## NMDS plot ----moisture

nmds_plot <-  make_nmds_plot_moisture(data.scores)

filename <- paste0(Dir.f, "fticr-meoh-nmds-moisture.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0(Dir.f, "fticr-meoh-nmds-moisture.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)



# PERMANOVA ----

#permanova
dist<- ecodist::distance(m.meoh_bin, "jaccard")

adonis2(dist ~ timepoint*moisture,  data.scores)
# adonis2(formula = dist ~ timepoint * moisture, data = data.scores)
# Df SumOfSqs      R2      F Pr(>F)  
# Model    11  0.76272 0.38651 1.2028  0.055 .
# Residual 21  1.21064 0.61349                
# Total    32  1.97336 1.00000   

adonis2(dist ~ moisture,  data.scores)
# adonis2(formula = dist ~ moisture, data = data.scores)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     1  0.12296 0.06231 2.0599  0.014 *
#   Residual 31  1.85040 0.93769                
# Total    32  1.97336 1.00000      

adonis2(dist ~ timepoint,  data.scores)
#adonis2(formula = dist ~ timepoint, data = data.scores)
#Df SumOfSqs      R2      F Pr(>F)  
#Model     5   0.5725 0.20274 1.3732  0.013 *
#  Residual 27   2.2513 0.79726                
#Total    32   2.8238 1.00000          


adonis2(dist ~ plot, data.scores)
# adonis2(formula = dist ~ timepoint, data = data.scores)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     5  0.39772 0.20155 1.3631  0.031 *
#   Residual 27  1.57564 0.79845                
# Total    32  1.97336 1.00000       


