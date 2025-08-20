##Linnea Honeker
##NMDS and PCA on metaT data from IMG
##2/26/24

library(ggfortify)
library(tidyverse)
library(factoextra)
library(ggnewscale)
library(ggrepel)
library(viridis)
library(vegan)
library(nlme)
library(ecodist)
library(ggplot2)
library(stats)


source("./qSIP_scripts/extra_functions_PCA_NMDS_plotting.R")

# parameters for plots
h = 6
w = 8
res = 300
size = 10

###inport data
#processed with molecular formulas
fticr.0 <- read.csv("./metabodirect/H2O/1_preprocessing_output/Report_processed_noNorm.csv", header = TRUE)


new.id.meta <- read.csv("./qSIP_data/FTICR/metadata-H20-no-extra.csv") %>%
  separate(SampleID, c("ID", "other"), sep = "_", remove = FALSE) %>%
  tidyr::unite(New_id, c("ID", "Zone", "Treatment", "Timepoint"), remove = FALSE) %>%
  select(-other) %>%
  column_to_rownames(var = "New_id") %>%
  mutate(Zone = str_replace(Zone, "Rhizosphere", "Rhizo"))

new.id <- new.id.meta %>%
  rownames_to_column(var = "New_id") %>%
  select(New_id, SampleID)
        

fticr <- fticr.0 %>%
  select(-c(C, H, O, N, S, P, NeutralMass, Error_ppm, El_comp, OC, HC, NOSC, GFE, DBE, DBE_O, AI, AI_mod, DBE_AI)) %>%
  #filter(MolecularFormula != "") %>%
  unite(Mass_class, c("Mass", "MolecularFormula", "Class"), sep = "_") %>%
  column_to_rownames(var = "Mass_class") 
  
#determine threshold
limit =0.005*ncol(fticr) 

#filter based on threshold
fticr.f <- fticr %>%
  mutate(sum.count = rowSums(.)) %>%
  mutate(gr.zero = rowSums(. > 0)) %>%
  filter(gr.zero > limit) %>%
  filter(gr.zero != ncol(fticr)) %>%
  mutate(sum.count.per = sum.count/sum(sum.count)) %>%
  filter(sum.count.per > 0.0000001) %>%
  select(-c(gr.zero, sum.count.per, sum.count)) 

fticr.final <- fticr.f %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  merge(new.id, by = "SampleID") %>%
  column_to_rownames(var = "New_id") %>%
  select(-SampleID)

fticr_bin <- fticr.final > 0  # Converts abundance to presence/absence

#metaT.log2 <- log2(metaT +1)

###determine specific compounds driving clustering 

km <- kmeans(fticr_bin, centers=3) ## Creates three clusters using kmeans

## Runs the combination analysis using IndVal.g as statistic
pt <- multipatt(fticr.final, km$cluster, control = how(nperm=999)) 

## Lists those species with significant association to one combination
summary(pt) 

# Extract results
results_df <- as.data.frame(pt$sign)
results_df$Feature <- rownames(results_df)
results_df$stat <- pt$stat[,1]
results_df$p.value <- pt$sign$p.value

#upset plot
results_df_upset <- results_df %>%
  select(-index, -p.value, -Feature) # %>%
  #rename(rhizo_normal = s.1, bulk = s.2, detritus_drought = s.3) 

upset(results_df_upset, order.by = "freq", mainbar.y.label = "Feature Intersections")
ggsave("./Figures/nmds/fticr/upset_plot_km.png", height = 2, width = 2, units = "in", dpi = 300)

#plot classes
all <- results_df %>%
  filter(s.1 == 1 & s.2 == 1 & s.3 == 1) %>%
  mutate(category1 = "all") %>%
  mutate(category2 = "1,2,3")

rhizo.normal.only <- results_df %>%
  filter(s.1 == 1 & s.2 == 0 & s.3 == 0) %>%
  mutate(category1 = "rhizo_normal_only") %>%
  mutate(category2 = "1")

bulk.only <- results_df %>%
  filter(s.2 == 1 & s.1 == 0 & s.3 == 0) %>%
  mutate(category1 = "bulk_only") %>%
  mutate(category2 = "2")

detritus.drought.only <- results_df %>%
  filter(s.1 == 0 & s.2 == 0 & s.3 ==1) %>%
  mutate(category1 = "detritus_drought_only") %>%
  mutate(category2 = "3")

rhizo.normal.detritus.drought <- results_df %>%
  filter(s.1 == 1 & s.2 == 0 & s.3 == 1) %>%
  mutate(category1 = "rhizo_normal_detritus_drought") %>%
  mutate(category2 = "1, 3")

bulk.detritus.drought <- results_df %>%
  filter(s.1 == 0 & s.2 == 1 & s.3 == 1) %>%
  mutate(category1 = "bulk_detritus_drought") %>%
  mutate(category2 = "2, 3")

#make pie 
comb <- rbind(all, rhizo.normal.only, bulk.only, detritus.drought.only, rhizo.normal.detritus.drought, bulk.detritus.drought)

drivers.pie <- comb %>%
  separate(Feature, c("Mass", "Formula", "Class"), sep = "_") %>%
  group_by(category2, Class) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = "", y = count, fill = Class)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c(col_list2)) +
  facet_wrap(~ category2, scales = "free") +
  theme_void()
drivers.pie
ggsave("./Figures/nmds/fticr/pie_drivers_all_compounds_kmeans_3_groups.png", height = 4, width = 4, units = "in", dpi = 300)


drivers.pie.unique <- comb %>%
  separate(Feature, c("Mass", "Formula", "Class"), sep = "_") %>%
  filter(category2 == "1" | category2 == "2" | category2 == "3") %>%
  group_by(category2, Class) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = "", y = count, fill = Class)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c(col_list2)) +
  facet_wrap(~ category2, scales = "free") +
  theme_void()
drivers.pie.unique
ggsave("./Figures/nmds/fticr/pie_drivers_all_compounds_kmeans_3_groups_unique.png", height = 4, width = 4, units = "in", dpi = 300)

drivers.pie.unique.no.other <- comb %>%
  separate(Feature, c("Mass", "Formula", "Class"), sep = "_") %>%
  filter(category2 == "1" | category2 == "2" | category2 == "3") %>%
  filter(Class != "Other") %>%
  group_by(category2, Class) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = "", y = count, fill = Class)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c(col_list2)) +
  facet_wrap(~ category2, scales = "free") +
  theme_void()
drivers.pie.unique.no.other
ggsave("./Figures/nmds/fticr/pie_drivers_all_compounds_kmeans_3_groups_unique_no_other.png", height = 4, width = 4, units = "in", dpi = 300)


drivers.bar <- comb %>%
  separate(Feature, c("Mass", "Formula", "Class"), sep = "_") %>%
  group_by(category2, Class) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = category2, y = count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(col_list2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
drivers.bar
ggsave("./Figures/nmds/fticr/bar_drivers_all_compounds_kmeans_3_groups.png", height = 4, width = 4, units = "in", dpi = 300)

drivers.bar.unique <- comb %>%
  separate(Feature, c("Mass", "Formula", "Class"), sep = "_") %>%
  filter(category2 == "1" | category2 == "2" | category2 == "3") %>%
  group_by(category2, Class) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = category2, y = count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(col_list2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
drivers.bar.unique
ggsave("./Figures/nmds/fticr/bar_drivers_all_compounds_kmeans_3_groups_unique.png", height = 4, width = 4, units = "in", dpi = 300)


drivers.bar.unique.no.other <- comb %>%
  separate(Feature, c("Mass", "Formula", "Class"), sep = "_") %>%
  filter(category2 == "1" | category2 == "2" | category2 == "3") %>%
  filter(Class != "Other") %>%
  group_by(category2, Class) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = category2, y = count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(col_list2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
drivers.bar.unique.no.other
ggsave("./Figures/nmds/fticr/bar_drivers_all_compounds_kmeans_3_groups_unique_no_other.png", height = 4, width = 4, units = "in", dpi = 300)


#add cluster to metadata in order to plot kmeans in nmds plot

cluster_df <- as.data.frame(km$cluster) %>%
  rownames_to_column(var = "new_ID") 

colnames(cluster_df) <- c("new_ID", "cluster")

new.id.meta.km <- new.id.meta %>%
  rownames_to_column(var = "new_ID") %>%
  merge(cluster_df, by = "new_ID", all = TRUE)

#add cluster to fticr.final new_id
fticr.final.c <- fticr.final %>%
  rownames_to_column(var = "new_ID") %>%
  merge(cluster_df, by = "new_ID") %>%
  tidyr::unite(col = "new_ID_cluster", new_ID, cluster, sep = "_") %>%
  column_to_rownames(var = "new_ID_cluster")



########################################## nmds on metaT-self assembly#####################
# Calculate nmds
fticr_bin <- fticr.final.c > 0  # Converts abundance to presence/absence

set.seed(123)
nmds = metaMDS(fticr_bin, distance = "jaccard")
plot(nmds)

scores <- scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  separate(., SampleID, into = c("SampleID", "Zone", "Treatment", "Timepoint", "cluster"), sep = "_", remove = T) %>%
  mutate(Zone = str_replace(Zone, "Rhizosphere", "Rhizo"))


#determine fit of features (how do features drive patterns, which are significant drivers?)

#fit <- envfit(nmds, fticr_bin, permutations = 1)

#plot(nmds)
#plot(fit, p.max = 0.05)  # Only show vectors with p ≤ 0.05


#determine environmetnal drivers (how do environmetnal variables drive patterns)

ef <- envfit(nmds ~ Zone + Treatment + Timepoint + cluster, data = new.id.meta.km,
             permutations = 999)
plot(nmds, type = "n")
points(nmds, pch = 21, bg = "lightgreen")
plot(ef, col = "blue", lwd = 2)  # draws arrows for each env var

#Goodness of fit:
#  r2 Pr(>r)    
#r2 Pr(>r)    
#Zone      0.5929  0.001 ***
#  Treatment 0.0117  0.309  
#Timepoint 0.0169  0.490  

#write.csv(data.scores,file=paste0("./figures/nmds/nmds_individual_coordinates.csv"),row.names=TRUE)

### nmds plot (timepoint and moisture)
make_nmds_plot(data.scores, Group1 =Zone, Group2 =Treatment)

filename <- paste0("./Figures/nmds/fticr/zone_treatment_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

### nmds plot (timepoint and moisture) with ellipses
make_nmds_plot_ellipse(data.scores, Group1 =Zone, Group2 =Treatment)

filename <- paste0("./Figures/nmds/fticr/zone_treatment_ellipse_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)


### nmds plot (timepoint and moisture)
make_nmds_plot(data.scores, Group1 =Treatment, Group2 =Zone)


filename <- paste0("./Figures/nmds/fticr/treatment_zone_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

### nmds plot (zone and timepoint)
make_nmds_plot(data.scores, Group1 =Zone, Group2 =Timepoint)

filename <- paste0("./Figures/nmds/fticr/zone_timepoint_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

### nmds plot (cluster and zone)
make_nmds_plot(data.scores, Group1 =cluster, Group2 =Zone)

filename <- paste0("./Figures/nmds/fticr/cluster_zone_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

### nmds plot (cluster and treatment)
make_nmds_plot(data.scores, Group1 =cluster, Group2 =Treatment)

filename <- paste0("./Figures/nmds/fticr/cluster_treatment_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

### nmds plot (cluster and timepoint)
make_nmds_plot(data.scores, Group1 =cluster, Group2 =Timepoint)

filename <- paste0("./Figures/nmds/fticr/cluster_timepoint_all_compounds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)



####statistical tests#####
##metaT###
#df.df <- as.data.frame(df)
df.df <- as.data.frame(fticr.final)

df.df <- df.df[sapply(df.df, is.numeric)]


df.stat <- df.df %>% 
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c("SampleID", "Zone", "Treatment", "Timepoint"), sep = "_", remove = T) 
  


####permanova on jaccard distance matrix####
##all time points##
dist <- vegdist(df.df > 0, method = "jaccard")  # binary matrix

adonis2(dist ~ Treatment*Zone,  df.stat)
#adonis2(formula = dist ~ Treatment * Zone, data = df.stat)
#Df SumOfSqs      R2      F Pr(>F)    
#Model     7   3.0108 0.46145 10.894  0.001 ***
#  Residual 89   3.5138 0.53855                  
#Total    96   6.5246 1.00000   

adonis2(dist ~ Zone,  df.stat)
#adonis2(formula = dist ~ Zone, data = df.stat)
#Df SumOfSqs      R2      F Pr(>F)    
#Model     3   2.6198 0.40152 20.798  0.001 ***
#  Residual 93   3.9048 0.59848                  
#Total    96   6.5246 1.00000 

     
adonis2(dist ~ Treatment,  df.stat)
#adonis2(formula = dist ~ Treatment, data = df.stat)
#Df SumOfSqs      R2      F Pr(>F)
#Model     1   0.1280 0.01962 1.9013  0.102
#Residual 95   6.3966 0.98038              
#Total    96   6.5246 1.00000   

adonis2(dist ~ Timepoint, df.stat)
#adonis2(formula = dist ~ Timepoint, data = df.stat)
#Df SumOfSqs      R2      F Pr(>F)
#Model     2   0.1534 0.02351 1.1317  0.298
#Residual 94   6.3712 0.97649              
#Total    96   6.5246 1.00000       

adonis2(dist ~ Treatment*Zone*Timepoint, df.stat)
#adonis2(formula = dist ~ Treatment * Zone * Timepoint, data = df.stat)
#Df SumOfSqs     R2      F Pr(>F)    
#Model    23   4.2481 0.6511 5.9229  0.001 ***
#  Residual 73   2.2765 0.3489  







  






