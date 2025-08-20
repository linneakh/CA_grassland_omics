# ***************************************************************
# 
# MetaboDirect
# Data Exploration step
# MetaboDirect version 1.0
# by Christian Ayala
#     based on scripts and functions by Nathalia Graf Grachet
# Licensed under the MIT license. See LICENSE.md file.
#
# Modified by Linnea Honeker
# ***************************************************************

# Loading libraries ----

library(ggvenn)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(vegan)
library(UpSetR)
library(lmerTest)
library(lme4)
library(emmeans)
library(afex)
library(tidyverse)


# Defining variables ----

output_dir <- './output/H2O-FTICR/'
figure_dir <- './figures/Fig3-FigS3-H2O-fticr'


# Loading custom functions ----
source("./scripts/metabodirect-custom_functions.R")

class_colors <- get_palette(palette = 'Set3', k = 9)
names(class_colors) <- c(classification$Class, 'Other')

colors.class = c("black", "#6D6E71","#C4A484",
                 "chartreuse3", "darkseagreen1",
                 "#F7941D", "yellow","darkorchid", "plum2",
                 "darkred","green4")

my_colors_moisture <- c("brown", "darkgreen")

# Import data ----

## Defining paths
my_data.file <- file.path(output_dir, 'metabodirect', 'Report_processed_noNorm_MolecFormulas.csv')
my_metadata.file <- file.path('./data/FTICR/metadata-H20-no-extra.csv')
my_classcomp.file <- file.path(output_dir,'metabodirect', 'class_composition.csv')

## Loading tables
df <- read_csv(my_data.file)
metadata <- read_csv(my_metadata.file)
class_comp <- read_csv(my_classcomp.file)

# Reformat data files ----

## Change all metadata columns to factors to avoid problems at plotting

metadata <- metadata %>% 
  mutate(across(!SampleID, as.factor)) %>%
  filter(!grepl("DI", SampleID)) %>%
  separate(SampleID, c("SampleInfo"), sep = "_", remove = FALSE) %>%
  separate(SampleInfo, c("Time_moisture", "Plot"), sep = "t") %>%
  mutate(Plot = str_remove(Plot, "a")) %>%
  mutate(Plot = str_remove(Plot, "b")) %>%
  dplyr::select(-Time_moisture) 
  
## Intensity data file
df_longer_orig <- df %>%
  select(-contains("DI")) %>%
  pivot_longer(metadata$SampleID, names_to = 'SampleID', values_to = 'Intensity') %>% 
  filter(Intensity > 0) %>% 
  left_join(metadata, by = 'SampleID')

## Class composition
class_comp <- class_comp %>%
  pivot_longer(!SampleID, names_to = 'Class', values_to = 'Count') %>% 
  left_join(metadata, by = 'SampleID') %>%
  filter(!grepl("DI", SampleID)) 
  

# Plotting ----

#boxplots of nosc and gfe across all compounds

my_nosc_plot <- df_longer_orig %>%
  group_by(Plot, Timepoint, Moisture) %>%
  summarise(NOSC = mean(NOSC)) %>%
  ggplot(aes(x = Timepoint, y = NOSC, fill = Moisture)) +

  geom_boxplot() +
  scale_fill_manual(values = c(my_colors_moisture)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),   # ← removes major grid lines
        panel.grid.minor = element_blank() )# ← removes minor  )
my_nosc_plot
ggsave(paste(figure_dir, "/NOSC_moisture_time.png", sep = ""), dpi = 300, height = 4, width = 5, unit = "in")
ggsave(paste(figure_dir, "/NOSC_moisture_time.pdf", sep = ""), dpi = 300, height = 4, width = 5, unit = "in")


my_gfe_plot <- df_longer_orig %>%
  group_by(Plot, Timepoint, Moisture) %>%
  summarise(GFE = mean(GFE)) %>%
  ggplot(aes(x = Timepoint, y = GFE, fill = Moisture)) +
  geom_boxplot() +
  scale_fill_manual(values = c(my_colors_moisture)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
my_gfe_plot
ggsave(paste(figure_dir, "/GFE_moisture_time.png", sep = ""), dpi = 300, height = 4, width = 5, unit = "in")
ggsave(paste(figure_dir, "/GFE_moisture_time.pdf", sep = ""), dpi = 300, height = 4, width = 5, unit = "in")



###test if nosc and gfe are significantly different
#df_longer_orig <- df_longer
df_longer <- df_longer_orig %>%
  group_by(Plot, Timepoint, Moisture) %>%
  summarise(NOSC = mean(NOSC))

#lmer
m <- lmer(NOSC ~ Moisture * Timepoint + (1|Plot), data = df_longer)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                      Sum Sq    Mean Sq NumDF DenDF F value   Pr(>F)   
#Moisture           0.0007209 0.00072090     1     4  1.8017 0.250628   
#Timepoint          0.0108236 0.00216472     5    20  5.4102 0.002629 **
#Moisture:Timepoint 0.0005442 0.00010883     5    20  0.2720 0.923087 

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

# imepoint = 0:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.01162 0.0196 16.4   0.594  0.5609
# 
# Timepoint = 3:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.01600 0.0196 16.4   0.817  0.4255
# 
# Timepoint = 24:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.01199 0.0196 16.4   0.612  0.5486
# 
# Timepoint = 48:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.00962 0.0196 16.4   0.491  0.6297
# 
# Timepoint = 72:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.03254 0.0196 16.4   1.662  0.1155
# 
# Timepoint = 168:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.02043 0.0196 16.4   1.044  0.3118


##GFE
df_longer <- df_longer_orig %>%
  group_by(Plot, Timepoint, Moisture) %>%
  summarise(GFE = mean(GFE))

#lme
m <- lmer(GFE ~ Moisture * Timepoint + (1|Plot), data = df_longer)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
#Moisture           0.5855 0.58555     1     4  1.8017 0.250628   
#Timepoint          8.7915 1.75830     5    20  5.4102 0.002629 **
#Moisture:Timepoint 0.4420 0.08840     5    20  0.2720 0.923087  

emmeans(m, pairwise ~ Moisture | Timepoint)
pairs(em, adjust="tukey") #sig different at t=168 h



## Plot - Class comp bar ----

#plot class as boxplot
#create second dataframe to enforce axes
df2 <- data.frame(Timepoint = rep("0",18), 
                  Class = c(unique(class_comp$Class),unique(class_comp$Class)),
                  Count = c(2.5,4,38,52, 3, 0.3, 3.5, 18, 0.8,2.5,4,38,52, 3, 0.3, 3.5, 18, 0.8),
                  Moisture = c("50", "100"))

class_comp$Class <- factor(class_comp$Class, c("Amino sugar", "Carbohydrate", "Protein", "Lipid",
                                               "Cond. HC", "Unsat. HC", "Lignin", "Tannin", "Other"))
df2$Class <- factor(df2$Class, c("Amino sugar", "Carbohydrate", "Protein", "Lipid",
                                               "Cond. HC", "Unsat. HC", "Lignin", "Tannin", "Other"))



class_box <- class_comp %>%
  ggplot(aes(x=Timepoint, y = Count, fill = Moisture)) +
  geom_boxplot(position = "dodge") +
  geom_point(data = df2, aes(x = Timepoint, y = Count), colour = "white") +
  scale_fill_manual(values = my_colors_moisture) +
  labs(y="Relative abundance", x = "Timepoint (h)") +
  facet_wrap(~Class , scales = "free_y", nrow = 9) +
  theme_bw () +
  theme(text = element_text(size = 16),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
class_box
filename <- file.path(figure_dir, 'FigS3-Composition_by_class_boxplot.png')
ggsave(filename, class_box, dpi = 300, width = 6, height = 15, units = "in")
filename <- file.path(figure_dir, 'FigS3-Composition_by_class_boxplot.pdf')
ggsave(filename, class_box, dpi = 300, width = 6, height = 15, units = "in")



#plot pi chart of all compounds
pie <- class_comp %>%
  group_by(Class) %>%
  summarise(mean = mean(Count)) %>%
  mutate(mean = round(mean, 2)) %>%
  drop_na() %>%
  ggplot(aes(x= "", y = mean, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(colors.class)) +
  coord_polar(theta = "y") +
  labs(title = "Class composition", fill = "Class") +
  theme_void() + 
  theme(legend.position = "bottom")
  
pie
filename <- file.path(figure_dir, "FigS3-Composition_by_class_piechart.png")
ggsave(filename, pie, dpi = 300, width = 4, height = 4, units = "in")
filename <- file.path(figure_dir, "FigS3-Composition_by_class_piechart.pdf")
ggsave(filename, pie, dpi = 300, width = 4, height = 4, units = "in")

#summarize mean counts of each class to label pie chart
class.summary <- class_comp %>%
  group_by(Class) %>%
  summarise(mean = mean(Count),
            sd = sd(Count))
class.summary
#Class          mean     sd
#<chr>         <dbl>  <dbl>
#1 Amino sugar   0.821 0.361 
#2 Carbohydrate  1.05  0.627 
#3 Protein       2.08  0.441 
#4 Lipid         1.02  0.272 
#5 Cond. HC     32.3   1.73  
#6 Unsat. HC     0.225 0.104 
#7 Lignin       47.9   1.08  
#8 Tannin       14.5   0.809 
#9 Other         0.195 0.0632


#class statistics
#amino sugars
AS <- class_comp %>%
  filter(Class == "Amino sugar")


#Amino sugars with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = AS)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)   
#Moisture           0.06312 0.06312     1     4  0.8913 0.39859   
#Timepoint          2.08165 0.41633     5    20  5.8784 0.00169 **
#Moisture:Timepoint 0.07002 0.01400     5    20  0.1977 0.95968   

em <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey") #no sig diff

em <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")
#100
# Timepoint0 - Timepoint48     0.7300 0.217 20   3.360  0.0319
# Timepoint0 - Timepoint72     0.6800 0.217 20   3.129  0.0514
# Timepoint0 - Timepoint168    0.7000 0.217 20   3.221  0.0426

#carbs
carb <- class_comp %>%
  filter(Class == "Carbohydrate")

m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = carb)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Moisture            0.0126 0.01264     1     4  0.2252    0.6599    
#Timepoint          11.1568 2.23135     5    20 39.7541 1.016e-09 ***
#Moisture:Timepoint  0.0146 0.00293     5    20  0.0521    0.9980    

em <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey") #no sig diff

em <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey") 
# 100
# Timepoint0 - Timepoint24     0.9967 0.193 20   5.152  0.0006
# Timepoint0 - Timepoint48     1.2700 0.193 20   6.565  <.0001
# Timepoint0 - Timepoint72     1.3333 0.193 20   6.893  <.0001
# Timepoint0 - Timepoint168    1.4767 0.193 20   7.634  <.0001
# Timepoint3 - Timepoint24     0.7567 0.193 20   3.912  0.0097
# Timepoint3 - Timepoint48     1.0300 0.193 20   5.325  0.0004
# Timepoint3 - Timepoint72     1.0933 0.193 20   5.652  0.0002
# Timepoint3 - Timepoint168    1.2367 0.193 20   6.393  <.0001

# 50
# Timepoint0 - Timepoint24     0.9967 0.193 20   5.152  0.0006
# Timepoint0 - Timepoint48     1.1733 0.193 20   6.066  0.0001
# Timepoint0 - Timepoint72     1.3433 0.193 20   6.944  <.0001
# Timepoint0 - Timepoint168    1.4733 0.193 20   7.616  <.0001
# Timepoint3 - Timepoint24     0.7267 0.193 20   3.757  0.0136
# Timepoint3 - Timepoint48     0.9033 0.193 20   4.670  0.0018
# Timepoint3 - Timepoint72     1.0733 0.193 20   5.549  0.0003
# Timepoint3 - Timepoint168    1.2033 0.193 20   6.221  0.0001

#CH
CH <- class_comp %>%
  filter(Class == "Cond. HC")


#CH with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = CH)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Moisture            0.013  0.0126     1     4  0.0090 0.9289331    
#Timepoint          55.604 11.1209     5    20  7.9434 0.0002923 ***
#Moisture:Timepoint  0.788  0.1576     5    20  0.1126 0.9882054    


# Post-hoc pairwise contrasts between each timeoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")
# 50
# Timepoint0 - Timepoint168    -3.557 0.966 20  -3.681  0.0160
# 100
# Timepoint0 - Timepoint48     -3.077 0.966 20  -3.185  0.0459
# Timepoint0 - Timepoint168    -3.237 0.966 20  -3.350  0.0325

#lignin (nothing significant)
lig <- class_comp %>%
  filter(Class == "Lignin")

#lignin with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = lig)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
#Moisture           0.00710 0.00710     1     4  0.0084 0.9315
#Timepoint          2.73126 0.54625     5    20  0.6437 0.6693
#Moisture:Timepoint 0.42299 0.08460     5    20  0.0997 0.9911


# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")



#Protein
pro <- class_comp %>%
  filter(Class == "Protein")


#protein with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = pro)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
#Moisture           0.24643 0.24643     1     4  2.9400 0.161557   
#Timepoint          2.49076 0.49815     5    20  5.9433 0.001592 **
#Moisture:Timepoint 0.21242 0.04248     5    20  0.5069 0.767664   



# Post-hoc pairwise contrasts between timepoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")
# 100
# Timepoint0 - Timepoint48     0.7600 0.236 20   3.215  0.0431

#unsaturated HC
uHC <- class_comp %>%
  filter(Class == "Unsat. HC")


#unsaturated hydrocarbon with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = uHC)
anova(m)
#Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq  Mean Sq NumDF DenDF F value   Pr(>F)    
#Moisture           0.001225 0.001225     1    24  0.2420   0.6272    
#Timepoint          0.218381 0.043676     5    24  8.6297 8.65e-05 ***
#Moisture:Timepoint 0.036225 0.007245     5    24  1.4315   0.2490  

em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")# Timepoint = 168:
#   contrast                 estimate     SE df t.ratio p.value
# Moisture50 - Moisture100  -0.1233 0.0581 24  -2.123  0.0442

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")
# 100
# Timepoint0 - Timepoint168  -0.33333 0.0581 20  -5.739  0.0002
# Timepoint3 - Timepoint168  -0.21000 0.0581 20  -3.615  0.0185
# Timepoint48 - Timepoint168 -0.26000 0.0581 20  -4.476  0.0027


#lipid
lip <- class_comp %>%
  filter(Class == "Lipid")


#lipids without t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = lip)
anova(m)
                         
#Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
#Moisture           0.09401 0.094008     1     4  3.3892 0.1394415    
#Timepoint          0.91618 0.183236     5    20  6.6060 0.0008806 ***
#Moisture:Timepoint 0.13265 0.026529     5    20  0.9564 0.4672782 

em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 72:
#   contrast                 estimate    SE   df t.ratio p.value
# Moisture50 - Moisture100  -0.4367 0.174 13.7  -2.513  0.0252

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")
# 50
# Timepoint0 - Timepoint72     0.4367 0.136 20   3.211  0.0435
# Timepoint0 - Timepoint168    0.5200 0.136 20   3.824  0.0117
# 100
# Timepoint3 - Timepoint168    0.4600 0.136 20   3.383  0.0304

#tannin
tan <- class_comp %>%
  filter(Class == "Tannin")



#tannin with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = tan)
anova(m)
# #Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# Moisture           1.52663 1.52663     1     4  4.3343 0.1058
# Timepoint          0.33822 0.06764     5    20  0.1921 0.9621
# Moisture:Timepoint 0.60966 0.12193     5    20  0.3462 0.8786 

em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")

# Comparisons ----
#new_colors <- c("purple", "seagreen", "grey80")
new_colors <- c("darkgreen", "brown", "grey80")


#create new column called "presence" to indicate if compound found in both or are unique
unique_compounds <- df_longer_orig %>%
  distinct(Mass, Plot, Moisture) %>%
  mutate(Moisture = paste("x", Moisture, sep = "")) %>%
  group_by(Mass, Moisture) %>%
  summarise(n_treat = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = Moisture,
    values_from = n_treat,
    values_fill = 0) %>%
  mutate(present_50 = case_when(
    x50 > 1 ~ 'present',
    TRUE ~ 'not_present'
  )) %>%
  mutate(present_100 = case_when(
    x100 > 1 ~ 'present',
    TRUE ~ 'not_present'
  )) %>%
  mutate(Presence = case_when(
    (present_50 == "present" &
       present_100 == "not_present") ~ "u50",
    (present_50 == "not_present" &
       present_100 == "present") ~ "u100",
    TRUE ~ 'shared'
  )) 

df_group_moisture <- df_longer_orig %>%
  select(Mass, HC, OC) %>% 
  distinct() %>% 
  merge(unique_compounds, by = "Mass") %>%
  select(Mass, HC, OC, Presence) %>%
  mutate(Presence = factor(Presence, levels = c("u100", "u50", "shared")))

#group by presence and summarise NOSC, GFE, and AI_mod
df_presence_table_summary_moisture <- df_group_moisture %>%
  merge(df_longer_orig, by = c("Mass", "HC", "OC")) %>%
  select(-Timepoint, SampleID, Moisture, Plot) %>%
  distinct() %>%
  group_by(Presence) %>%
  summarise(NOSC = mean(NOSC),
            GFE = mean(GFE),
            AI_mod = mean(AI_mod)) #%>%

#group by presence and determine total counts of each compound class in each category
df_presence_table_summary_moisture_class <- df_group_moisture %>%
  merge(df_longer_orig, by = c("Mass", "HC", "OC")) %>%
  select(-Timepoint, -SampleID, -Moisture, -Plot, -Intensity) %>%
  distinct() %>%
  group_by(Presence, Class) %>%
  summarise(count = n()) 

#group by presence and plot to determine total counts of each compound class in each category per plot
df_presence_table_summary_moisture_class_plot <- df_group_moisture %>%
  merge(df_longer_orig, by = c("Mass", "HC", "OC")) %>%
  select(-Timepoint, -SampleID, -Moisture, -Intensity) %>%
  distinct() %>%
  group_by(Presence, Class, Plot) %>%
  summarise(count = n()) 

# Presence Class        count
# <chr>    <chr>        <int>
#   1 shared   Amino sugar     97
# 2 shared   Carbohydrate   167
# 3 shared   Cond. HC      2443
# 4 shared   Lignin        3795
# 5 shared   Lipid          128
# 6 shared   Other           23
# 7 shared   Protein        216
# 8 shared   Tannin        1215
# 9 shared   Unsat. HC       33
# 10 u100     Amino sugar     41
# 11 u100     Carbohydrate    23
# 12 u100     Cond. HC        36
# 13 u100     Lignin          91
# 14 u100     Lipid           27
# 15 u100     Protein         38
# 16 u100     Tannin          14
# 17 u100     Unsat. HC       10
# 18 u50      Amino sugar      6
# 19 u50      Carbohydrate    20
# 20 u50      Cond. HC       170
# 21 u50      Lignin         288
# 22 u50      Lipid           10
# 23 u50      Other            4
# 24 u50      Protein         13
# 25 u50      Tannin         131
# 26 u50      Unsat. HC        1


#merge df_group_moisture with df_longer_orig to add thermodynamics and class composition
df_group_moisture_timepoint <- df_group_moisture %>%
  merge(df_longer_orig, by = c("Mass", "HC", "OC")) 

#group by time, plot, class, and presence in order to determine mean thermodynamics per timepoint per plot per class
df_group_moisture_timepoint_class_stat <- df_group_moisture_timepoint %>%
  group_by(Plot, Class, Presence, Timepoint) %>%
  summarise(NOSC = mean(NOSC),
            GFE = mean(GFE),
            AI_mod = mean(AI_mod))

#group by time, plot, and presence in order to determine mean thermodynamics per timepoint per plot
df_group_moisture_timepoint_stat <- df_group_moisture_timepoint %>%
  group_by(Plot, Presence, Timepoint) %>%
  summarise(NOSC = mean(NOSC),
            GFE = mean(GFE),
            AI_mod = mean(AI_mod))


#vk plots
vk_plot <- plot_van_krevelen(df_group_moisture, color_by = 'Presence') +
  scale_color_manual(labels = c("100", "50", "shared"),
                     values = new_colors) +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),   # ← removes major grid lines
        panel.grid.minor = element_blank() )# ← removes minor  )
vk_plot
filename <- file.path(figure_dir, paste0('FigS3-vk_', g, '_all_data.png'))
ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)
filename <- file.path(figure_dir, paste0('FigS3_vk_', g, '_all_data.pdf'))
ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)



#plot vk_plot without shared points
df_group_unique <- df_group_moisture %>%
  filter(Presence != "shared")

vk_plot <- plot_van_krevelen(df_group_unique, color_by = 'Presence') +
  scale_color_manual(labels = c("100", "50", "shared"),
                     values = new_colors) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),   # ← removes major grid lines
        panel.grid.minor = element_blank() )# ← removes minor  )
vk_plot
filename <- file.path(figure_dir, paste0('Fig3-vk_', g, '_unique.png'))
ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)
filename <- file.path(figure_dir, paste0('Fig3-vk_', g, '_unique.pdf'))
ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)


#plot overage NOSC and GFE across shared and unique compounds at each timepoint
plot_nosc <- df_group_moisture_timepoint %>%
  ggplot(aes(x=Presence, y = NOSC, fill = Presence)) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  
  geom_boxplot(width = 0.2) +
  scale_fill_manual(labels = c("100", "50", "shared"),
                    values = c(new_colors)) +
  scale_x_discrete(labels = c("100", "50", "shared")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),   # ← removes major grid lines
        panel.grid.minor = element_blank() )# ← removes minor  )
plot_nosc
filename <- file.path(figure_dir, paste0('FigS4-NOSC-unique-violin.png'))
ggsave(filename, plot_nosc, dpi = 300, width = 4, height = 4)
filename <- file.path(figure_dir, paste0('FigS4-NOSC-unique-violin.pdf'))
ggsave(filename, plot_nosc, dpi = 300, width = 4, height = 4)


plot_gfe <- df_group_moisture_timepoint %>%
  ggplot(aes(x=Presence, y = GFE, fill = Presence)) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  
  geom_boxplot(width = 0.2) +
  scale_fill_manual(labels = c("100", "50", "shared"),
                    values = c(new_colors)) +
  scale_x_discrete(labels = c("100", "50", "shared")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),   # ← removes major grid lines
        panel.grid.minor = element_blank() )# ← removes minor  )
plot_gfe
filename <- file.path(figure_dir, paste0('FigS4-GFE-unique-violin.png'))
ggsave(filename, plot_gfe, dpi = 300, width = 4, height = 4)
filename <- file.path(figure_dir, paste0('FigS4-GFE-unique-violin.pdf'))
ggsave(filename, plot_gfe, dpi = 300, width = 4, height = 4)


#plot unique classes at each timepoint
plot_classes <- df_presence_table_summary_moisture_class %>%
  ggplot(aes(x=Presence, y = count, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(colors.class)) +
  scale_x_discrete(labels = c("100", "50", "shared")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),   # ← removes major grid lines
        panel.grid.minor = element_blank() )# ← removes minor  )
plot_classes
filename <- file.path(figure_dir, paste0('Fig3-class-unique.png'))
ggsave(filename, plot_classes, dpi = 300, width = 3, height = 4)
filename <- file.path(figure_dir, paste0('Fig3-class-unique.pdf'))
ggsave(filename, plot_classes, dpi = 300, width = 3, height = 4)


##lme test nosc
df_group_moisture_timepoint_unique <- df_group_moisture_timepoint %>%
  distinct(Mass, NOSC, GFE, Plot, Presence) %>%
  filter(Presence != "shared") 

m <- lmerTest::lmer(NOSC ~ Presence  + (1|Plot), data = df_group_moisture_timepoint_unique)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Presence  102.9   102.9     1 2375.9  340.43 < 2.2e-16 ***



##lme test gfe
m <- lmerTest::lmer(GFE ~ Presence  + (1|Plot), data = df_group_moisture_timepoint_unique)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Presence  83580   83580     1 2375.9  340.43 < 2.2e-16 ***


###g-test
library(DescTools)


# 3) Build the contingency table of Class × Moisture for those masses
tab_unique_thresh <- df_longer_orig %>%
  distinct(Mass, Moisture, Class) %>%
  # keep only Mass/Moisture combos that passed both filters
  merge(unique_compounds, by = "Mass") %>%
  distinct(Mass, Presence, Class) %>%
  filter(Presence != "shared") %>%
  dplyr::count(Presence, Class, name = "n") %>%
  pivot_wider(
    names_from  = Presence,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("Class") %>%
  as.matrix() 

# 4) Run the G‐test and look at standardized residuals
GTest(tab_unique_thresh)

# Log likelihood ratio (G-test) test of independence without correction
# 
# data:  tab_unique_thresh
# G = 229.01, X-squared df = 8, p-value < 2.2e-16

chi <- chisq.test(tab_unique_thresh, correct = FALSE)
chi$stdres
#                   u100       u50
# Amino sugar   8.709903 -8.709903
# Carbohydrate  3.382276 -3.382276
# Cond. HC     -4.555516  4.555516
# Lignin       -3.489152  3.489152
# Lipid         5.758241 -5.758241
# Protein       7.060103 -7.060103
# Tannin       -5.900343  5.900343
# Unsat. HC     4.396405 -4.396405
# Other        -1.322655  1.322655



# pairwise G-tests: join the Class for each Mass, keep only u50/u100
uc <- unique_compounds %>%
  inner_join(df_longer_orig %>% distinct(Mass, Class),
             by = "Mass") %>%
  filter(Presence %in% c("u50","u100"))

# totals per category
totals <- uc %>% count(Presence) %>% deframe()

# for each Class, build a 2×2 and run GTest()
g_results <- uc %>%
  count(Class, Presence) %>%
  pivot_wider(names_from = Presence,
              values_from = n,
              values_fill  = 0) %>%
  rowwise() %>%
  mutate(
    a = u50,
    b = u100,
    c = totals["u50"]  - a,
    d = totals["u100"] - b,
    mat = list(matrix(c(a,b,c,d), nrow = 2)),
    G     = GTest(mat)$statistic,
    p.val = GTest(mat)$p.value
  ) %>%
  ungroup() %>%# **add this** FDR step:
  mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
  
  select(Class, G, p.val, p.adj)


g_results


