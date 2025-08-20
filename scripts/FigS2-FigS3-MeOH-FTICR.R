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

output_dir <- './output/MeOH-FTICR/'
figure_dir <- './figures/FigS2-FigS3-MeOH-fticr/'


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
my_metadata.file <- file.path('./data/FTICR/metadata-MeOH-no-extra-no-out.csv')
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
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq    Mean Sq NumDF   DenDF F value Pr(>F)
# Moisture           0.00079178 0.00079178     1  4.1693  0.7528 0.4327
# Timepoint          0.00264694 0.00052939     5 18.1695  0.5034 0.7699
# Moisture:Timepoint 0.00254056 0.00050811     5 18.1695  0.4831 0.7843

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey") #no significant diff


##GFE
df_longer <- df_longer_orig %>%
  group_by(Plot, Timepoint, Moisture) %>%
  summarise(GFE = mean(GFE))

#lme
m <- lmer(GFE ~ Moisture * Timepoint + (1|Plot), data = df_longer)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# Moisture           0.64312 0.64312     1  4.1693  0.7528 0.4327
# Timepoint          2.14997 0.42999     5 18.1695  0.5034 0.7699
# Moisture:Timepoint 2.06357 0.41271     5 18.1695  0.4831 0.7843

emmeans(m, pairwise ~ Moisture | Timepoint)
pairs(em, adjust="tukey") #no significant diff



## Plot - Class comp bar ----

#plot class as boxplot
#create second dataframe to enforce axes
df2 <- data.frame(Timepoint = rep("0",18), 
                  Class = c(unique(class_comp$Class),unique(class_comp$Class)),
                  Count = c(15,20,23,45,10, 0.125, 20, 4, 0.5,15,20,23,45,10, 0.125, 20, 4, 0.5),
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
# Class          mean      sd
# <fct>         <dbl>   <dbl>
#   1 Amino sugar   8.68   0.982 
# 2 Carbohydrate 11.8    1.81  
# 3 Protein      15.2    1.48  
# 4 Lipid         5.71   0.920 
# 5 Cond. HC     17.0    2.04  
# 6 Unsat. HC     0.315  0.0790
# 7 Lignin       39.0    2.43  
# 8 Tannin        2.14   0.528 
# 9 Other        NA     NA     


#class statistics
#amino sugars
AS <- class_comp %>%
  filter(Class == "Amino sugar")


#Amino sugars with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = AS)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
# Moisture            0.5578 0.55775     1    22  0.9801 0.332937   
# Timepoint          13.4174 2.68349     5    22  4.7155 0.004455 **
# Moisture:Timepoint  4.5539 0.91078     5    22  1.6005 0.201537   

em <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey") 
# Timepoint = 0:
#   contrast                 estimate    SE df t.ratio p.value
# Moisture50 - Moisture100  -1.6033 0.616 22  -2.603  0.0162

em <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey") #no sig diff
# Timepoint0 - Timepoint3      2.1867 0.616 18.0   3.550  0.0235

#carbs
carb <- class_comp %>%
  filter(Class == "Carbohydrate")

m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = carb)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# Moisture            0.915  0.9152     1    22  0.5636 0.4607863    
# Timepoint          56.462 11.2925     5    22  6.9533 0.0004948 ***
# Moisture:Timepoint 13.576  2.7152     5    22  1.6719 0.1832572    

em <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey") #T0
# Timepoint = 0:
#   contrast                 estimate   SE df t.ratio p.value
# Moisture50 - Moisture100  -2.2167 1.04 22  -2.130  0.0446

em <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey") #T0
# Timepoint0 - Timepoint3      3.8267 1.04 18.0   3.678  0.0181

#CH
CH <- class_comp %>%
  filter(Class == "Cond. HC")


#CH with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = CH)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# Moisture            0.338  0.3378     1  4.6254  0.1221 0.74209  
# Timepoint          40.270  8.0540     5 18.6126  2.9114 0.04143 *
# Moisture:Timepoint 27.054  5.4108     5 18.6126  1.9559 0.13300  


# Post-hoc pairwise contrasts between each timeoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 0:
#   contrast                 estimate   SE   df t.ratio p.value
# Moisture50 - Moisture100    3.093 1.44 20.8   2.145  0.0440

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")
# 100
#Timepoint0 - Timepoint48     -4.910 1.36 18.0  -3.616  0.0206
# Timepoint0 - Timepoint72     -4.670 1.36 18.0  -3.439  0.0296
# Timepoint0 - Timepoint168    -4.507 1.36 18.0  -3.319  0.0378

#lignin (nothing significant)
lig <- class_comp %>%
  filter(Class == "Lignin")

#lignin with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = lig)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# Moisture            0.3156  0.3156     1  4.5964  0.0486 0.8349
# Timepoint          26.2772  5.2554     5 18.6500  0.8100 0.5571
# Moisture:Timepoint 14.9084  2.9817     5 18.6500  0.4596 0.8012


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
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# Moisture           0.0005 0.00053     1  4.3031  0.0003 0.9866
# Timepoint          4.1836 0.83673     5 18.2880  0.4951 0.7758
# Moisture:Timepoint 7.3061 1.46122     5 18.2880  0.8646 0.5232 



# Post-hoc pairwise contrasts between timepoint
em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")

#unsaturated HC
uHC <- class_comp %>%
  filter(Class == "Unsat. HC")


#unsaturated hydrocarbon with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = uHC)
anova(m)

#Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq   Mean Sq NumDF   DenDF F value Pr(>F)
# Moisture           0.008632 0.0086324     1  4.6378  1.5530 0.2719
# Timepoint          0.014888 0.0029776     5 18.6144  0.5357 0.7467
# Moisture:Timepoint 0.034475 0.0068949     5 18.6144  1.2405 0.3301

em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")# Timepoint = 168:
# Timepoint = 168:
#   contrast                 estimate     SE   df t.ratio p.value
# Moisture50 - Moisture100  0.16667 0.0652 20.4   2.555  0.0187

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")

#lipid
lip <- class_comp %>%
  filter(Class == "Lipid")


#lipids without t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = lip)
anova(m)
                         
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# Moisture           0.0017 0.00170     1  4.242  0.0025 0.9623
# Timepoint          5.1909 1.03818     5 18.238  1.5349 0.2281
# Moisture:Timepoint 0.9804 0.19607     5 18.238  0.2899 0.9125

em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

em    <- emmeans(m, ~ Timepoint | Moisture)
pairs(em, adjust="tukey")



#tannin
tan <- class_comp %>%
  filter(Class == "Tannin")



#tannin with t0
m <- lmer(Count ~ Moisture * Timepoint + (1|Plot), data = tan)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# Moisture           0.51118 0.51118     1  3.9607  2.8464 0.1676
# Timepoint          1.55233 0.31047     5 17.9734  1.7288 0.1791
# Moisture:Timepoint 0.58080 0.11616     5 17.9734  0.6468 0.6675

em    <- emmeans(m, ~ Moisture | Timepoint)
pairs(em, adjust="tukey")

# imepoint = 24:
#   contrast                 estimate    SE   df t.ratio p.value
# Moisture50 - Moisture100    0.967 0.402 16.8   2.403  0.0281


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
filename <- file.path(figure_dir, paste0('FigS2-vk_', g, '_unique.png'))
ggsave(filename, vk_plot, dpi = 300, width = 8, height = 8)
filename <- file.path(figure_dir, paste0('FigS2-vk_', g, '_unique.pdf'))
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
filename <- file.path(figure_dir, paste0('FigS2-NOSC-unique-violin.png'))
ggsave(filename, plot_nosc, dpi = 300, width = 4, height = 4)
filename <- file.path(figure_dir, paste0('FigS2-NOSC-unique-violin.pdf'))
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
filename <- file.path(figure_dir, paste0('FigS2-GFE-unique-violin.png'))
ggsave(filename, plot_gfe, dpi = 300, width = 4, height = 4)
filename <- file.path(figure_dir, paste0('FigS2-GFE-unique-violin.pdf'))
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
filename <- file.path(figure_dir, paste0('FigS2-class-unique.png'))
ggsave(filename, plot_classes, dpi = 300, width = 3, height = 4)
filename <- file.path(figure_dir, paste0('FigS2-class-unique.pdf'))
ggsave(filename, plot_classes, dpi = 300, width = 3, height = 4)


##lme test nosc
df_group_moisture_timepoint_unique <- df_group_moisture_timepoint %>%
  distinct(Mass, NOSC, GFE, Plot, Presence) %>%
  filter(Presence != "shared") 

m <- lmerTest::lmer(NOSC ~ Presence  + (1|Plot), data = df_group_moisture_timepoint_unique)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Presence 7.8342  7.8342     1 2059.5   26.72 2.578e-07 ***



##lme test gfe
m <- lmerTest::lmer(GFE ~ Presence  + (1|Plot), data = df_group_moisture_timepoint_unique)
anova(m)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Presence 6363.4  6363.4     1 2059.5   26.72 2.578e-07 ***


###g-test
library(DescTools)


# 3) Build the contingency table of Class × Moisture for those masses
tab_unique_thresh <- df_longer_orig %>%
  distinct(Mass, Moisture, Class) %>%
  # keep only Mass/Moisture combos that passed both filters
  merge(unique_compounds, by = "Mass") %>%
  distinct(Mass, Presence, Class) %>%
  dplyr::filter(Presence != "shared") %>%
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
# G = 30.362, X-squared df = 8, p-value = 0.0001824

chi <- chisq.test(tab_unique_thresh, correct = FALSE)
chi$stdres
# u100        u50
# Amino sugar   1.6845314 -1.6845314
# Carbohydrate  4.1704501 -4.1704501
# Cond. HC     -2.3241605  2.3241605
# Lignin       -0.2035640  0.2035640
# Lipid        -1.3675079  1.3675079
# Protein      -0.1167391  0.1167391
# Tannin       -1.7268813  1.7268813
# Unsat. HC     1.0424522 -1.0424522
# Other        -1.2739422  1.2739422



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




