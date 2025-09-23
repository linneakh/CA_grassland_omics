### FUNCTIONS FOR nmds

# Helper: split SampleID into factors used above
split_sample_id <- function(df, id_col = "SampleID") {
  df %>% tidyr::separate({{id_col}}, into = c("plot", "moisture", "timepoint"), sep = "_", remove = FALSE) %>%
    mutate(
      timepoint = factor(timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168")),
      moisture  = factor(moisture,  levels = c("x50", "x100"))
    )
}


make_nmds_plot_time <- function(nmds_matrix) {
  nmds_plot <-  ggplot(mapping = aes(x, y)) +
    geom_point(data=nmds_matrix, aes(x=NMDS1, y=NMDS2, col=timepoint, shape=moisture), 
               size=size/2, show.legend = TRUE) +
    theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
    scale_color_viridis(discrete = TRUE, direction = -1, labels = c("0h", "3h", "24h", "48h", "72h", "168h")) +
    scale_shape_manual(values= list_of_shapes, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
    
    theme( text = element_text(size=size+5),
           legend.text = element_text(size=size+5, face="bold"),
           legend.title = element_blank(),
           #legend.key.size = unit(0.6, "cm"),
           #legend.key.width = unit(0.6,"cm"),
           legend.position = "bottom",
           axis.title.x = element_text(size=size+3,face="bold"),
           axis.title.y = element_text(size=size+3,face="bold"),
           plot.title = element_text(size=size+3,face="bold"),
           panel.grid.major = element_blank(),   # ← removes major grid lines
           panel.grid.minor = element_blank() )# ← removes minor  
  return(nmds_plot)
  
}

make_nmds_plot_time_centroid <- function(nmds_object, nmds_matrix) {
  # calculate fit of environmental data (moisture and timepoint)
  fit <- envfit(nmds_object, nmds_matrix, permutations = 999)
  factors <- as.data.frame(vegan::scores(fit, display = "factors")) %>%
    filter(grepl("timepoint", rownames(.)))
  
  #add significant factor (moisture) to plot as psuedo vectors
  nmds_plot_vector <- nmds_plot + geom_segment(data = factors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                                               arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1) +
    geom_text(data = factors, aes(x = NMDS1, y = NMDS2, label = rownames(factors)),
              color = "red", vjust = -1, fontface = "italic") 
  
  #add significant factor (moisture) to plot as centroids
  nmds_plot_centroid <- nmds_plot + geom_point(data = factors, aes(x = NMDS1, y = NMDS2), color = "red") +
    geom_text(data = factors, aes(x = NMDS1, y = NMDS2, label = rownames(factors)),
              color = "red", vjust = -1, fontface = "italic") 
  print(fit)
  return(nmds_plot_centroid)
}

make_nmds_plot_moisture <- function(nmds_matrix) {
  nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=nmds_matrix, aes(x=NMDS1, y=NMDS2, col=moisture), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = c(col_list_moisture), breaks = c("x50", "x100"),labels =c("50%", "100%")) +
theme( text = element_text(size=size+5),
       legend.text = element_text(size=size+5, face="bold"),
       legend.title = element_blank(),
       legend.position = "bottom",
       axis.title.x = element_text(size=size+3,face="bold"),
       axis.title.y = element_text(size=size+3,face="bold"),
       plot.title = element_text(size=size+3,face="bold"),
       panel.grid.major = element_blank(),   # ← removes major grid lines
       panel.grid.minor = element_blank() )# ← removes minor  
return(nmds_plot)
}
  
make_nmds_plot_moisture_centroid <- function(nmds_object, nmds_matrix) {
  # calculate fit of environmental data (moisture and timepoint)
  fit <- envfit(nmds_object, nmds_matrix, permutations = 999)
  factors <- as.data.frame(vegan::scores(fit, display = "factors")) %>%
    filter(grepl("moisture", rownames(.)))
  
  
  #add significant factor (moisture) to plot as psuedo vectors
  nmds_plot_vector <- nmds_plot + geom_segment(data = factors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                           arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1) +
    geom_text(data = factors, aes(x = NMDS1, y = NMDS2, label = rownames(factors)),
              color = "red", vjust = -1, fontface = "italic") 
  
  #add significant factor (moisture) to plot as centroids
  nmds_plot_centroid <- nmds_plot + geom_point(data = factors, aes(x = NMDS1, y = NMDS2), color = "red") +
    geom_text(data = factors, aes(x = NMDS1, y = NMDS2, label = rownames(factors)),
              color = "red", vjust = -1, fontface = "italic") 
  print(fit)
  return(nmds_plot_centroid)
}

make_nmds_plot_plot <- function(nmds_matrix) {
  nmds_plot <-  ggplot(mapping = aes(x, y)) +
    geom_point(data=nmds_matrix, aes(x=NMDS1, y=NMDS2, col=plot), 
               size=size/2, show.legend = TRUE) +
    theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
    scale_color_manual(values = c(col_list_subject)) +
    theme( text = element_text(size=size+5),
           legend.text = element_text(size=size+5, face="bold"),
           legend.title = element_blank(),
           legend.position = "bottom",
           axis.title.x = element_text(size=size+3,face="bold"),
           axis.title.y = element_text(size=size+3,face="bold"),
           plot.title = element_text(size=size+3,face="bold"),
           panel.grid.major = element_blank(),   # ← removes major grid lines
           panel.grid.minor = element_blank() )# ← removes minor  
  return(nmds_plot)
}

make_nmds_plot_plot_centroid <- function(nmds_object, nmds_matrix) {
  # calculate fit of environmental data (moisture and timepoint)
  fit <- envfit(nmds_object, nmds_matrix, permutations = 999)
  factors <- as.data.frame(vegan::scores(fit, display = "factors")) %>%
    filter(grepl("plot", rownames(.)))
  
  
  #add significant factor (moisture) to plot as psuedo vectors
  nmds_plot_vector <- nmds_plot + geom_segment(data = factors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                                               arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1) +
    geom_text(data = factors, aes(x = NMDS1, y = NMDS2, label = rownames(factors)),
              color = "red", vjust = -1, fontface = "italic") 
  
  #add significant factor (moisture) to plot as centroids
  nmds_plot_centroid <- nmds_plot + geom_point(data = factors, aes(x = NMDS1, y = NMDS2), color = "red") +
    geom_text(data = factors, aes(x = NMDS1, y = NMDS2, label = rownames(factors)),
              color = "red", vjust = -1, fontface = "italic") 
  print(fit)
  return(nmds_plot_centroid)
}

make_nmds_plot_time_no_t0 <- function(nmds_matrix) {
  nmds_plot <-  ggplot(mapping = aes(x, y)) +
    geom_point(data=nmds_matrix, aes(x=NMDS1, y=NMDS2, col=timepoint, shape=moisture), 
               size=size/2, show.legend = TRUE) +
    theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
    scale_color_viridis(discrete = TRUE, direction = -1, labels = c("3h", "24h", "48h", "72h", "168h")) +
    scale_shape_manual(values= list_of_shapes, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
    
    theme( text = element_text(size=size+5),
           legend.text = element_text(size=size+5, face="bold"),
           legend.title = element_blank(),
           #legend.key.size = unit(0.6, "cm"),
           #legend.key.width = unit(0.6,"cm"),
           legend.position = "bottom",
           axis.title.x = element_text(size=size+3,face="bold"),
           axis.title.y = element_text(size=size+3,face="bold"),
           plot.title = element_text(size=size+3,face="bold"),
           panel.grid.major = element_blank(),   # ← removes major grid lines
           panel.grid.minor = element_blank() )# ← removes minor  
  return(nmds_plot)
  
}

