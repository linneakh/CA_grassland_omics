### FUNCTIONS FOR PCA and nmds

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

make_screeplot <- function(eigen, dimensions){
  scree <- ggplot(data=eigen, aes(x=as.factor(dimensions), y=round(variance.percent,digits=2))) +
    geom_bar(stat="identity", fill='steelblue', color="black") +
    theme_linedraw(base_size = size) +
    xlab('Dimensions') +
    ylab('Explained variance (%)') + ylim(0,max(eigen$variance.percent)*1.3) +
    geom_text(data=eigen, aes(label=paste0(round(variance.percent,digits=2),'%')), size=4, vjust=-0.7, hjust=-0.1, angle=30) +
    geom_point(data=eigen, aes(x=as.factor(dimensions), y=variance.percent), color='black') +
    geom_path(data=eigen, aes(x=as.numeric(dimensions), y=variance.percent), color='black') +
    theme(axis.title.x = element_text(size=size+2,face="bold"),
          axis.title.y = element_text(size=size+2,face="bold"),
          axis.text = element_text(size=size-2),
          plot.title = element_text(size=size+2,face="bold"))
  return(scree)
}

make_cumvar <- function(eigen, dimensions){
  cumvar <- ggplot(data=eigen, aes(x=as.factor(dimensions), y=cumulative.variance.percent/100, group=1)) +
    geom_line(size=1.2, color='black') +
    geom_point(size=3, color='black') +
    xlab('PC axes') + ylab('Amount of explained variance') +
    theme_linedraw(base_size = size) + ylim(0, 1.05) +
    theme(axis.title.x = element_text(size=size+2,face="bold"),
          axis.title.y = element_text(size=size+2,face="bold"),
          axis.text = element_text(size=size-2),
          plot.title = element_text(size=size+2,face="bold"))
  return(cumvar)
}

make_pca_plot <- function(pca_coordinates){
  pca_plot <-  ggplot(mapping = aes(x, y)) +
    geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Treatment, shape=Depth), 
               size=size/3, show.legend = TRUE) +
    # stat_conf_ellipse(data=pca_coordinates, aes(x=PC1, y=PC2, color= XX), 
    #                   alpha=0, geom='polygon', show.legend = FALSE) +
    theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
    scale_color_manual(values= list_of_colors) +
    scale_shape_manual(values= list_of_shapes) +
    theme( legend.text = element_text(face="bold"),
           legend.title = element_blank(),
           legend.key.size = unit(0.6, "cm"),
           legend.key.width = unit(0.6,"cm"),
           axis.title.x = element_text(size=size+2,face="bold"),
           axis.title.y = element_text(size=size+2,face="bold"),
           plot.title = element_text(size=size+2,face="bold")) 
  return(pca_plot)
}

get_arrows <- function(pca, pca_coordinates){
  # prepare arrows
  # extract variable coordinates
  vars <- facto_summarize(pca, element = "var", result=c("coord","contrib","cos2"), axes=c(1,2))
  colnames(vars)[2:3] <-  c("xend", "yend")
  # rescale variable coordinates
  r <- min( (max(pca_coordinates[,"PC1"])-min(pca_coordinates[,"PC1"])/(max(vars[,"xend"])-min(vars[,"xend"]))),
            (max(pca_coordinates[,"PC2"])-min(pca_coordinates[,"PC2"])/(max(vars[,"yend"])-min(vars[,"yend"]))) )
  # multiply vars data by r before biplotting
  # in the original code (https://github.com/kassambara/factoextra/blob/master/R/fviz_pca.R), it calculates r for rescalling and
  # multiplies by 0.7 for plotting... IDK why...
  vars[,c("xend", "yend")] <- vars[,c("xend", "yend")]*r
  arrows <- as.data.frame(vars)
  return(arrows) # columns names, xend, yend for biplot :)
}

make_pca_biplot <- function(pca_coordinates, arrows){
  pca_biplot <- make_pca_plot(pca_coordinates) +
    new_scale_color() +
    geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend, color=contrib),
                 arrow=arrow(length = unit(0.1,"cm")), size=0.6, alpha = 0.5) +
    scale_color_gradient(low="#D9F0A3",high="#005A32") +
    geom_text_repel(data=arrows, aes(x=xend, y=yend, color=contrib),
                    label=arrows$name, size=size/4,show.legend = FALSE) +
    guides(color = guide_colourbar(barwidth = 1, barheight = 5))
  return(pca_biplot)
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

make_correlation_plots <- function(df) {
p <- ggplot(data = df, aes_string(x="sum.vst", y="CUE.metric.log", color = "Timepoint", shape = "Moisture")) +
  geom_point(size = 3) +
  scale_color_viridis(discrete = TRUE, 
                      direction = -1, 
                      labels = c("3h", "24h", "48h", "72h", "168h")) +
  scale_shape_manual(values = shapes_moisture, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
  ylab("CGE metric (log)") +
  xlab("transcript abundance (sum vst)") +
  geom_smooth(aes(group = "Moisture"), method = "lm", se = FALSE, color = "red4", size = 0.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

}

make_correlation_plots_with_legend <- function(df) {
  p <- ggplot(data = df, aes_string(x="sum.vst", y="CUE.metric.log", color = "Timepoint", shape = "Moisture")) +
    geom_point(size = 3) +
    scale_color_viridis(discrete = TRUE, 
                        direction = -1, 
                        labels = c( "3h", "24h", "48h", "72h", "168h")) +
    scale_shape_manual(values = shapes_moisture, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
    ylab("CGE metric (log)") +
    xlab("transcript abundance (sum vst)") +
    geom_smooth(aes(group = "Moisture"), method = "lm", se = FALSE, color = "red4", size = 0.5) +
    guides(color = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size=11, face="bold"),
      legend.title = element_blank(),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
}
