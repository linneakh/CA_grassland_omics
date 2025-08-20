

#create function to recreate Steve's bar plots of birth rates
bar_plot_b_fixed <- function(rate_file){
  rate_file <- as.data.frame(rate_file)
  
  df <- rate_file |>
    drop_na() |>
    filter(type == "b.boot.median") |>
    separate(trt.code.2, c("Timepoint", "Isotope"), remove = FALSE) |>
    mutate(Timepoint = factor(Timepoint, levels = c("t3", "t24", "t48", "t72", "t168"))) |>
    tidyr::separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") |>
    mutate(Phylum2 = case_when(
      Phylum == "p_Armatimonadetes" ~ "Phylum < 2%",
      Phylum == "p_candidate division WPS-1" ~ "Phylum < 2%",
      Phylum == "p_candidate division WPS-2" ~ "Phylum < 2%",
      Phylum == "p_Deinococcus-Thermus" ~ "Phylum < 2%",
      Phylum == "p_Thaumarchaeota" ~ "Phylum < 2%",
      Phylum == "p_Fusobacteria" ~ "Phylum < 2%",
      Phylum == "p_Chloroflexi" ~ "Phylum < 2%",
      Phylum == "p_Nitrospirae" ~ "Phylum < 2%",
      TRUE ~ Phylum
    )) |>
    mutate(Phylum2 = str_remove(Phylum2, "p_"))
  
  rate_plot <- ggplot(df, aes(x=Timepoint, y=median, fill = Phylum2)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~type + moisture,dir = "h",  ncol=2) +
    scale_fill_manual(values = col_list) + 
    scale_x_discrete(labels = c("3", "24", "48", "72", "168")) +
    labs(fill = "Phylum") +
    theme_bw() +
    theme(text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(rate_plot)
}

#create function to recreate Steve's bar plots of birth rates
bar_plot_d_fixed <- function(rate_file){
  rate_file <- as.data.frame(rate_file)
  
  df <- rate_file |>
    drop_na() |>
    filter(type == "d.boot.median") |>
    separate(trt.code.2, c("Timepoint", "Isotope"), remove = FALSE) |>
    mutate(Timepoint = factor(Timepoint, levels = c("t3", "t24", "t48", "t72", "t168"))) |>
    tidyr::separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") |>
    mutate(Phylum2 = case_when(
      Phylum == "p_Armatimonadetes" ~ "Phylum < 2%",
      Phylum == "p_candidate division WPS-1" ~ "Phylum < 2%",
      Phylum == "p_candidate division WPS-2" ~ "Phylum < 2%",
      Phylum == "p_Deinococcus-Thermus" ~ "Phylum < 2%",
      Phylum == "p_Thaumarchaeota" ~ "Phylum < 2%",
      Phylum == "p_Fusobacteria" ~ "Phylum < 2%",
      Phylum == "p_Chloroflexi" ~ "Phylum < 2%",
      Phylum == "p_Nitrospirae" ~ "Phylum < 2%",
      TRUE ~ Phylum
    )) |>
    mutate(Phylum2 = str_remove(Phylum2, "p_"))
  
  rate_plot <- ggplot(df, aes(x=Timepoint, y=median, fill = Phylum2)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~type + moisture,dir = "h",  ncol=2) +
    scale_fill_manual(values = col_list) + 
    scale_x_discrete(labels = c("3", "24", "48", "72", "168")) +
    labs(fill = "Phylum") +
    theme_bw() +
    theme(text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(rate_plot)
}

#create function to recreate Steve's bar plots of birth rates
bar_plot_d <- function(rate_file){
  rate_file <- as.data.frame(rate_file)
  
  df <- rate_file |>
    drop_na() |>
    filter(type == "d.boot.median") |>
    separate(trt.code.2, c("Timepoint", "Isotope"), remove = FALSE) |>
    mutate(Timepoint = factor(Timepoint, levels = c("t3", "t24", "t48", "t72", "t168"))) |>
    tidyr::separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") |>
    mutate(Phylum2 = case_when(
      Phylum == "p_Armatimonadetes" ~ "Phylum < 2%",
      Phylum == "p_candidate division WPS-1" ~ "Phylum < 2%",
      Phylum == "p_candidate division WPS-2" ~ "Phylum < 2%",
      Phylum == "p_Deinococcus-Thermus" ~ "Phylum < 2%",
      Phylum == "p_Thaumarchaeota" ~ "Phylum < 2%",
      Phylum == "p_Fusobacteria" ~ "Phylum < 2%",
      Phylum == "p_Chloroflexi" ~ "Phylum < 2%",
      Phylum == "p_Nitrospirae" ~ "Phylum < 2%",
      TRUE ~ Phylum
    )) |>
    mutate(Phylum2 = str_remove(Phylum2, "p_"))
  
  rate_plot <- ggplot(df, aes(x=Timepoint, y=median, fill = Phylum2)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~type + moisture, dir = "h", scales = "free", ncol=2) +
    scale_fill_manual(values = col_list) + 
    scale_x_discrete(labels = c("3", "24", "48", "72", "168")) +
    labs(fill = "Phylum") +
    theme_bw() +
    theme(text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(rate_plot)
}


  