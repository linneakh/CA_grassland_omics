# User-supplied plotting functions (kept verbatim) -------------------------
make_correlation_plots <- function(df, x_label="sum VST") {
  p <- ggplot(data = df, aes_string(x = "sum.vst", y = "CUE.metric.log", color = "Timepoint", shape = "Moisture")) +
    geom_point(size = 3) +
    scale_color_viridis(discrete = TRUE,
                        direction = -1,
                        labels = c("3h", "24h", "48h", "72h", "168h")) +
    scale_shape_manual(values = shapes_moisture, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
    ylab("CGE metric (log)") +
    xlab(x_label) +
    geom_smooth(aes(group = "Moisture"), method = "lm", se = FALSE, color = "red4", size = 0.5) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}


make_correlation_plots_with_legend <- function(df) {
  p <- ggplot(data = df, aes_string(x = "sum.vst", y = "CUE.metric.log", color = "Timepoint", shape = "Moisture")) +
    geom_point(size = 3) +
    scale_color_viridis(discrete = TRUE,
                        direction = -1,
                        labels = c("3h", "24h", "48h", "72h", "168h")) +
    scale_shape_manual(values = shapes_moisture, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
    ylab("CGE metric (log)") +
    xlab("transcript abundance (sum vst)") +
    geom_smooth(aes(group = "Moisture"), method = "lm", se = FALSE, color = "red4", size = 0.5) +
    guides(color = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11, face = "bold"),
      legend.title = element_blank(),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.title.justification = "center"
    )
}

# Helpers to add stats and titles ------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  sprintf("p = %.1g", p)  # always 1 significant figure
}



sig_stars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

add_r_p_annot <- function(p, dat, x = "sum.vst", y = "CUE.metric.log", fdr_correct = FALSE) {
  # Prepare data and guard against missing/constant cases
  xv <- dat[[x]]; yv <- dat[[y]]
  ok <- is.finite(xv) & is.finite(yv)
  n_ok <- sum(ok)
  
  safe_range <- function(v) {
    r <- range(v[is.finite(v)], na.rm = TRUE)
    if (any(!is.finite(r))) c(0, 1) else if (r[1] == r[2]) r + c(-0.5, 0.5) else r
  }
  
  xrng <- safe_range(xv)
  yrng <- safe_range(yv)
  x_pos <- xrng[1] + 0.98 * diff(xrng)   # 95% to the right
  y_pos <- yrng[1] + 0.025 * diff(yrng)   # 15% up from bottom for visibility
  
  # Default label if not enough points for cor.test
  if (n_ok < 3) {
    return(
      p + annotate(
        geom = "text", x = x_pos, y = y_pos, hjust = 1, vjust = 0,
        label = "r = NA\np = NA (n<3)", size = 3, color = "black",
        fontface = "plain"
      )
    )
  }
  
  # Pearson r and p
  ct <- suppressWarnings(cor.test(xv[ok], yv[ok], method = "pearson"))
  r  <- unname(ct$estimate)
  pv <- unname(ct$p.value)
  if (fdr_correct) pv <- p.adjust(pv, method = "fdr")
  
  stars <- sig_stars(pv)
  label_text <- sprintf("r = %.2f\n%s %s", r, fmt_p(pv), stars)
  
  p + annotate(
    geom = "text", x = x_pos, y = y_pos, hjust = 1, vjust = 0,
    label = label_text, size = 3, color = "black",
    fontface = ifelse(stars != "", "bold", "plain")
  )
}



# Data wrangling helpers ----------------------------------------------------
# Intersect curated KOs with VIP list, sum VST per-sample for those KOs
make_sum_table <- function(cue_timepoint, vip_kos, ko_subset) {
  keep <- intersect(vip_kos, ko_subset)
  if (length(keep) == 0) {
    return(tibble(SampleID = character(), CUE.metric.log = numeric(), sum.vst = numeric(), Timepoint = factor(), Moisture = factor()))
  }
  cue_timepoint %>%
    dplyr::select(SampleID, CUE.metric.log, all_of(keep), Timepoint, Moisture) %>%
    pivot_longer(cols = all_of(keep), names_to = "KO", values_to = "VST") %>%
    group_by(SampleID, Timepoint, Moisture) %>%
    summarise(sum.vst = sum(VST, na.rm = TRUE), CUE.metric.log = mean(CUE.metric.log, na.rm = TRUE), .groups = "drop") %>%
    mutate(Timepoint = factor(Timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168")),
           Moisture = factor(Moisture, levels = c("x50", "x100")))
}

save_plot_duo <- function(plot, file_stub, width = 3.25, height = 3, dpi = 300) {
  ggsave(paste0(file_stub, ".png"), plot = plot, width = width, height = height, units = "in", dpi = dpi)
  ggsave(paste0(file_stub, ".pdf"), plot = plot, width = width, height = height, units = "in", dpi = dpi)
}

# Data wrangling helpers ----------------------------------------------------
# Intersect curated KOs with VIP list, sum VST per-sample for those KOs
make_sum_table_no_intersect <- function(cue_timepoint, ko_subset) {
  if (length(keep) == 0) {
    return(tibble(SampleID = character(), CUE.metric.log = numeric(), sum.vst = numeric(), Timepoint = factor(), Moisture = factor()))
  }
  cue_timepoint %>%
    dplyr::select(SampleID, CUE.metric.log, ko_subset, Timepoint, Moisture) %>%
    pivot_longer(cols = ko_subset, names_to = "KO", values_to = "VST") %>%
    group_by(SampleID, Timepoint, Moisture) %>%
    summarise(sum.vst = sum(VST, na.rm = TRUE), CUE.metric.log = mean(CUE.metric.log, na.rm = TRUE), .groups = "drop") %>%
    mutate(Timepoint = factor(Timepoint, levels = c("t0", "t3", "t24", "t48", "t72", "t168")),
           Moisture = factor(Moisture, levels = c("x50", "x100")))
}

save_plot_duo <- function(plot, file_stub, width = 3.25, height = 3, dpi = 300) {
  ggsave(paste0(file_stub, ".png"), plot = plot, width = width, height = height, units = "in", dpi = dpi)
  ggsave(paste0(file_stub, ".pdf"), plot = plot, width = width, height = height, units = "in", dpi = dpi)
}


# New: correlation plot with two regression lines (one per Moisture) and separate r/p labels
make_correlation_plots_by_moisture <- function(df, show_legend = FALSE) {
  base <- ggplot(df, aes(x = sum.vst, y = CUE.metric.log, color = Timepoint, shape = Moisture)) +
    geom_point(size = 3) +
    scale_color_viridis(discrete = TRUE, direction = -1,
                        labels = c("3h", "24h", "48h", "72h", "168h")) +
    scale_shape_manual(values = shapes_moisture, breaks = c("x50", "x100"), labels = c("50%", "100%")) +
    ylab("CGE metric (log)") +
    xlab("transcript abundance (sum vst)") +
    theme_bw() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  base +
    geom_smooth(
      inherit.aes = FALSE,
      data = df,
      aes(x = sum.vst, y = CUE.metric.log, linetype = Moisture, group = Moisture),
      method = "lm", se = FALSE, color = "red4", linewidth = 0.6
    ) +
    scale_linetype_manual(values = c(x50 = "solid", x100 = "dashed"))
}



add_r_p_annot_by_moisture <- function(p, dat, group_var = "Moisture", x = "sum.vst", y = "CUE.metric.log", fdr_correct = FALSE) {
  # Compute group-wise stats (one row per moisture)
  stats <- dat %>%
    dplyr::filter(!is.na(.data[[group_var]])) %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      .groups = "drop",
      n = sum(is.finite(.data[[x]]) & is.finite(.data[[y]])),
      r = {
        xv <- .data[[x]]; yv <- .data[[y]]; ok <- is.finite(xv) & is.finite(yv)
        if (sum(ok) < 3) NA_real_ else suppressWarnings(unname(cor(xv[ok], yv[ok], method = "pearson")))
      },
      p = {
        xv <- .data[[x]]; yv <- .data[[y]]; ok <- is.finite(xv) & is.finite(yv)
        if (sum(ok) < 3) NA_real_ else suppressWarnings(unname(cor.test(xv[ok], yv[ok], method = "pearson")$p.value))
      }
    )
  
  # FDR (optional)
  if (fdr_correct) stats$p <- p.adjust(stats$p, method = "fdr")
  
  # Vector-safe formatting of p and stars (sig_stars/fmt_p are scalar helpers)
  stats <- stats %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      grp = as.character(.data[[group_var]]),
      stars = sig_stars(p),          # safe now because rowwise -> scalar
      p_fmt = fmt_p(p),              # formatted p-value per row
      label = ifelse(
        is.na(p),
        paste0(grp, ": r = NA
p = NA (n<3)"),
        sprintf("%s: r = %.2f
%s %s", grp, r, p_fmt, stars)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(grp)
  
  # Build positions from panel ranges (use true rendered panel limits)
  b <- ggplot_build(p %+% dat)
  xr <- b$layout$panel_scales_x[[1]]$range$range
  yr <- b$layout$panel_scales_y[[1]]$range$range
  if (!is.finite(diff(xr)) || diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  if (!is.finite(diff(yr)) || diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
  x_left  <- xr[1] + 0.02 * diff(xr)
  x_right <- xr[2] - 0.02 * diff(xr)
  y_base  <- yr[1] + 0.02 * diff(yr)
  
  # Assign corners: first group left/bottom, second right/bottom
  stats <- stats %>%
    dplyr::mutate(
      row_id = dplyr::row_number(),
      x = dplyr::case_when(row_id == 1 ~ x_left,
                           row_id == 2 ~ x_right,
                           TRUE ~ x_left),
      hjust = dplyr::case_when(row_id == 1 ~ 0,
                               row_id == 2 ~ 1,
                               TRUE ~ 0),
      y = y_base,
      is_sig = !is.na(p) & p < 0.05
    )
  
  p +
    geom_text(data = dplyr::filter(stats, !is_sig),
              aes(x = x, y = y, label = label, hjust = hjust),
              vjust = 0, size = 3, color = "black", inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(stats, is_sig),
              aes(x = x, y = y, label = label, hjust = hjust),
              vjust = 0, size = 3, color = "black", fontface = "bold", inherit.aes = FALSE)
}

# get r and p from a data.frame used in a panel
compute_r_p <- function(dat, x = "sum.vst", y = "CUE.metric.log") {
  ok <- is.finite(dat[[x]]) & is.finite(dat[[y]])
  if (sum(ok) < 3) return(list(r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(dat[[x]][ok], dat[[y]][ok], method = "pearson"))
  list(r = unname(ct$estimate), p = unname(ct$p.value))
}

# same as your add_r_p_annot but lets you inject an adjusted p (and/or r)
add_r_p_annot_fixed <- function(p, dat, x = "sum.vst", y = "CUE.metric.log",
                                pv_override = NULL, r_override = NULL) {
  xv <- dat[[x]]; yv <- dat[[y]]
  ok <- is.finite(xv) & is.finite(yv)
  n_ok <- sum(ok)
  
  safe_range <- function(v) {
    r <- range(v[is.finite(v)], na.rm = TRUE)
    if (any(!is.finite(r))) c(0,1) else if (r[1]==r[2]) r + c(-0.5,0.5) else r
  }
  xrng <- safe_range(xv); yrng <- safe_range(yv)
  x_pos <- xrng[1] + 0.95*diff(xrng)
  y_pos <- yrng[1] + 0.15*diff(yrng)
  
  if (n_ok < 3) {
    return(p + annotate("text", x = x_pos, y = y_pos, hjust = 1, vjust = 0,
                        label = "r = NA\np = NA (n<3)", size = 3, color = "black"))
  }
  
  ct <- suppressWarnings(cor.test(xv[ok], yv[ok], method = "pearson"))
  r_use <- if (is.null(r_override)) unname(ct$estimate) else r_override
  p_use <- if (is.null(pv_override)) unname(ct$p.value) else pv_override
  
  stars <- sig_stars(p_use)
  lab <- sprintf("r = %.2f\n%s %s", r_use, fmt_p(p_use), stars)
  
  p + annotate("text", x = x_pos, y = y_pos, hjust = 1, vjust = 0,
               label = lab, size = 3, color = "black",
               fontface = ifelse(stars != "", "bold", "plain"))
}

# Collect p-values from multiple datasets (across all panels), adjust with FDR, and pass them to annotation
annotate_panels_with_fdr <- function(p_list, data_list, x = "sum.vst", y = "CUE.metric.log", bottom_row_indices = NULL) {
  # Extract raw p-values from each dataset (per panel)
  get_p <- function(dat) {
    xv <- dat[[x]]; yv <- dat[[y]]
    ok <- is.finite(xv) & is.finite(yv)
    if (sum(ok) < 3) return(NA_real_)
    suppressWarnings(cor.test(xv[ok], yv[ok], method = "pearson")$p.value)
  }
  raw_p <- vapply(data_list, get_p, numeric(1))
  
  # Adjust all p-values together across ALL panels
  adj_p <- p.adjust(raw_p, method = "fdr")
  
  # Annotate each panel with its adjusted p-value and hide x-axis labels except bottom row
  annotated_plots <- mapply(function(p, dat, pval, idx) {
    xv <- dat[[x]]; yv <- dat[[y]]
    ok <- is.finite(xv) & is.finite(yv)
    n_ok <- sum(ok)
    safe_range <- function(v) {
      r <- range(v[is.finite(v)], na.rm = TRUE)
      if (any(!is.finite(r))) c(0, 1) else if (r[1] == r[2]) r + c(-0.5, 0.5) else r
    }
    xrng <- safe_range(xv)
    yrng <- safe_range(yv)
    x_pos <- xrng[1] + 0.95 * diff(xrng)
    y_pos <- yrng[1] + 0.15 * diff(yrng)
    if (n_ok < 3) {
      p <- p + annotate("text", x = x_pos, y = y_pos, hjust = 1, vjust = 0,
                        label = "r = NA\np = NA (n<3)", size = 3, color = "black")
    } else {
      r_val <- unname(cor(xv[ok], yv[ok], method = "pearson"))
      stars <- sig_stars(pval)
      label_text <- sprintf("r = %.2f\n%s %s", r_val, fmt_p(pval), stars)
      p <- p + annotate("text", x = x_pos, y = y_pos, hjust = 1, vjust = 0,
                        label = label_text, size = 3, color = "black",
                        fontface = ifelse(stars != "", "bold", "plain"))
    }
    if (!is.null(bottom_row_indices) && !(idx %in% bottom_row_indices)) {
      p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    p
  }, p_list, data_list, adj_p, seq_along(p_list), SIMPLIFY = FALSE)
  
  annotated_plots
}




