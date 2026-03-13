# pattern analysis plot

select_colomns_with_distinct_rows <- function(df, mode) {
  if (mode == "gpep") {
    col_ID <- c(
      "First_acc",
      "First_entry",
      "GN",
      "Description",
      "Glycosite_in_first_protein",
      "Nsite",
      "Nsite_seq",
      "Atypical",
      "Nsite_have",
      "GlycanComposition",
      "Hex",
      "HexNAc",
      "NeuAc",
      "NeuGc",
      "Fuc",
      "HighMan",
      "CoreFucosed",
      "Bisecting",
      "Category",
      "site_ID",
      "acc_site_glycan_ID",
      "acc_category_ID",
      "glycan_ID"
    )
    
    sample_col <- c("_imp$")
  } else if(mode == "igot") {
    col_ID <- c(
      "First_acc",
      "First_entry",
      "GN",
      "Description",
      "Glycosite_in_first_protein",
      "Nsite",
      "Nsite_seq",
      "Atypical",
      "Nsite_have",
      "acc_site_glycan_ID",
      "site_ID"
    )
    
    sample_col <- c("_igot$")
  } else if (mode == "ratio") {
    col_ID <- c(
      "First_acc",
      "First_entry",
      "GN",
      "Description",
      "Glycosite_in_first_protein",
      "Nsite",
      "Nsite_seq",
      "Atypical",
      "Nsite_have",
      "GlycanComposition",
      "Hex",
      "HexNAc",
      "NeuAc",
      "NeuGc",
      "Fuc",
      "HighMan",
      "CoreFucosed",
      "Bisecting",
      "Category",
      "site_ID",
      "acc_site_glycan_ID",
      "acc_category_ID",
      "glycan_ID"
    )
    
    sample_col <- c("ratio")
  }
  
  df <- select(df, all_of(col_ID), matches(sample_col))
  df %>% distinct()
}

transform_input_zscore <- function(df, fraction, mode) {

  time_order <- c( "p0","p8", "2w", "4w", "3M")
  replicates <- c("_1", "_2", "_3")
   
  if(mode == "gpep"){
    ordered_cols <- lapply(time_order, function(tp) {
      paste0("mean_", tp, replicates, "_imp")}) %>% unlist
    rowname_column <- "acc_site_glycan_ID"
  } else if(mode == "igot"){
    ordered_cols <- lapply(time_order, function(tp) {
      paste0(fraction, "_", tp, replicates, "_log2m_imp_igot")}) %>% unlist
    rowname_column <- "site_ID"
  } else if(mode == "ratio"){
    ordered_cols <- lapply(time_order, function(tp) {
      paste0("ratio_", tp, replicates)}) %>% unlist
    rowname_column <- "acc_site_glycan_ID"
  } else if(mode == "igot_cor"){
    ordered_cols <- lapply(time_order, function(tp) {
      paste0(fraction, "_", tp, replicates, "_log2m_imp_igot")}) %>% unlist
    rowname_column <- "acc_site_glycan_ID"
  }
  
  # select columns and set rownames
  rowname <- df[[rowname_column]]
  df <- df[, ordered_cols] 
  df <- as.matrix(df)
  rownames(df) <- rowname
  if(mode == "igot") {
    unique_rownames <- unique(rowname)
    df <- df[unique_rownames, ]
  } 
  
  # obtain z-score for each row
  df_zscore_colnames <- colnames(df)
  df_zscore <- t(apply(df, 1, scale))
  colnames(df_zscore) <- df_zscore_colnames
  
  df_zscore
}

extract_sig_vars <- function(df, fraction_mode) {
  
  
  design_sheet <- data.frame(
    Time       = c(0, 0, 0, 8, 8, 8, 14, 14, 14, 28, 28, 28, 84, 84, 84),
    Replicates = c(1, 1, 1, 2, 2, 2,  3,  3,  3,  4,  4,  4,  5,  5,  5),
    value      = c(1, 1, 1, 1, 1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
  )
  colnames(design_sheet)[3] <- fraction_mode
  rownames(design_sheet) <- colnames(df)
  
  d <- make.design.matrix(design_sheet, degree = 2)
  p <- p.vector(df, d, Q = 0.05)
  t <- T.fit(p)
  res <- get.siggenes(t, vars = "all") 
  
  res
}

plot_pattern_cluster <- function(pattern_res, k_value, fraction, mode, graph_mode, color_pattern) {
    time_order <- c(0, 8, 14, 28, 84)
    replicates <- 3
    
    sig_genes_data <- pattern_res$sig.genes  # extract significant genes
    sig_genes_df <- as.data.frame(sig_genes_data$sig.profiles)
    cluster_res <- see.genes(pattern_res$sig.genes, k = k_value)
    
    # rename significant genes colnames
    colnames(sig_genes_df) <- c("0_1", "0_2", "0_3",
                                "8_1", "8_2", "8_3",
                                "14_1", "14_2", "14_3",
                                "28_1", "28_2", "28_3",
                                "84_1", "84_2", "84_3")
    
    sig_genes_df$peptide <- rownames(sig_genes_df)
    
    # transform to long format
    long_data <- sig_genes_df %>%
      pivot_longer(
        cols = -peptide,
        names_to = c("Time", "Replicate"),
        names_sep = "_",
        values_to = "score"
      )
    long_data$Time <- as.numeric(long_data$Time)
    
    # add clustering result from cluster_res
    long_data$Cluster <-
      factor(cluster_res$cut[match(long_data$peptide, names(cluster_res$cut))])
    long_data_individual <- long_data
    
    # for each time point、grouping by cluster and replicate and taking the mean
    mean_data_per_time_individual <- long_data_individual %>%
      group_by(Cluster, Time, Replicate) %>%
      summarise(Individual_Expression = mean(score, na.rm = TRUE), .groups = 'drop')
    
    # for each time point, grouping by cluster and taking the mean
    mean_data_per_time <- long_data_individual %>%
      group_by(Cluster, Time) %>%
      summarise(Average_Expression = mean(score, na.rm = TRUE),
                .groups = 'drop')
    
    # count each cluster size
    cluster_counts <- long_data %>%
      group_by(Cluster) %>%
      summarise(n = n_distinct(peptide))
    cluster_labels <- paste0(cluster_counts$Cluster, " (n=", cluster_counts$n, ")")
    names(cluster_labels) <- cluster_counts$Cluster
    
    # Title
    cluster_pattern_title <- paste(fraction, mode, "ratio cluster pattern")
    all_pattern_title <- paste("All expression patterns of significant", 
                               fraction,
                               mode,"ratio")
    cluter_all_pattern_title <- paste("Expression patterns by cluster in", 
                                      fraction, mode, "ratio")
    
    # Plot the mean and individual abundance
    if (graph_mode == "cluster") {
      if (color_pattern == "incre") {
        color <- c("#ff3366" , "#4662D7FF",   "#ffa500", "#3cb371")
      } else {
        color <- c("#4662D7FF", "#ff3366" ,  "#ffa500", "#3cb371")
      }
      names(color) <- cluster_counts$Cluster
      
      p <-
        ggplot(mean_data_per_time_individual,
               aes(x = Time, y = Individual_Expression, color = Cluster, group = Cluster)) +
        geom_line(data = mean_data_per_time, 
                  aes(x = Time, y = Average_Expression, color = Cluster), 
                  linewidth = 1.5) +
        geom_vline(xintercept = c(0, 8, 14, 28, 84),color = "dimgray", linewidth = 0.5, linetype = "dotted") +
        geom_point(size = 3, shape = 16, alpha = 0.6) +
        labs(title = cluster_pattern_title, x = "",  y = "scale abundance") +
        scale_x_continuous(breaks = c(0, 8, 14, 28, 84),
                           labels = c("P0", "P8", "2w", "1M", "3M")) +
        scale_color_manual(values = color, labels = cluster_labels) +
        theme_minimal() +
        # theme(axis.text = element_text(size = 14, face = "bold", color = "black"),
        #       axis.title = element_text(size = 16, face = "bold", color = "black"),
        #       axis.line = element_line(color = "black", size = 1.2))
        # 
      theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background  = element_rect(fill = "white", color = NA),  
      panel.grid = element_blank(), 
      axis.text.y  = element_text(color = "black", size = 20, face = "bold"),
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),
      axis.text.x  = element_text(color = "black",size = 20, face = "bold"),                                          axis.line.x = element_line(color = "black", size = 1),
      axis.line.y = element_line(color = "black", size = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.length.y = unit(-3, "pt"),
      legend.text  = element_text(size = 15))
    }
    
    if (graph_mode == "all") {
      #全遺伝の発現変化のプロット
      long_data_avg <- long_data %>%
        group_by(peptide, Time, Cluster) %>%
        summarise(Expression = mean(score, na.rm = TRUE),.groups = 'drop')
      
      # ggplot2ですべての発現パターンをプロット
      p <-
        ggplot(long_data_avg,
               aes(
                 x = Time,
                 y = Expression,
                 group = peptide,
                 color = peptide
               )) +
        geom_line() +
        theme_minimal() +
        labs(title = all_pattern_title,
             x = "Time",
             y = "Scale abundance") +
        scale_x_continuous(
          breaks = c(0, 8, 14, 28, 84),
          labels = c("P0", "P8", "2w", "1M", "3M")
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position = "none")
      
    }
    
    if (graph_mode == "cluster_all") {
      #全遺伝の発現変化のプロット
      long_data_avg <- long_data %>%
        group_by(peptide, Time, Cluster) %>%
        summarise(Expression = mean(score, na.rm = TRUE),
                  .groups = 'drop')
      
      # cluster別に発現パターンをプロット
      p <-
        ggplot(long_data_avg,
               aes(
                 x = Time,
                 y = Expression,
                 group = peptide,
                 color = peptide
               )) +
        geom_line() +
        facet_wrap( ~ Cluster, scales = "free_y") +
        theme_minimal() +
        labs(title = cluter_all_pattern_title,
             x = "Time",
             y = "Scale abundance") +
        scale_x_continuous(
          breaks = c(0, 8, 14, 28, 84),
          labels = c("P0", "P8", "2w", "1M", "3M")
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    }
    
    p
}

