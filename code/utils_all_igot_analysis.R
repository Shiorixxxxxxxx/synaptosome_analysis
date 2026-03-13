count_of_glycan<-function(df){
  
  count_df<-df%>%
    group_by(site_ID)%>%
    mutate(nGlycan = n())
  
  count_df
}

count_of_site<-function(df){
  
  count_df<-df%>%
    group_by(acc)%>%
    mutate(nsite = n())
  
  count_df
}

insite_glycan_violine_plot<-function(plot_df, x_col, values_color, col){
  
  ggplot(plot_df, aes(x=!!sym(x_col), y=!!sym(col), fill = !!sym(x_col))) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.8, width = 0.6) + 
    #geom_point(color = "#dcdcdc", size = 0.8, width = 0.15, alpha = 0.8) +
    scale_fill_manual(values = values_color) +
    stat_summary(fun = mean, geom = "point", color = "#515A5A", size = 2, alpha =0.8) +
    #stat_compare_means(method = "t.test", label = "p.signif") +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background  = element_rect(fill = "white", color = NA),  
      panel.grid = element_blank(), 
      axis.text.y  = element_text(color = "black", size = 20, face = "bold"),
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),
      axis.text.x  = element_text(color = "black",size = 20, face = "bold"),                               
      axis.ticks.x = element_blank())                             
}

select_col<-function(df){
  
  target_col<-c("First_acc",
                "First_entry",
                "GN",
                "Description",
                "Glycosite_in_first_protein",
                "Nsite",                     
                "Nsite_seq",
                "Atypical", 
                "Nsite_have",                
                "site_ID",
                "igot_cluster")
  
  ref_df<-select(df, all_of(target_col))
  
}

col_name_change<-function(df, fraction){
  
  
  old_name<-colnames(df)
  
  if(fraction == "sps"){
    new_name<-old_name%>%
      sub(".*sps_([^_]+_[^_]+)_.*", "\\1", .) |>
      gsub("p0", "p0", x = _) |>
      gsub("p8", "p8", x = _) |>
      gsub("2w", "p14", x = _) |>
      gsub("4w", "p28", x = _) |>
      gsub("3M", "p84", x = _)}else{
        new_name<-old_name%>%
          sub(".*whole_([^_]+_[^_]+)_.*", "\\1", .) |>
          gsub("p0", "p0", x = _) |>
          gsub("p8", "p8", x = _) |>
          gsub("2w", "p14", x = _) |>
          gsub("4w", "p28", x = _) |>
          gsub("3M", "p84", x = _)
      }
  
  colnames(df)<-new_name
  
  df
  
}

venn_analysis<-function(list,color){
  venn<-venndir(list,
                proportional=TRUE,
                set_colors = color,
                show_labels="cs",
                show_segments=FALSE,
                inside_percent_threshold=0,
                overlap_type="overlap",
                font_cex = c(1.5, 1.5, 0.8))
  venn}

select_igot_position<-function(df){
  
  extract_col<-colnames(df)
  df_expanded <- df %>%
    mutate(site_split = str_split(Nigot_seq, "/")) %>%   # "/" で分割
    unnest(site_split) %>%
    group_by(`Master.Protein.Accessions`)%>%
    mutate(pos = str_extract(site_split, "^[0-9]+"),
           is_single = !str_detect(Nigot_seq, "/"))%>%ungroup()
  
  df_flag <- df_expanded %>%
    group_by(`Master.Protein.Accessions`, pos)%>%
    mutate(single_exists = any(is_single),
           Atypical = str_detect(site_split, "N.S|N.T|N.C")) %>%
    ungroup()
  
  df_filtered <- df_flag %>%
    group_by(`Master.Protein.Accessions`)%>%
    filter(!(single_exists & !is_single)) %>% 
    filter(!(!Atypical)) %>% 
    distinct(Nigot_seq, .keep_all = TRUE) %>% ungroup()%>%
    select(all_of(extract_col))
  
  df_filtered
}

heat_mat<-function(df){
  
  acc_GN_list<-select(df, acc, GN, acc_site_ID)
  df_sorted <- df[order(df$acc), ]
  mat <- df[, grep("^p", colnames(df_sorted))]
  rownames(mat)<-df$acc_site_ID
  
  mat
  
}

transform_input_zscore <- function(df) {
  
  rownames(df)<-df$acc_site_ID
  df<-select(df, -acc_site_ID)
  
  time_order <- c( "p0","p8", "p14", "p28", "p84")
  replicates <- c("_1", "_2", "_3")
  
  ordered_cols <- lapply(time_order, function(tp) {
    paste0(tp, replicates)}) %>% unlist
  
  # select columns and set rownames
  rowname <- rownames(df)
  df <- df[, ordered_cols] 
  df <- as.matrix(df)
  rownames(df) <- rowname
  unique_rownames <- unique(rowname)
  df <- df[unique_rownames, ]
  
  # obtain z-score for each row
  df_zscore_colnames <- colnames(df)
  df_zscore <- t(apply(df, 1, scale))
  colnames(df_zscore) <- df_zscore_colnames
  
  df_zscore_clean <- df_zscore %>%
    as.data.frame() %>%
    filter(!if_all(everything(), ~ is.na(.) | is.nan(.)))
  
  df_zscore_clean
}

add_pro_info<-function(mat, origine_df){
  
  mat_df <- data.frame(acc_site_ID = rownames(mat), mat)
  
  acc_GN_list<-select(origine_df, "acc", "GN", "acc_site_ID")
  mat_joined <- mat_df %>%
    left_join(acc_GN_list, by = "acc_site_ID")
  mat_joined <- mat_joined %>%
    relocate(acc, GN, .after = acc_site_ID)
  
  rownames(mat_joined)<-mat_joined$acc_site_ID
  
  mat_joined<-select(mat_joined, -acc_site_ID)
  
  mat_joined
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

silhouette_analysis<-function(res){
  
  k_range <- 2:10
  res<-sps_igot_pattern
  
  silhouette_scores <- sapply(k_range, function(k) {
    pam_result <- pam(res$sig.genes$sig.profiles, k = k)
    pam_result$silinfo$avg.width
  })
  
  plot(k_range, silhouette_scores, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of Clusters", ylab = "Average Silhouette Width",
       main = "Silhouette Method for Optimal Clusters")
  
  optimal_k <- k_range[which.max(silhouette_scores)]
  
  optimal_k}

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
  cluster_pattern_title <- paste(fraction, mode, "cluster pattern")
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
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background  = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(), 
        axis.text.y  = element_text(color = "black", size = 15, face = "bold"),
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        axis.text.x  = element_text(color = "black",size = 15, face = "bold"),                                          
        axis.line.x = element_line(color = "black", size = 1),
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

obtain_cluster_with_ID <- function(cluster, cluster_name) {
  cluster_num_df <- cluster$cut %>% data.frame() %>% setNames(cluster_name)
  cluster_num_df$site_ID <- rownames(cluster_num_df)
  rownames(cluster_num_df) <- NULL
  cluster_num_df
}

extract_items <- function(label_df, target_label) {
  key<-label_df %>%
    dplyr::filter(venn_label == target_label) %>%
    dplyr::pull(items) %>%
    unlist()
  
  if(target_label == "overlap"){
    total_both_acc <- label_df %>%
      dplyr::filter(is.na(venn_label)) %>%
      dplyr::pull(items) %>%
      unlist()
    
    key<-c(total_both_acc, key)
  }
  key
}

GO_analysis <- function(sel_df, GO_mode) {
  
  # 1. GNリスト（例: Uniprot ID → Entrez ID マッピング）
  GN_list <- unique(sel_df$GN)
  entrez_ids <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = GN_list,
    columns = c("ENTREZID"),
    keytype = "SYMBOL"
  )
  
  # 2. enrichGO 実行
  ego <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = GO_mode,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  ego_df <- as.data.frame(ego)
  
  ego_df }

GO_selfbg_analysis <- function(sel_df,bg_df, GO_mode) {
  
  GN_list <- unique(sel_df$GN)
  entrez_ids <- AnnotationDbi::select( org.Mm.eg.db, keys = GN_list, columns = c("ENTREZID"), keytype = "SYMBOL")
  
  bg_GN <- unique(bg_df$GN)
  bg_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = bg_GN, columns = c("ENTREZID"), keytype = "SYMBOL")
  
  ego <- enrichGO(
    gene = entrez_ids$ENTREZID,
    universe = bg_entrez$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = GO_mode,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  ego_df <- as.data.frame(ego)
  
  ego_df }

GO_plot<-function(ego_df, GO_mode){
  
  plot_df <- ego_df %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(
      Description = factor(Description, levels = rev(Description)),
      log10_p = -log10(p.adjust))
  
  p<-ggplot(plot_df, aes(x = log10_p, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(name = "FDR", low = "#1e90ff", high = "#ff6347", trans = "reverse") +
    labs(
      title = paste0("GO enrichment (", GO_mode, ")"),
      x = "-log10(FDR)",
      y = "GO term") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text.y  = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x  = element_text(color = "black",size = 10, face = "bold"),                                          
      axis.ticks.y = element_line(color = "black", size = 1))
      
  p}

GO_res_transGN<-function(GO_res, tar_ID){
  
  GN<-GO_res%>%
    filter(ID %in% tar_ID) %>%
    separate_rows(geneID, sep = "/") %>%
    transmute(ENTREZID = geneID) %>%
    distinct() %>%
    left_join(bitr(.$ENTREZID,
                   fromType = "ENTREZID",
                   toType   = "SYMBOL",
                   OrgDb    = org.Mm.eg.db),
              by = "ENTREZID") %>%
    pull(SYMBOL) %>%
    unique()
  
  GN} 

trans_longdata<-function(plot_df){
  
  
  time_point <- c("0_1", "0_2", "0_3",
                  "8_1", "8_2", "8_3",
                  "14_1", "14_2", "14_3",
                  "28_1", "28_2", "28_3",
                  "84_1", "84_2", "84_3")
  
  plot_df<-select(plot_df, "GN", "acc_site_ID", "igot_cluster", starts_with("p"))
  names(plot_df)[-(1:3)] <- time_point
  
  long_data <- plot_df %>%
    pivot_longer(
      cols = -c("GN", "acc_site_ID", "igot_cluster"),
      names_to = c("Time", "replicate"),
      names_sep = "_",
      values_to = "Score")%>%
    group_by(acc_site_ID, Time)%>%
    mutate(expression_pattern = mean(Score, na.rm = TRUE),
           expression_group = case_when(igot_cluster == 1 ~"decre",
                                        igot_cluster == 2 ~ "incre",
                                        is.na(igot_cluster) ~ "no_significant"),
           group_key = paste0(GN,"_" ,expression_group),
           Time = as.numeric(Time))%>%
    group_by(GN, Time, expression_group)%>%
    mutate(GN_mean_expression_pattern = mean(Score, na.rm = TRUE))%>%ungroup()
  
  long_data
  
}

igot_ave_plot<-function(df, GN_color){
  
  #dot = 各siteの動き、line = 各siteの動きの平均値
  plot_df <- df %>%select(GN, expression_group, group_key,  Time,  Score, GN_mean_expression_pattern)%>% distinct()
  
  p <-  ggplot(plot_df,aes(x = Time, y = Score, color = GN, group = expression_group)) +
    geom_line(data = df, aes(x = Time, y = GN_mean_expression_pattern, color = GN, group = group_key , linetype = expression_group), linewidth = 1.0) +
    scale_linetype_manual(values = c("decre" = "dotted", "incre" = "solid", "no_significant" = "solid")) +
    scale_alpha_manual(values = c("nosignificant" = 0.3, "decre" = 1.0, "incre" = 1.0))+
    geom_vline(xintercept = c(0, 8, 14, 28, 84),color = "dimgray", linewidth = 0.5, linetype = "dotted") +
    geom_point(size = 1, shape = 16, alpha = 0.6, position = position_jitter(width = 1, height = 0)) +
    labs(title = "" , x = "",  y = "scale abundance") +
    scale_x_continuous(breaks = c(0, 8, 14, 28, 84),
                       labels = c("P0", "P8", "2w", "1M", "3M")) +
    scale_color_manual(values = GN_color)  +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background  = element_rect(fill = "white", color = NA),  
      panel.grid = element_blank(), 
      axis.text.y  = element_text(color = "black", size = 20, face = "bold"),
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),
      axis.text.x  = element_text(color = "black",size = 20, face = "bold"),                                          
      axis.line.x = element_line(color = "black", size = 1),
      axis.line.y = element_line(color = "black", size = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.length.y = unit(-3, "pt"),
      legend.text  = element_text(size = 15))
  
  
  
  p}

trans_longdata_for_lineplot<-function(plot_df){
  
  time_point <- c("0_1", "0_2", "0_3",
                  "8_1", "8_2", "8_3",
                  "14_1", "14_2", "14_3",
                  "28_1", "28_2", "28_3",
                  "84_1", "84_2", "84_3")
  
  names(plot_df)[-(1:8)] <- time_point
  
  long_data <- plot_df %>%
    mutate(GN_site = str_replace(acc_site_ID, pattern = acc, GN))%>%
    pivot_longer(
      cols = -c("GN", "acc","GN_site", "acc_site_ID", "igot_cluster", "cluster_num_in_site", "dif_cluster", "site_num", "cluster_occupancy" ),
      names_to = c("Time", "replicate"),
      names_sep = "_",
      values_to = "Score")%>%
    mutate(group_key = ifelse(igot_cluster == 1, "decre", "incre"))%>%
    mutate(GN_group_key = paste0(GN, "_" ,group_key))%>%
    group_by(GN, group_key, Time)%>%
    mutate(GN_mean_expression_pattern = mean(Score, na.rm = TRUE),
           Time = as.numeric(Time))%>%ungroup()
  
  long_data
  
}

extract_GN<-function(ego_df, desc){
  
  ego_df<-overlap_24_GO
  desc<-target_high_GO
  
  target<-ego_df%>%
    filter(Description %in% desc)
  
  target_ID<-unlist(strsplit(target$geneID, "/"))
  
  target_GN <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = target_ID,
    keytype = "ENTREZID",
    columns = "SYMBOL") %>% pull(SYMBOL)
  
  target_GN
  
}

igot_GN_plot<-function(df, GN_color , target_GN ){
  
  df<-df%>%
    filter(GN %in% target_GN)
  
  #dot = 各マウスの動き、line = 各マウスの平均値
  
  plot_df <- df %>%
    select(GN_site, GN, group_key, GN_group_key,  Time,  Score, GN_mean_expression_pattern)%>% 
    distinct()
  
  n_site <- dplyr::n_distinct(plot_df$GN_site)
  
  cols <- if (is.function(GN_color)) GN_color(n_site) else GN_color
  
  
  p <-  ggplot(plot_df,aes(x = Time, y = Score, color = GN_site, group = group_key)) +
    geom_line(data = df, aes(x = Time, y = GN_mean_expression_pattern, color = GN_site, group = GN_group_key ), linewidth = 1.0) +
    geom_vline(xintercept = c(0, 8, 14, 28, 84),color = "dimgray", linewidth = 0.5, linetype = "dotted") +
    geom_point(size = 1, shape = 16, alpha = 0.6, position = position_jitter(width = 1, height = 0)) +
    labs(title = "" , x = "",  y = "scale abundance") +
    scale_x_continuous(breaks = c(0, 8, 14, 28, 84),
                       labels = c("P0", "P8", "2w", "1M", "3M")) +
    scale_color_manual(values = cols)  +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background  = element_rect(fill = "white", color = NA),  
      panel.grid = element_blank(), 
      axis.text.y  = element_text(color = "black", size = 20, face = "bold"),
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),
      axis.text.x  = element_text(color = "black",size = 20, face = "bold"),                                          
      axis.line.x = element_line(color = "black", size = 1),
      axis.line.y = element_line(color = "black", size = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.length.y = unit(-3, "pt"),
      legend.text  = element_text(size = 15))
  
  
  
  p}
