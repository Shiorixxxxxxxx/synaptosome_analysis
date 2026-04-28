venn_analysis<-function(list){
  venn<-venndir(list,
                proportional=TRUE,
                set_colors = c( "#ff6347", "#6495ed", "#ffff00"),
                show_labels="cs",
                show_segments=FALSE,
                inside_percent_threshold=0,
                overlap_type="overlap",
                font_cex = c(1.5, 1.5, 0.8))
  venn}

extract_items <- function(label_df, target_label) {
  key<-label_df %>%
    dplyr::filter(overlap_set == target_label) %>%
    dplyr::pull(items) %>%
    unlist()
  
  key
}

trans_longdata<-function(plot_df){
  
  time_point <- c("0_1", "0_2", "0_3",
                  "8_1", "8_2", "8_3",
                  "14_1", "14_2", "14_3",
                  "28_1", "28_2", "28_3",
                  "84_1", "84_2", "84_3")
  
  names(plot_df)[-(1:3)] <- time_point
  
  long_data <- plot_df %>%
    pivot_longer(
      cols = -c("GN", "site_ID", "igot_cluster"),
      names_to = c("Time", "replicate"),
      names_sep = "_",
      values_to = "Score")%>%
    group_by(site_ID, Time)%>%
    mutate(expression_pattern = mean(Score, na.rm = TRUE),
           Time = as.numeric(Time))%>%
    group_by(GN, Time)%>%
    mutate(GN_mean_expression_pattern = mean(Score, na.rm = TRUE))%>%ungroup()
  
  long_data
  
}

trans_longdata_gpep<-function(plot_df, target_glycan_ID){
  
  
  time_point <- c("0_1", "0_2", "0_3",
                  "8_1", "8_2", "8_3",
                  "14_1", "14_2", "14_3",
                  "28_1", "28_2", "28_3",
                  "84_1", "84_2", "84_3")
  
  names(plot_df)[-(1:5)] <- time_point
  
  long_data <- plot_df %>%
    pivot_longer(
      cols = -c("GN", "site_ID", "glycan_ID","acc_site_glycan_ID", any_of(c("ratio_cluster", "gpep_cluster"))),
      names_to = c("Time", "replicate"),
      names_sep = "_",
      values_to = "Score")%>%
    mutate(group_key = ifelse(glycan_ID %in% target_glycan_ID, "decre", "incre"))%>%
    mutate(GN_group_key = paste0(GN, "_" ,group_key))%>%
    group_by(GN, group_key, Time)%>%
    mutate(GN_mean_expression_pattern = mean(Score, na.rm = TRUE),
           Time = as.numeric(Time))%>%ungroup()
  
  long_data
  
}

igot_ave_plot<-function(df, GN_color){
  
  #dot = 各siteの動き、line = 各siteの動きの平均値
  plot_df <- df %>%select(GN, site_ID, Time, igot_cluster, expression_pattern, GN_mean_expression_pattern)%>% distinct()
  
  p <-  ggplot(plot_df,
               aes(x = Time, y = expression_pattern, color = GN, group = GN)) +
    geom_line(data = df, aes(x = Time, y = GN_mean_expression_pattern, color = GN), linewidth = 1.0) +
    geom_vline(xintercept = c(0, 8, 14, 28, 84),color = "dimgray", linewidth = 0.5, linetype = "dotted") +
    geom_point(size = 1, shape = 16, alpha = 0.6, position = position_jitter(width = 2, height = 0)) +
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

gpep_ave_plot<-function(df, GN_color){
  
  #dot = 各マウスの動き、line = 各マウスの平均値
  
  plot_df <- df %>%select(GN, group_key, GN_group_key,  Time,  Score, GN_mean_expression_pattern)%>% distinct()
  
  n_GN <- dplyr::n_distinct(plot_df$GN)
  
  cols <- if (is.function(GN_color)) GN_color(n_GN) else GN_color
  
  p <-  ggplot(plot_df,aes(x = Time, y = Score, color = GN, group = group_key)) +
    geom_line(data = df, aes(x = Time, y = GN_mean_expression_pattern, color = GN, group = GN_group_key , linetype = group_key), linewidth = 1.0) +
    scale_linetype_manual(values = c("decre" = "dotted", "incre" = "solid")) +
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