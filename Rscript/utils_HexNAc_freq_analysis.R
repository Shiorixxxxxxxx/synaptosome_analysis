Glycan_class<-function(gpep_df, mode){
  
  if(mode == "Branch_HexNAc"){
    
    HN_class_df<-gpep_df%>%
      filter(Category != "Pauci", HexNAc>=3)%>%
      mutate(Branch_HexNAc = as.numeric(Branch_HexNAc)-2,
             pattern = case_when(ratio_cluster == 1 ~ "Decreased",
                                 ratio_cluster == 2 ~ "Increased & V-shaped",
                                 ratio_cluster == 3 ~ "Increased & V-shaped",
                                 ratio_cluster == NA ~ "no_cluster"))
  }else if(mode == "TermFuc"){
    
    HN_class_df<-gpep_df%>%
      filter(Category != "Pauci", HexNAc>=3)%>%
      mutate(TermFuc = as.numeric(TermFuc),
             pattern = case_when(ratio_cluster == 1 ~ "Decreased",
                                 ratio_cluster == 2 ~ "Increased & V-shaped",
                                 ratio_cluster == 3 ~ "Increased & V-shaped",
                                 ratio_cluster == NA ~ "no_cluster"))
  }else if(mode == "Bisecting"){
    
    HN_class_df<-gpep_df%>%
      filter(!Category %in% c("HM","Pauci"), HexNAc>=3)%>%
      mutate(all_Bisect = ifelse(Bisecting == 1 | CF_Bisect == 1, 1, 0),
             pattern = case_when(ratio_cluster == 1 ~ "Decreased",
                                 ratio_cluster == 2 ~ "Increased & V-shaped",
                                 ratio_cluster == 3 ~ "Increased & V-shaped",
                                 ratio_cluster == NA ~ "no_cluster"))
    
  }else if(mode =="CoreFucose"){
    
    HN_class_df<-gpep_df%>%
      filter(!Category %in% c("HM","Pauci"), HexNAc>=3)%>%
      mutate(all_CF = ifelse(CoreFucose == 1 | CF_Bisect == 1, 1, 0),
             pattern = case_when(ratio_cluster == 1 ~ "Decreased",
                                 ratio_cluster == 2 ~ "Increased & V-shaped",
                                 ratio_cluster == 3 ~ "Increased & V-shaped",
                                 ratio_cluster == NA ~ "no_cluster"))
  }
  HN_class_df
  
}

plot_data_make<-function(df, col, all_mode, count_col){
  #col = acc_site_ID, GN, pattern...
  
  if(all_mode == FALSE){
    count_HN_df<-df%>%
      filter(!is.na(pattern))%>%
      select(all_of(col), GN, all_of(count_col))
    #mutate(GN = factor(GN, level = Aliment_level))
  }else{
    count_HN_df<-df%>%
      select(all_of(col), GN, all_of(count_col))
  }
  count_HN_df
}

freq_cal<-function(df, count_col){
  
  freq_df <- df %>%
    filter(!is.na(pattern))%>%
    count(GN, !!sym(count_col), pattern)%>%
    group_by(GN, pattern) %>%
    mutate(freq = n / sum(n)) %>%       
    ungroup()
  
  freq_df
  
}

glycan_dif_analysis<-function(freq_df, Glycan){
  
  
  Glycan_change<-freq_df%>%
    filter(freq >= 0.4)%>%
    group_by(GN) %>%
    mutate(n_rows = n())%>%
    group_by(GN, !!sym(Glycan))%>%
    slice_max(freq, n = 1, with_ties = TRUE)%>%ungroup()
  
  extract_Glycan_change_GN<- Glycan_change%>%
    group_by(GN, pattern)%>%
    mutate(average_num = mean(!!sym(Glycan), na.rm = FALSE))%>%
    slice_max(average_num, n = 1, with_ties = FALSE)%>%
    group_by(GN)%>%
    mutate(new_nrows = n())%>%
    filter(n_rows != 1, new_nrows != 1)%>%
    mutate(target_glyco_num = ifelse(pattern == "Decreased", -average_num, average_num))%>%
    group_by(GN)%>%
    mutate(dif_num = sum(target_glyco_num))%>%
    ungroup()
  
  extract_Glycan_change_GN
}

pairwise_res<-function(pattern_df, count_col, pattern_pair){
  
  if(pattern_pair == FALSE){
    pairwise_res_df <- pattern_df %>%
      pairwise_t_test(formula = as.formula(paste(count_col, "~ GN")), p.adjust.method = "BH", detailed = T)
    
    means <- pattern_df %>%
      group_by(GN) %>%
      summarise(mean_val = mean(.data[[count_col]], na.rm = TRUE))
    
    pairwise_res_df <- pairwise_res_df %>%
      left_join(means, by = c("group1" = "GN")) %>%
      rename(mean1 = mean_val) %>%
      left_join(means, by = c("group2" = "GN")) %>%
      rename(mean2 = mean_val) %>%
      mutate(
        estimate = mean1 - mean2,
        direction = ifelse(estimate > 0,
                           paste(group1, ">", group2),
                           paste(group1, "<", group2)))
    
  }else{
    pairwise_res_df <- pattern_df %>%
      group_by(GN)%>%
      pairwise_t_test(formula = as.formula(paste(count_col, "~ pattern")), p.adjust.method = "BH", detailed = T)
    
    means <- pattern_df %>%
      group_by(GN, pattern) %>%
      summarise(mean_val = mean(.data[[count_col]], na.rm = TRUE),.groups = "drop")
    
    pairwise_res_df <- pairwise_res_df %>%
      left_join(means, by = c("GN", "group1" = "pattern")) %>%
      rename(mean1 = mean_val) %>%
      left_join(means, by = c("GN", "group2" = "pattern")) %>%
      rename(mean2 = mean_val) %>%
      mutate(
        estimate = mean1 - mean2,
        direction = ifelse(estimate > 0,
                           paste(group1, ">", group2),
                           paste(group1, "<", group2)))
  }
  pairwise_res_df
}

count_heatmap<-function(df, count_glycan, group, all_mode){
  
  if(all_mode == FALSE){
    p<-ggplot(df, aes(x = !!sym(count_glycan) , y = !!sym(group) , fill = freq)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "plasma") +
      facet_wrap(~ pattern, strip.position = "top") +
      theme_classic() +
      theme(
        strip.placement = "outside",    
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background  = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(), 
        axis.text.y  = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        axis.text.x  = element_text(color = "black",size = 20, face = "bold"),                               
        axis.ticks.x = element_blank())
  }else{
    p<-ggplot(df, aes(x = !!sym(count_glycan) , y = !!sym(group) , fill = freq)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "plasma") +
      theme_classic() +
      theme(
        strip.placement = "outside",    
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold", color = "black"),
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background  = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(), 
        axis.text.y  = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        axis.text.x  = element_text(color = "black",size = 20, face = "bold"),                               
        axis.ticks.x = element_blank())
  }
  p
  
}

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
