test_cor_for_site <- function(gpep_df, igot_df, fraction) {
  
  gpep_mat <- transform_input_zscore(gpep_df, fraction, "gpep")
  igot_mat <- transform_input_zscore(igot_df, fraction, "igot_cor")
  
  n <- nrow(gpep_mat)
  cor_vals <- p_vals <- c()
  
  for (i in 1:n) {
    x <- gpep_mat[i, ]
    y <- igot_mat[i, ]
    
    if (sum(complete.cases(x, y)) >= 3) {
      test <- suppressWarnings(cor.test(x, y, method = "pearson"))  
      cor_vals[i] <- test$estimate
      p_vals[i] <- test$p.value
    } else {
      cor_vals[i] <- NA
      p_vals[i] <- NA
    }
  }
  
  q_vals <- p.adjust(p_vals, method = "BH")
  cor_site <- data.frame(
    acc_site_glycan_ID = rownames(gpep_mat),
    gpep_igot_correlation = cor_vals,
    gpep_igot_cor_p_value = p_vals,
    gpep_igot_cor_q_value = q_vals
  )
  cor_site
}

find_sig_correlated_site <- function(cor_df, cut_off_column, cut_off) {
  
  
  sig_cor <- cor_df %>% 
    dplyr::select("acc_site_glycan_ID", 
                  all_of(cut_off_column), 
                  "gpep_igot_correlation") %>%
    filter( !is.na(.data[["gpep_igot_correlation"]]), 
            .data[[cut_off_column]] < cut_off, 
            .data[["gpep_igot_correlation"]]> 0) %>%
    mutate(significant_cor = 1) %>%
    dplyr::select("acc_site_glycan_ID", "significant_cor")
  
  sig_cor_df <- full_join(cor_df, sig_cor, by = "acc_site_glycan_ID")
  sig_cor_df
}


cluster_number_output_list <-
  function(df,
           gpep_cluster,
           igot_cluster,
           ratio_cluster,
           gpep_zscore,
           igot_zscore,
           ratio_zscore) {
    
  
  obtain_zscore_with_ID <- function(mat, mode){
    if(mode == "igot"){
      scale_df <- mat %>% data.frame() %>% rownames_to_column(var = "site_ID")
    } else {
      scale_df <- mat %>% data.frame() %>% rownames_to_column(var = "acc_site_glycan_ID")
    }
    rownames(scale_df) <- NULL
    scale_df
  }
  
  obtain_cluster_with_ID <- function(cluster, cluster_name, mode) {
    cluster_num_df <- cluster$cut %>% data.frame() %>% setNames(cluster_name)
    if(mode == "igot") {
      cluster_num_df$site_ID <- rownames(cluster_num_df)
    } else {
      cluster_num_df$acc_site_glycan_ID <- rownames(cluster_num_df)
    }
    rownames(cluster_num_df) <- NULL
    cluster_num_df
  }
  
  gpep_cluster  <- obtain_cluster_with_ID(gpep_cluster, "gpep_cluster", "gpep")
  igot_cluster  <- obtain_cluster_with_ID(igot_cluster, "igot_cluster", "igot")
  ratio_cluster <- obtain_cluster_with_ID(ratio_cluster, "ratio_cluster", "ratio")
  
  df_output <- df %>%
    full_join(gpep_cluster, by = "acc_site_glycan_ID" )%>%
    full_join(ratio_cluster, by = "acc_site_glycan_ID")%>%
    full_join(igot_cluster, by ="site_ID")
  
  scale_gpep_col <- c("First_acc","First_entry", "GN" ,
                      "Description","Glycosite_in_first_protein", "Nsite", 
                      "Nsite_seq", "Atypical", "Nsite_have",                
                      "GlycanComposition", "Hex",  "HexNAc",
                      "NeuAc", "NeuGc", "Fuc",
                      "HighMan", "CoreFucosed", "Bisecting",
                      "Category", "site_ID", "acc_site_glycan_ID",
                      "acc_category_ID", "glycan_ID")
  
  scale_igot_col <- c("First_acc","First_entry", "GN" ,
                      "Description","Glycosite_in_first_protein", "Nsite", 
                      "Nsite_seq", "Atypical", "Nsite_have",                
                      "site_ID")
  
  gpep_scale  <- obtain_zscore_with_ID(gpep_zscore, "gpep")
  igot_scale  <- obtain_zscore_with_ID(igot_zscore, "igot")
  ratio_scale <- obtain_zscore_with_ID(ratio_zscore, "ratio")
  
  gpep_scale_output <- select(df, all_of(scale_gpep_col)) %>%
    full_join(gpep_cluster, by = "acc_site_glycan_ID") %>%
    full_join(gpep_scale, by = "acc_site_glycan_ID")
  
  ratio_scale_output <- select(df, all_of(scale_gpep_col)) %>%
    full_join(ratio_cluster, by = "acc_site_glycan_ID") %>%
    full_join(ratio_scale, by = "acc_site_glycan_ID")
  
  igot_scale_output <- select(df, all_of(scale_igot_col)) %>%
    full_join(igot_cluster, by = "site_ID") %>%
    full_join(igot_scale, by = "site_ID") %>% distinct()
  
  list(
    "pattern_result" = df_output,
    "gpep_zscore" = gpep_scale_output,
    "ratio_zscore" = ratio_scale_output,
    "igot_zscore" = igot_scale_output
  )
  
}

