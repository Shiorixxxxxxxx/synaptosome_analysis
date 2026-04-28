select_column <- function(gpep) {
  
   gpep_common_col <- c(
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
  
  # add acc_category_ID column
  gpep_df <- gpep %>% 
    mutate(acc_category_ID = paste0(First_acc, "_", Nsite_seq, "_", Category), .after = "glycan_ID") %>% 
    select(all_of(gpep_common_col), matches("Ec[0-9]+$"))
  
  gpep_df
}

site_sum<-function(gpep){
  #acc_seq_site_Glycan_ID→acc_site_Glycan_IDへデータサマライズをここで実行するよう変更 
    gpep_common_col <- c(
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
  
gpep <- gpep %>%
  group_by(acc_site_glycan_ID) %>%
  mutate(across(contains("_Ec"), 
                list(sum = ~ sum(.)), 
                .names = "{.col}_site_sum")) %>%
  ungroup() %>% 
  select(all_of(gpep_common_col), matches("site_sum$"))%>%
  distinct()
colnames(gpep)<-str_replace(colnames(gpep), "_site_sum", "")  
gpep

}

log2_median_normal <- function(data) {
  
  data <- data[, grep("Ec", colnames(data), value = TRUE)]
  for(i in 1:ncol(data)) {
    # set 0 as NA
    data[[i]][data[[i]] == 0] <- NA
    # find median of each column
    col_median <- median(data[[i]], na.rm = TRUE)
    # calculate log2-median
    data[[i]] <- log2(data[[i]]) - log2(col_median)
  }
  colnames(data) <- paste0(colnames(data), "_log2m")
  data
}

rename_column <- function(data) {
  # change the column names
  colnames(data) <- gsub(".*_ctx_(.*)_(p0|p8|2w|4w|3M)-([1-3]).*_(Ec[0-9]+).*",
                         "\\1_\\2_\\3_\\4", colnames(data))
  data
} 



get_mean_column <- function(data){
  # group by acc_site_glycan_ID, and sum accrose the group
  
  time_order <- c("p0", "p8", "2w", "4w", "3M")
  shots <- 1:3
  # integrate 3 shots for each time point
  # RULE: 
  #  - if 2 non-missing v.s. 1 missing, take the mean of 2 non-missing values
  #  - if 1 non-missing v.s. 2 missing, take the non-missing value 
  for (tp in time_order) {
    for (rep in shots) {
      # get 3 shots
      pattern <- paste0("_", tp, "_", rep, "_Ec[0-9]+_log2m$")
      cols <- grep(pattern, names(data), value = TRUE)
      
      # set the integrated column name as "mean" and add to the current table
      new_col <- paste0("mean_", tp, "_", rep)
      tmp_func <- function(x) {
        valid_x <- x[!is.na(x)]
        if (length(valid_x) >= 2) {
          mean(valid_x, na.rm = TRUE)
        } else if (length(valid_x) == 1) {
          valid_x
        } else {
          NA
        }
      }
      data[[new_col]] <- apply(data[, cols, drop=FALSE], 1, tmp_func)
    }
  }
  data
}


##acc_site_glycan_ID_listの選択、>0を削除
select_rows <- function(data) {
 
  # filter 1: select the rows whose acc_site_glycan_ID has >= 3 counts
  data_group_by_site_ID <- data %>% group_by(site_ID)
  site_ID_list1 <- data_group_by_site_ID %>%
                   summarise(site_count = n_distinct(acc_site_glycan_ID)) %>% 
                   filter(site_count >= 3) %>% 
                   pull(site_ID) 
  
  
  data_sel <- data %>% filter(site_ID %in% site_ID_list1)
  
  
  # filter 2: select the rows that has at least two non-zero values 
  #           at any time point ("p0", "p8", "2w", "4w", "3M")
  agg_timepoint_cols <- list(
    "p0" = c("mean_p0_1", "mean_p0_2", "mean_p0_3"),
    "p8" = c("mean_p8_1", "mean_p8_2", "mean_p8_3"),
    "2w" = c("mean_2w_1", "mean_2w_2", "mean_2w_3"),
    "4w" = c("mean_4w_1", "mean_4w_2", "mean_4w_3"),
    "3M" = c("mean_3M_1", "mean_3M_2", "mean_3M_3")
  )
  
  acc_site_glycan_ID_list <- data_sel %>% 
    group_by(site_ID) %>%
    summarise(
      valid_time = any(sapply(agg_timepoint_cols, function(cols) {
        rowSums(!is.na(across(all_of(cols)))) >= 2
      })),
      .groups = "drop")
  
  acc_site_glycan_ID_list <-  acc_site_glycan_ID_list %>% 
    filter(valid_time) %>%
    pull(site_ID)
  
  data_sel <- data_sel %>% filter(site_ID %in% acc_site_glycan_ID_list) 
  
  # filter 3: remove all NA rows
  
  sample_cols<-c("mean_p0_1", "mean_p0_2", "mean_p0_3", 
                 "mean_p8_1", "mean_p8_2", "mean_p8_3",
                 "mean_2w_1", "mean_2w_2", "mean_2w_3", 
                 "mean_4w_1", "mean_4w_2", "mean_4w_3",
                 "mean_3M_1", "mean_3M_2","mean_3M_3")
  
  all_na <- rowSums(!is.na(data_sel[, sample_cols])) == 0
  
  data_sel <- data_sel[!all_na, ]
  
  data_sel
}

rule_based_imputation <- function(data) {
  data <- as.matrix(data)
  # find global minimum
  global_min <- min(data, na.rm = TRUE)
  
  unique_samples <- c("p0", "p8", "2w", "4w", "3M")
  for (sample in unique_samples) {
    cols <- grep(sample, colnames(data), value = TRUE)
    
    for(i in 1:nrow(data)) {
      n <- sum(!is.na(data[i,cols]))
      # n = 0: all missing, impute global_min 
      if (n == 0) {
        data[i,cols] <- global_min
      } else if(n == 2) {
        # n = 2: one missing, impute mean
        data[i,cols][is.na(data[i,cols])] <- mean(data[i,cols], na.rm = TRUE)  
      }
    }
  }
  data
}

PCA_imputation <- function(data) {
  # PCA imputation is only applied for after the rule_based_imputationo 
  pca_res <- pcaMethods::pca(t(data))
  imp_data <- t(pca_res@completeObs)
  colnames(imp_data) <- paste0(colnames(imp_data), "_imp")
  imp_data
}

merged_file <- function(gpep, igot){
 
  gpep_for_ratio <- gpep %>% dplyr::rename("acc_site_ID"="site_ID")%>%
                    select("acc_site_ID", contains("imp"))
  igot_for_ratio <- igot %>% select("acc_site_ID", contains("_log2m_imp"))
  colnames(igot_for_ratio) <- gsub("_Ec[0-9]+", "", colnames(igot_for_ratio))
  
  colnames(gpep_for_ratio)[-1] <- paste0(colnames(gpep_for_ratio)[-1], "_gpep")
  colnames(igot_for_ratio)[-1] <- paste0(colnames(igot_for_ratio)[-1], "_igot")
  
  # Merge two files
  merged_df <- left_join(
    gpep_for_ratio,
    igot_for_ratio,
    by = "acc_site_ID")
  
  # Calculate ratio
  samples  <- c(
    "p0_1", "p0_2", "p0_3", 
    "p8_1", "p8_2", "p8_3",
    "2w_1", "2w_2", "2w_3",
    "4w_1", "4w_2", "4w_3",
    "3M_1", "3M_2", "3M_3"
  )
  for(sample in samples) {
    
    col_ratio <- paste0("ratio_", sample)
    idx_gpep <- grep(paste0(sample, ".*gpep"), colnames(merged_df), value = TRUE)
    idx_igot <- grep(paste0(sample, ".*igot"), colnames(merged_df), value = TRUE)
    
    merged_df[[col_ratio]] <- merged_df[[idx_gpep]] - merged_df[[idx_igot]]
    
  }
  
  merged_df
}


output_merge_file<-function(gpep, igot, ratio_df){
  
  ratio_df<-bind_cols(gpep, select(ratio_df,(matches("_igot")), (matches("ratio"))))
  
  output_temp<-list("igot"=igot, "gpep"=gpep, "ratio"=ratio_df)
  output_temp
}


