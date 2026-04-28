select_cols <- function(df) {
  # select the corresponding columns
  cols <- c(
    "acc_site_glycan_ID",
    "GN",
    "gpep_igot_cor_p_value",
    "significant_cor",
    "gpep_cluster",
    "ratio_cluster",
    "igot_cluster"
  )
  df %>% filter(!is.na(gpep_igot_correlation))%>% select(all_of(cols))
}

add_glycan_structural_info <- function(df) {
  
  # extract acc / site / glycan_comp / category from rownames
  extract_glycan_info <- function(rowname) {
    parts <- strsplit(rowname, "_")[[1]]
    list(acc = parts[1],
         site = parts[2],
         glycan_comp = parts[3],
         Category = paste(parts[4:length(parts)], collapse = "_"),
         acc_site_glycan_ID = rowname)
  }
  glycan_list <- lapply(df$acc_site_glycan_ID, extract_glycan_info) %>% bind_rows()
  
  # find and add glycan category info as columns
  glycan_list <-  glycan_list %>% 
    mutate(
    acc_site_ID = paste0(glycan_list$acc,  "_", glycan_list$site),
    GlycanComposition_ID = paste0(glycan_list$glycan_comp, "_", glycan_list$Category),
    # (bool) catergory info
    HM = if_else(glycan_list$Category == "HM", 1, 0),
    CoreFucose = if_else(glycan_list$Category == "CF_noBisect", 1, 0),
    Bisecting = if_else(glycan_list$Category == "noCF_Bisect", 1, 0),
    CF_Bisect = if_else(glycan_list$Category == "CF_Bisect", 1, 0),
    Pauci = if_else(glycan_list$Category == "Pauci", 1, 0),
    # (integer) composition info
    Hex = str_extract(glycan_list$glycan_comp, "(?<=Hex\\()\\d+") %>% 
      as.numeric() %>% replace_na(replace = 0),
    HexNAc = str_extract(glycan_list$glycan_comp, "(?<=HexNAc\\()\\d+") %>% 
      as.numeric() %>% replace_na(replace = 0),
    Fuc = str_extract(glycan_list$glycan_comp, "(?<=Fuc\\()\\d+") %>% 
      as.numeric() %>% replace_na(replace = 0),
    NeuAc = str_extract(glycan_list$glycan_comp, "(?<=NeuAc\\()\\d+") %>% 
      as.numeric() %>% replace_na(replace = 0),
    NeuGc = str_extract(glycan_list$glycan_comp, "(?<=NeuGc\\()\\d+") %>% 
      as.numeric() %>% replace_na(replace = 0), 
    # (integer) terminal fucose
    TermFuc = if_else(CoreFucose == 1| CF_Bisect ==1 , Fuc - 1, Fuc),
    Branch_HexNAc = if_else(Bisecting == 1| CF_Bisect ==1, pmax(HexNAc - 1, 0), HexNAc),
    # (bool) if the glycan has both sialic acid and terminal fucose
    has_sialic_fuc = if_else(NeuAc > 0 | NeuGc > 0 | TermFuc > 0, 1, 0) %>% replace_na(replace = 0),
    has_sia = if_else(NeuAc > 0 | NeuGc > 0, 1, 0) %>% replace_na(replace = 0),
    # (integer) sialic acid category fucose
    fuc_sia_category = case_when(TermFuc > 0 & has_sia == 0 ~ 1,
                                 TermFuc > 0 & has_sia > 0 ~ 2,
                                 TermFuc == 0 & has_sia > 0 ~ 3,
                                 .default = 4)
    ) 
  
  df_res <- full_join(df, glycan_list, by = "acc_site_glycan_ID")
  df_res
}


add_cor_ratio_by_site <- function(df) {
  # find the correlated glyco-peptide (igot - gpep) ratio for each site 
  df_summary <- df %>%
    group_by(acc_site_ID) %>%
    summarise(
      n_total = sum(!is.na(gpep_igot_cor_p_value)),
      n_sig = sum(significant_cor == 1, na.rm = TRUE),
      ratio_sig = n_sig / n_total
    )
  
  df_res <- full_join(df, df_summary, by = "acc_site_ID")
  df_res
}

select_rows <- function(df,
                        cor ,
                        cut_off ,
                        data_level ,
                        cluster ,
                        sig_Glycan_least_mode ) {
  
    # cor = TRUE: find correlated site
    # cor = FALSE: find uncorrelated site
    # cluster = c("gpep_cluster", "ratio_cluster", NULL) 
    # cut_off: correlated gpep ratio for each site (default from 0.4~0.6)
    
    if (cor) {
      sel_key <- df %>% filter(ratio_sig > cut_off) %>%
        pull(.data[[data_level]]) %>%
        unique()  
    } else {
      sel_key <- df %>% filter(ratio_sig <= cut_off) %>%
        pull(.data[[data_level]]) %>%
        unique()  
    }
    df_sel <- df %>% filter(.data[[data_level]] %in% sel_key)
    
    # filter site 
    # selecting sites that have at least one significantly correlated pair (gpep, igot)
    if (sig_Glycan_least_mode) {
      sig_site <- df %>% 
        filter(is.na(significant_cor)) %>%
        group_by(acc_site_ID) %>%
        mutate(n_clustered_glycans = sum(!is.na(.data[[cluster]]))) %>%
        ungroup() %>%
        filter(n_clustered_glycans >= 1) %>%
        pull(.data[[data_level]]) %>%
        unique()
      
      df_sel <- df_sel %>% filter(.data[[data_level]] %in% sig_site)
    }
    df_sel
}

GO_term <- function(df, GO_mode) {
  
  GO_enrichment <- function(sel_df) {
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
    
    # 上位GO term取得
    top_terms <- ego_df %>%
      arrange(p.adjust) %>%
      dplyr::slice(1:20) %>%
      pull(ID)  # または pull(Description)
    
    # p.adjust付きの gene-term ペアに展開
    gene2go_full <- ego_df %>%
      filter(ID %in% top_terms) %>%
      dplyr::select(ID, Description, p.adjust, geneID) %>%
      separate_rows(geneID, sep = "/")
    
    # 各geneIDごとに、最も小さいp.adjustのGO termを1つだけ選ぶ
    best_gene2go <- gene2go_full %>%
      group_by(geneID) %>%
      arrange(p.adjust) %>%
      dplyr::slice(1) %>%
      ungroup()
    unique(best_gene2go$Description)
    
    # 5. タンパク質ごとに機能グループ化
    protein_function_df <- entrez_ids %>%
      left_join(best_gene2go, by = c("ENTREZID" = "geneID")) %>%
      mutate(function_group = if_else(is.na(ID), "Other", ID)) %>%
      dplyr::rename("GN" = "SYMBOL") %>%
      dplyr::select(GN, Description, function_group, p.adjust)
    
  }
  
  #igotの変動ごとにGO enrichment
  igot_increased_site <- filter(df, igot_cluster == 1)
  igot_decreased_site <- filter(df, igot_cluster == 2)
  igot_nochanged_site <- filter(df, igot_cluster != 1 | 2)
  
  if (length(igot_increased_site) != 0) {
    igot_increase_GO <- GO_enrichment(igot_increased_site)
  }
  if (length(igot_decreased_site) != 0) {
    igot_decrease_GO <- GO_enrichment(igot_decreased_site)
  }
  if (length(igot_nochanged_site) != 0) {
    igot_nochanged_GO <- GO_enrichment(igot_nochanged_site)
  }
  
  if (length(igot_increased_site) != 0) {
    igot_nochanged_GO <- bind_rows(igot_nochanged_GO, igot_increase_GO)
  }
  if (length(igot_decreased_site) != 0) {
    igot_nochanged_GO <- bind_rows(igot_nochanged_GO, igot_decrease_GO)
  }
  
  #重複した場合、p.adjastが低い方を優先
  GO_selected <- igot_nochanged_GO %>%
    group_by(GN) %>%
    slice_min(order_by = if_else(is.na(p.adjust), Inf, p.adjust),
              with_ties = FALSE) %>%
    ungroup()
  
  #元データに適応
  df <- left_join(df, GO_selected, by = "GN")
  
  df
}

sort_key <- function(df, Alignment) {
  
  if (Alignment == "Category") {
    acc_sort <- df %>%
      mutate(
        major_group = case_when(
          HM == 1 ~ 1,
          Pauci == 1 ~ 2,
          CoreFucose == 1 ~ 3,
          CF_Bisect == 1 ~ 4,
          Bisecting == 1 ~ 5,
          .default = 6
        ),
        sort_key = paste0(
          major_group,
          "_",
          sprintf("%02d", Branch_HexNAc),
          "_",
          fuc_sia_category,
          "_",
          sprintf("%d", NeuAc),
          "_",
          sprintf("%d", NeuGc),
          "_",
          sprintf("%02d", TermFuc),
          "_",
          sprintf("%02d", Hex)
        )
      )
  }
  
  if (Alignment == "HexNAc") {
    acc_sort <- df %>%
      mutate(
        major_group = case_when(HM == 1 ~ 1,
                                Pauci == 1 ~ 2,
                                TRUE ~ 3),
        
        category_group = case_when(CoreFucose == 1 ~ 1,
                                   CF_Bisect == 1 ~ 2,
                                   Bisecting == 1 ~ 3,
                                   .default = 4),
        
        sort_key = paste0(
          sprintf("%d", major_group),
          "_",
          sprintf("%02d", Branch_HexNAc),
          "_",
          has_sia,
          "_",
          sprintf("%d", category_group),
          "_",
          sprintf("%d", fuc_sia_category),
          "_",
          sprintf("%d", NeuAc),
          "_",
          sprintf("%d", NeuGc),
          "_",
          sprintf("%02d", TermFuc),
          "_",
          sprintf("%02d", Hex)
        )
      )
  }
  
  if (Alignment == "Sia") {
    acc_sort <- df %>%
      mutate(
        major_group = case_when(HM == 1 ~ 1,
                                Pauci == 1 ~ 2,
                                .default = 3),
        
        category_group = case_when(CoreFucose == 1 ~ 1,
                                   CF_Bisect == 1 ~ 2,
                                   Bisecting == 1 ~ 3,
                                   .default = 4),
        
        sort_key = paste0(
          major_group,
          "_",
          has_sia,
          "_",
          sprintf("%02d", Branch_HexNAc),
          "_",
          category_group,
          "_",
          fuc_sia_category,
          "_",
          sprintf("%d", NeuAc),
          "_",
          sprintf("%d", NeuGc),
          "_",
          sprintf("%02d", TermFuc),
          "_",
          sprintf("%02d", Hex)
        )
      )
  }
  
  if (Alignment == "Fuc") {
    acc_sort <- acc_sort %>%
      mutate(
        major_group = case_when(HM == 1 ~ 1,
                                Pauci == 1 ~ 2,
                                .default = 3),
        
        category_group = case_when(CoreFucose == 1 ~ 1,
                                   CF_Bisect == 1 ~ 2,
                                   Bisecting == 1 ~ 3,
                                   .default = 4),
        
        sort_key = paste0(
          major_group,
          "_",
          sprintf("%02d", TermFuc),
          "_",
          has_sia,
          "_",
          fuc_sia_category,
          "_",
          sprintf("%d", NeuAc),
          "_",
          sprintf("%d", NeuGc),
          "_",
          sprintf("%02d", Branch_HexNAc),
          "_",
          category_group,
          "_",
          sprintf("%02d", Hex)
        )
      )
  }
  
  acc_sort
}

extract_GN<-function(ego_df, desc){
  
  target<-ego_df%>%
    filter(ID %in% desc)
  
  target_ID<-unlist(strsplit(target$geneID, "/"))
  
  target_GN <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = target_ID,
    keytype = "ENTREZID",
    columns = "SYMBOL") %>% pull(SYMBOL)
  
  target_GN
  
}

heatmaply_fig <- function(acc_sort, cluster) {
  # obtain the row and column names
  mat_col_names <- acc_sort %>%　
    arrange(sort_key) %>%　
    pull(GlycanComposition_ID) %>%　
    unique()
  mat_row_names <- acc_sort %>% 
    arrange(function_group, acc, igot_cluster) %>% 
    pull(acc_site_ID) %>% 
    unique()
  # first construct the matrix
  n <- length(mat_row_names)
  p <- length(mat_col_names)
  mat <- matrix(NA, nrow = n, ncol = p)
  rownames(mat) <- mat_row_names
  colnames(mat) <- mat_col_names
  # fill the matrices
  for(i in 1:n) for(j in 1:p) {
    acc_site_id <- mat_row_names[i]
    glycan_comp_id <- mat_col_names[j]
    cond1 <- acc_sort$GlycanComposition_ID == glycan_comp_id
    cond2 <- acc_sort$acc_site_ID == acc_site_id
    val <- acc_sort[[cluster]][cond1 & cond2]
    # add to gpep_mat
    if(length(val) != 0) {
      if(is.na(val)) {
        mat[i, j] <- 0
      } else if(val == 2) {
        mat[i, j] <- 2
      } else if(val == 1) {
        mat[i, j] <- 1
      }else if(val == 3) {          #ratio_clusterに備えて3の場合を作成
        mat[i, j] <- 3
      }
    }
  }
  # cluster 1: "#4662D7FF"
  # cluster 2: "#ff3366"
  # cluster 3: "#ffa500"
  # NA: "#f8f8ff"
  if(cluster == "gpep_cluster"){
  colors <- c("#e2e8f8",  "#4169e1", "#ff3366")}
  
  if(cluster == "ratio_cluster"){
    colors <- c("#e2e8f8", "#4169e1", "#ff3366", "#ffa500")}  #NA からスタートするので色の並びを修正
  
  # annotation row color bars
  acc_sort_row <- acc_sort %>% arrange(function_group, acc, igot_cluster)
  ind <- match(rownames(mat), acc_sort_row$acc_site_ID)
  acc_list <- unique(acc_sort_row$acc)
  acc_palette <- brewer.pal(12, "Set3")
  acc_color_map <- setNames(acc_palette[(seq_along(acc_list) - 1) %% length(acc_palette) + 1], acc_list)
  acc_colors <- acc_color_map[acc_sort_row$acc[ind]]
  func_list <- unique(acc_sort_row$function_group)
  func_palette <- brewer.pal(min(length(func_list), 12), "Paired")
  func_color_map <- setNames(func_palette[(seq_along(func_list) - 1) %% length(func_palette) + 1], func_list)
  func_colors <- func_color_map[acc_sort_row$function_group[ind]]
  igot_cluster <- acc_sort_row$igot_cluster
  igot_cluster[is.na(igot_cluster)] <- "nocluster"
  igot_list <- unique(igot_cluster)
  igot_color_map <- c("1" = "#4169e1", "2" = "#ff3366", "nocluster" = "#f5f5f5")
  igot_colors <- igot_color_map[igot_cluster[ind]]
  
  row_side_colors <- data.frame(Acc = acc_colors,
                                Function = func_colors,
                                IGOT = igot_colors)
  rownames(row_side_colors) <- rownames(mat)
  # annotation col color bars
  acc_sort_col <- acc_sort %>% arrange(sort_key)
  ind <- match(colnames(mat), acc_sort_col$GlycanComposition_ID)
  category_list <- unique(acc_sort_col$Category)
  category_color_map <- c(
    "HM" = "#40e0d0",
    "CF_noBisect" = "#ff6347",
    "CF_Bisect" = "#ffa500",
    "noCF_noBisect" = "#3cb371",
    "Pauci" = "#e0ffff",
    "noCF_Bisect" = "#6495ed",
    "nocategory" = "#f8f8ff"
  )
  category_colors <- category_color_map[acc_sort_col$Category[ind]]
  col_side_colors <- category_colors
  # show the heatmaply
  p<-heatmaply(x = mat, Rowv = NA, Colv = NA, 
            colors = colors, hide_colorbar = FALSE, 
            RowSideColors = row_side_colors, 
            ColSideColors = col_side_colors)
  
  p <- p %>% layout(
    xaxis = list(showticklabels = TRUE),
    yaxis = list(showticklabels = FALSE),
    showlegend = FALSE)
  p
}

#張先生修正後削除
data_cor_select<-function(df,cluster, cor, cut_off, data_level, sig_Glycan_least_mode, num_Glycan){
  #cor =TRUE; 相関するSite→使用するClusterはgpep_cluster
  #cor =FALSE; gpepとigotの動きが非相関。→使用するClusterはratio_cluster
  #cuto_offは1siteあたりの相関するgpepの割合。0.4~0.6ぐらいが適当
  
  sel_key<- df %>%
    filter(if(cor) ratio_sig > cut_off else ratio_sig <= cut_off)%>%
    pull(.data[[data_level]])%>%
    unique()
  
  df_sel<-df %>%
    filter(.data[[data_level]] %in% sel_key)
  
  #acc_SIte上にGlycanが少なくとも指定数あるSiteのみ抽出 
  if(sig_Glycan_least_mode == TRUE){
    sig_site <- df %>%
      group_by(acc_site_ID) %>%
      mutate(n_clustered_glycans = sum(!is.na(.data[[cluster]]))) %>%
      ungroup()%>%
      filter(n_clustered_glycans >= num_Glycan)%>%
      pull(.data[[data_level]])%>%
      unique()
    
    df_sel<-df_sel%>%
      filter(.data[[data_level]] %in% sig_site)}
  
  df_sel}

venn_analysis<-function(df_cor, df_noncor, comp_ID){
  
  cor_ID<-pull(df_cor, comp_ID)%>%unique()
  noncor_ID<-pull(df_noncor, comp_ID)%>%unique()
  comp_list<-list("correlation" = cor_ID, "non_correlation"=noncor_ID)
  
  venn<-venndir(comp_list,
                proportional=TRUE,
                set_colors = c("#7fffd4", "#ffb6c1"),
                show_labels="cs",
                show_segments=FALSE,
                inside_percent_threshold=0,
                overlap_type="overlap",
                font_cex = c(1.5, 1.5, 0.8))
  venn
}

extract_items <- function(label_df, target_label) {
  key<-label_df %>%
    dplyr::filter(venn_label == target_label) %>%
    dplyr::pull(items) %>%
    unlist()
  
  if(target_label == "non_correlation"){
    total_both_acc <- label_df %>%
      dplyr::filter(is.na(venn_label) & overlap_set == "correlation&non_correlation") %>%
      dplyr::pull(items) %>%
      unlist()
    
    key<-c(total_both_acc, key)
  }
  key
}