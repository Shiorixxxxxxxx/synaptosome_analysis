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

get_uniprot_features <- function(acc){
  
  acc <- sub("-.*$", "", acc)
  url <- paste0("https://rest.uniprot.org/uniprotkb/", acc, ".json")
  r <- httr::GET(url)
  
  sc <- httr::status_code(r)
  
  # --- 404 はここでスキップ ---
  if (sc == 404) {
    message("404 Not found: ", acc)
    return(dplyr::tibble())
  }
  
  # --- それ以外で 200 じゃなければ、とりあえず飛ばす ---
  if (sc != 200) {
    message("HTTP ", sc, " for ", acc, " → スキップします")
    return(dplyr::tibble())
  }
  
  # ここまで来たら OK なレスポンス
  jj <- jsonlite::fromJSON(
    httr::content(r, "text", encoding = "UTF-8"),
    flatten = TRUE
  )
  
  feats <- jj$features
  
  if (is.null(feats) || nrow(feats) == 0) {
    message("No features for: ", acc)
    return(dplyr::tibble())
  }
  
  feats <- dplyr::as_tibble(feats) %>%
    dplyr::mutate(
      acc   = acc,
      start = as.numeric(location.start.value),
      end   = as.numeric(location.end.value)
    ) %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::select(acc, type, description, start, end)
  
  feats
}


get_feature<-function(data){
  
  feature<-NULL
  
  for(num in 1:length(unique(data$acc))){
    
    acc<-unique(data$acc)[num]
    
    acc_feature<-get_uniprot_features(acc)
    
    feature<-bind_rows(feature, acc_feature)
    
  }
  
  feature
}


domain_enrichment<-function(sample_gn, back_ground){
  
  uniprot_feature<-get_feature(back_ground)
  
  target_acc <- filter(back_ground, GN %in% sample_gn)%>%
    pull(acc)%>%unique()
  
  uniprot_feature <- uniprot_feature %>%
    mutate(domain = if_else(!is.na(description) & description != "",description,type))%>%
    distinct(acc, domain, .keep_all = TRUE)
  
  target_feature<-filter(uniprot_feature, acc %in% target_acc)
  
  bg_proteins <- unique(uniprot_feature$acc)
  
  N<-length(bg_proteins)
  n<-length(target_acc)
  
  domain_enrich <- uniprot_feature %>%
    
    distinct(acc, domain) %>%
    group_by(domain) %>%
    summarise(
      all = n_distinct(acc),                               # 背景でそのドメインを持つタンパク質数
      count = n_distinct(acc[acc %in% target_acc]),     # ターゲットでそのドメインを持つタンパク質数
      .groups = "drop"
    ) %>%
    filter(count > 0) %>%
    mutate(
      N = N,
      n = n,
      p_value = phyper(count - 1, all, N - all, n, lower.tail = FALSE),
      p_adj   = p.adjust(p_value, method = "BH"),
      bg_freq     = all / N,           # 背景におけるドメイン頻度
      target_freq = count / n,           # ターゲットにおけるドメイン頻度
      enrichment  = target_freq / bg_freq  # enrichment 値（>1 で多い）
    ) %>%
    arrange(p_adj, p_value)
  
  list<-list("res"= domain_enrich, "all_feature" =uniprot_feature, "target_feature"=target_feature)
  
  list
  
}

domain_enrich_plot<-function(enrichment_res){
  sig_domains <- enrichment_res %>%
    filter(
      all >= 3,        # 背景であまりにレアなドメインは除外（必要に応じて調整）
      count >= 2,      # ターゲットでも2つ以上出ているものに絞るなど
      p_adj < 0.05     # FDR 5%
    )%>%
    mutate(log10_p = -log10(p_adj))%>%
    arrange(log10_p)%>%
    mutate(domain = factor(domain, levels = unique(domain)))
  
  p<-ggplot(sig_domains, aes(x = log10_p, y = domain, size = count, color = p_adj)) +
    geom_point() +
    scale_color_gradient(name = "p_adjust", low = "#1e90ff", high = "#ff6347", trans = "reverse") +
    labs(
      title = paste0("Domain enrichment"),
      x = "-log10(p_adjust)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text.y  = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x  = element_text(color = "black",size = 10, face = "bold"),                                          
      axis.ticks.y = element_line(color = "black", size = 1))
  p}

domain_plot<-function(gpep, target_domain, plot_class, label_site, gn_order){
  
  gpep<-gpep%>%
    mutate(site_label = acc_site_ID %in% label_site)
  
  #最初のtransmembrane領域を抽出
  tm_anchor <- target_domain %>%
    filter(type %in% c("Transmembrane")) %>% 
    arrange(acc, start) %>%
    group_by(acc) %>%
    slice(1) %>%                               
    ungroup() %>%
    transmute(acc, type, dm1_start = start, dm1_end = end)
  
  #ALL
  ALL <- target_domain %>%
    arrange(acc, start) %>%
    group_by(acc) %>%
    slice(1) %>%                               
    ungroup() %>%
    transmute(acc, type, dm1_start = start, dm1_end = end)
  
  
  #GPIを抽出
  GPI_anchor <- target_domain %>%
    filter(type == "Lipidation", description == "GPI-anchor amidated serine" ) %>% 
    arrange(acc, start) %>%
    group_by(acc) %>%
    slice(1) %>%                               
    ungroup() %>%
    transmute(acc, type, dm1_start = start, dm1_end = end)
  
  #その他(分泌系)を抽出
  Signal<-target_domain %>%
    filter(type == "Signal") %>% 
    arrange(acc, start) %>%
    group_by(acc) %>%
    slice(1) %>%                               
    ungroup() %>%
    transmute(acc, type, dm1_start = start, dm1_end = end)%>%
    filter(!(acc %in% c(tm_anchor$acc, GPI_anchor$acc)))
  
  if(plot_class == "TM"){
    target_df<-tm_anchor
  }
  
  if(plot_class == "GPI"){
    target_df<-GPI_anchor
  }
  
  if(plot_class == "Secreted"){
    target_df<-Signal
  }
  
  if(plot_class == "All"){
    target_df<-ALL
  }
  
  
  #plot用feature
  feature_map<-target_domain%>%
    filter(type %in% c("Domain", "Signal", "Lipidation", "Binding site", "Active site", "Disulfide bond"))%>%
    left_join(select(target_df, acc, dm1_start, dm1_end), by = "acc") %>%
    filter(!is.na(dm1_start))
  
  
  #糖鎖修飾位置の抽出
  extract_NXS_TSC <- function(x) {
    parts <- strsplit(x, "/", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    keep <- grepl("\\(N[^P][STC][^)]*\\)", parts)
    out <- parts[keep]
    if (!any(keep)) return(NA_character_)
    
    parts[which(keep)[1]]
  }
  
  target_gcategory<-select(gpep, acc, GN, Category, site, ratio_cluster, site_label)%>% 
    mutate(first_site = vapply(site, extract_NXS_TSC, character(1)),
           first_site = ifelse(first_site == "NA_character_", site, first_site),
           first_position = as.integer(str_extract(first_site, "\\d+(?=\\()")),
           x_site = first_position)%>%
    distinct()%>%
    left_join(protein_length, by =c("GN", "acc"))
  
  CFB_cluster_site<-target_gcategory %>%
    filter(GN %in% unique(feature_map$GN), !is.na(ratio_cluster), Category == "CF_Bisect")%>%
    mutate(ratio_cluster = ifelse(ratio_cluster == 3, 2, ratio_cluster))%>%
    distinct()
  
  detected_site<-target_gcategory %>%
    filter(GN %in% unique(feature_map$GN))%>%
    select(acc, GN, site, site_label, first_position, first_site, x_site)%>%
    distinct()
  
  
  feature_tm <- target_domain %>%
    inner_join(target_df %>% select(acc, dm1_start, dm1_end), by = "acc") %>%
    mutate(
      domain_label = ifelse(is.na(description) | description=="", type, description),
      domain_family = str_replace(domain_label, "\\s+\\d+$", ""),
      domain_family = str_trim(domain_family),
      x_start = start,
      x_end   = end) %>%
    select(acc, GN, type, domain_label, domain_family, x_start, x_end)
  
  tm_df <- target_domain %>%
    filter(type == "Transmembrane") %>%
    inner_join(target_df %>% select(acc, dm1_start), by = "acc") %>%  # TM1 start基準を付与
    mutate( x_tm_start = start,
            x_tm_end   = end)
  
  domain_df <- feature_tm %>%
    filter(type %in% c("Domain", "Signal", "Binding site", "Active site")) %>%
    filter(domain_family %in% c("Ig-like C2-type", "Fibronectin type-III"))%>%
    transmute(acc, GN, domain_family, x_start, x_end)%>%
    filter(str_detect(domain_family,"(?i)egf|ig|fibronectin|binding site|signal"))
  
  #ploteinの全長
  prot_df <- target_domain %>%
    distinct(acc, GN, LEN) %>%
    inner_join(target_df %>% select(acc, dm1_start, dm1_end), by = "acc") %>%
    mutate(x_min = 1,          # N末端（TM1 startより左/右はscale_x_reverseで反転）
           x_max = LEN)
  
  
  # 並び順（GN順）
  
  row_step <- 2.0 #GNの間隔
  stack_step <- 0.5  # dotの間隔
  
  prot_df  <- prot_df  %>% mutate(GN = factor(GN, levels = gn_order))
  domain_df<- domain_df %>%
    left_join(distinct(prot_df, acc, GN), by=c("acc", "GN")) %>%
    mutate(GN = factor(GN, levels = gn_order))
  
  
  # GN順は既に gn_order を作っている前提
  detected_site <- detected_site %>%
    mutate(GN = factor(GN, levels = gn_order),
           y_base = as.integer(GN) * row_step,
           y_site = y_base + 0.35)
  
  CFB_cluster_site <- CFB_cluster_site %>%
    mutate(
      GN = factor(GN, levels = gn_order),
      ratio_cluster = factor(ratio_cluster, levels = c(1, 2)),
      y_base = as.integer(GN) * row_step,
      y_site = y_base + 0.7 + (as.integer(ratio_cluster) - 1) * stack_step
    )
  
  
  p <- ggplot() +
    # protein full length
    geom_segment(
      data = prot_df %>% mutate(GN = factor(GN, levels = gn_order),
                                y_base = as.integer(GN) * row_step),
      aes(x = x_min, xend = x_max, y = y_base, yend = y_base),
      linewidth = 2.2,
      color = "grey85"
    ) +
    # all TM segments
    geom_segment(
      data = tm_df %>% mutate(GN = factor(GN, levels = gn_order),
                              y_base = as.integer(GN) * row_step),
      aes(x = x_tm_start, xend = x_tm_end, y = y_base, yend = y_base),
      linewidth = 2,
      color = "black",
      lineend = "butt"
    ) +
    # domains (below)
    geom_rect(
      data = domain_df %>% mutate(GN = factor(GN, levels = gn_order),
                                  y_base = as.integer(GN) * row_step),
      aes(xmin = x_start, xmax = x_end,
          ymin = y_base - 0.45, ymax = y_base - 0.15,
          fill = domain_family),
      inherit.aes = FALSE,
      alpha = 0.7
    ) +
    # --- 背景：すべての糖ペプチドsite ---
    geom_point(
      data = detected_site,
      aes(x = x_site, y = y_site, shape = site_label),
      color = "#778899",
      size = 1.0,
      alpha = 0.6,
      na.rm = TRUE
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 17, `TRUE` = 8),
      name = "442_452 glycan change site"
    )+
    # --- 上書き：CFB (CF_Bisect) の動き ---
    geom_point(
      data = CFB_cluster_site,
      aes(x = x_site, y = y_site, color = ratio_cluster),
      size = 1.6,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = c(`1` = "#4169e1", `2` = "#ff6347"),
      name = "CFB ratio_cluster"
    ) +
    scale_y_continuous(
      breaks = seq_along(gn_order) * row_step,
      labels = gn_order,
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    # scale_x_reverse() +
    labs(x = "amino acid position", y = NULL) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),  
      plot.background  = element_rect(fill = "white", color = NA),  
      panel.grid = element_blank(), 
      axis.text.y  = element_text(color = "black", size =10, face = "bold"),
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank()) 
  
  p
  
}


domain_position<-function(gpep, target_domain){
  
  target_df <- target_domain %>%
    arrange(acc, start) %>%
    transmute(acc, type,description, dm1_start = start, dm1_end = end)%>%
    mutate(domain_clean = str_remove(description, " [0-9]+$"))
  
  #plot用feature
  feature_map<-gpep%>%left_join(select(target_df, acc, dm1_start, dm1_end), by = "acc") %>% filter(!is.na(dm1_start))
  
  
  #糖鎖修飾位置の抽出
  extract_NXS_TSC <- function(x) {
    parts <- strsplit(x, "/", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    keep <- grepl("\\(N[^P][STC][^)]*\\)", parts)
    out <- parts[keep]
    if (!any(keep)) return(NA_character_)
    
    parts[which(keep)[1]]
  }
  
  select_cluster<-c("ratio_cluster | igot_cluster")
  target_gcategory<-select(gpep, acc, GN, site, any_of(select_cluster))%>% 
    mutate(first_site = vapply(site, extract_NXS_TSC, character(1)),
           first_site = ifelse(first_site == "NA_character_", site, first_site),
           first_position = as.integer(str_extract(first_site, "\\d+(?=\\()")),
           x_site = first_position)%>%
    distinct()%>%
    left_join(protein_length, by =c("GN", "acc"))
  
  domain_info <- target_gcategory %>%
    left_join(target_df, by = "acc") %>%
    mutate(in_domain = x_site >= dm1_start & x_site <= dm1_end)
  
  domain_info
}

fisher_analaysis<-function(all_site_domain_df, tar_domain, tar_site){
  
  igot_site_domain <- all_site_domain_df %>%
    filter(grepl(tar_domain, domain_clean),in_domain)%>%
    distinct(acc, x_site, GN, site, .keep_all = TRUE)
  
  domain_match_tar<-list( "all_tar_site" = unique(paste0(igot_site_domain$acc, "_",igot_site_domain$site)), 
                          "tar_site" =    tar_site)
  domain_match_venn<-venn_analysis(domain_match_tar)
  match_site<-extract_items(domain_match_venn@label_df, "all_tar_site&tar_site")%>%unname()  
  
  all_sites <- all_site_domain_df %>%distinct(acc, x_site)
  domain_sites <- igot_site_domain %>%distinct(acc, x_site)
  
  a <- length(match_site)
  b <- length(tar_site)-a
  c <- nrow(domain_sites) - a
  d <- (nrow(all_sites) - nrow(domain_sites)) - b
  
  mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  res<- fisher.test(mat)
  
  # res<-tibble(match_site_num = a,
  #             other_target_glycan_change_site_num = b,
  #             all_target_domain_site_num = c,
  #             other_domain_site_num =d)%>%
  #      mutate()
  
  res
}


feature_position<-function(df, tar_site_df,  feat_type = "Transmembrane", desc = "Helical"){

    tm_anchor_tmp <- df %>%
      filter(type == feat_type, description == desc) %>%  
      transmute(acc, tm_start = as.integer(start), tm_end = as.integer(end)) %>%
      filter(!is.na(tm_start), !is.na(tm_end)) %>%
      mutate(tm_s = pmin(tm_start, tm_end),
             tm_e = pmax(tm_start, tm_end)) %>%
      select(acc, tm_s, tm_e)
    
  site_tm_min <- tar_site_df %>%
    select(acc, GN, site_pos, is_hit) %>%
    inner_join(tm_anchor_tmp, by="acc") %>%
    mutate(dist_to_seg = case_when(site_pos < tm_s ~ tm_s - site_pos,
                                   site_pos > tm_e ~ site_pos - tm_e,
                                   TRUE ~ 0L)) %>%
    group_by(acc, GN, site_pos, is_hit) %>%
    summarise(dist_TM_min = min(dist_to_seg, na.rm = TRUE),.groups="drop")
  
  site_tm_min 
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

GO_plot<-function(ego_df, GO_mode){
  
  plot_df <- ego_df %>%
    arrange(p.adjust) %>%
    head(10) %>%
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
