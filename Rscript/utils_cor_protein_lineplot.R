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
    theme_bw(base_size = 12)
  
  p}

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

trans_longdata_gpep_lineplot<-function(plot_df){
  
  time_point <- c("0_1", "0_2", "0_3",
                  "8_1", "8_2", "8_3",
                  "14_1", "14_2", "14_3",
                  "28_1", "28_2", "28_3",
                  "84_1", "84_2", "84_3")
  
  names(plot_df)[-(1:5)] <- time_point
  
  long_data <- plot_df %>%
    pivot_longer(
      cols = -c("GN", "site_ID", "glycan_ID","acc_site_glycan_ID", "gpep_cluster"),
      names_to = c("Time", "replicate"),
      names_sep = "_",
      values_to = "Score")%>%
    mutate(group_key = ifelse( gpep_cluster == 1, "decre", "incre"))%>%
    mutate(GN_group_key = paste0(GN, "_" ,group_key))%>%
    group_by(GN, group_key, Time)%>%
    mutate(GN_mean_expression_pattern = mean(Score, na.rm = TRUE),
           Time = as.numeric(Time))%>%ungroup()
  
  long_data
  
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

expression_plot <- function(df) {
  
  plot_df <- df %>%
    select(GN, site_ID, group_key, pep_group_key,Time, GN_mean_expression_pattern) %>%
    distinct()
  
  pep_info <- plot_df %>% distinct(GN,pep_group_key, group_key)
  
  GN_decre <- pep_info %>% filter(group_key == "decre") %>% pull(GN)
  GN_incre <- pep_info %>% filter(group_key == "incre") %>% pull(GN)
  
  ## グラデーション色定義
  cols_decre <- colorRampPalette(
    c("#08306B", "#4292C6", "#DEEBF7")
  )(length(GN_decre))
  
  cols_incre <- colorRampPalette(
    c("#dc143c", "#ff8c00", "#ffd700")
  )(length(GN_incre))
  
  names(cols_decre) <- GN_decre
  names(cols_incre) <- GN_incre
  
  cols <- c(cols_decre, cols_incre)
  
  p <- ggplot(
    plot_df,
    aes(
      x = Time,
      y = GN_mean_expression_pattern,
      color = GN,
      group = pep_group_key
    )
  ) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = cols) +
    geom_vline(
      xintercept = c(0, 8, 14, 28, 84),
      color = "dimgray",
      linewidth = 0.5,
      linetype = "dotted"
    ) +
    scale_x_continuous(
      breaks = c(0, 8, 14, 28, 84),
      labels = c("P0", "P8", "2w", "1M", "3M")
    ) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 20, face = "bold", color = "black"),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.length.y = unit(-3, "pt"),
      legend.text = element_text(size = 15)
    )
  
  p
}
