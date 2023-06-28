# Corheatmap --------------------------------------------------------------

normalize_importances <- function(original, control) {
  normalized_mat <-
    node_importance %>%
    filter(experiment %in% c(original, control)) %>%
    unite(experiment, seed, col = "exp_seed", sep = "/") %>%
    select(!layer) %>%
    pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
    column_to_rownames("reactome_id") %>%
    as.matrix() %>%
    normalizeQuantiles()

  idx_original <- str_detect(colnames(normalized_mat), original)

  tibble(
    experiment = original,
    reactome_id = rownames(normalized_mat),
    biased = rowMeans(normalized_mat[, idx_original]),
    control = rowMeans(normalized_mat[, !idx_original]),
    corrected = biased - control
  )
}

corr_mat <-
  tribble(
    ~original, ~control,
    "mskimpact_nsclc_original", "mskimpact_nsclc_shuffled",
    "mskimpact_bc_original", "mskimpact_bc_shuffled",
    "mskimpact_cc_original", "mskimpact_cc_shuffled",
    "mskimpact_pc_original", "mskimpact_pc_shuffled"
  ) %>%
  pmap(normalize_importances) %>%
  list_rbind() %>%
  select(!control) %>%
  pivot_wider(
    names_from = experiment,
    values_from = c(biased, corrected),
    names_glue = "{experiment}_{.value}"
  ) %>%
  select(!reactome_id) %>%
  cor(use = "pairwise.complete.obs")

color_limits <- range(corr_mat)

plot_heatmap <- function(corr_mat,
                         color_limits,
                         show_legends = TRUE,
                         heatmap_size = 40) {
  cor_range <- range(corr_mat[lower.tri(corr_mat)]) %>% round(2)
  message("Correlation range: ",
          cor_range[1], " to ", cor_range[2])

  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(color_limits[1], color_limits[2], length.out = 9),
      color("davos", reverse = TRUE)(9),
    ),

    show_heatmap_legend = show_legends,
    name = "correlation of\nnode importance\nscores",
    heatmap_legend_param = list(
      at = round(color_limits, 2),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),

    clustering_distance_rows = function(m) as.dist(1 - m[, rownames(m)]),
    clustering_distance_columns = function(m) as.dist(1 - m[, colnames(m)]),
    row_dend_gp = gpar(lwd = 0.5),

    width = unit(heatmap_size, "mm"),
    height = unit(heatmap_size, "mm"),
    border = FALSE,

    show_column_dend = FALSE,
    show_column_names = FALSE,

    show_row_names = TRUE,
    # row_labels = EXPERIMENT_NAMES[rownames(corr_mat)],
    row_title = "cancer"
  )
}

(p <- plot_heatmap(corr_mat, color_limits))
ggsave_publication("S6b_msk_correlation_alternative",
                   plot = p, height = 6, width = 12)
