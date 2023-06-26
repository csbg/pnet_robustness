predictions %>% distinct(experiment)



plot_roc("mskimpact")
ggsave_publication("X_msk_ROC", width = 4, height = 4)

plot_robustness(experiment = "mskimpact")
ggsave_publication("X_msk_robustness", width = 18, height = 4)

plot_bias("mskimpact")
ggsave_publication("X_msk_scores", width = 14, height = 4)

plot_bias_comparison("mskimpact")
ggsave_publication("X_msk_comparison", width = 14, height = 4)




# ROC curves MSK-IMPACT ---------------------------------------------------

plot_roc("mskimpact")
ggsave_publication("roc_msk_nsclc", width = 4, height = 4, type = "png")

plot_roc("mskimpact_nsclc")
ggsave_publication("roc_msk_nsclc", width = 4, height = 4, type = "png")
plot_roc("mskimpact_nsclc_deterministic")
plot_roc("mskimpact_nsclc_shuffled")
ggsave_publication("roc_msk_nsclc_shuffled", width = 4, height = 4, type = "png")


plot_roc("mskimpact_bc")
ggsave_publication("roc_msk_bc", width = 4, height = 4, type = "png")
plot_roc("mskimpact_bc_deterministic")
plot_roc("mskimpact_bc_shuffled")
ggsave_publication("roc_msk_bc_shuffled", width = 4, height = 4, type = "png")

plot_roc("mskimpact_cc")
ggsave_publication("roc_msk_cc", width = 4, height = 4, type = "png")
plot_roc("mskimpact_cc_deterministic")
plot_roc("mskimpact_cc_shuffled")
ggsave_publication("roc_msk_cc_shuffled", width = 4, height = 4, type = "png")

plot_roc("mskimpact_pc")
ggsave_publication("roc_msk_pc", width = 4, height = 4, type = "png")
plot_roc("mskimpact_pc_deterministic")
ggsave_publication("roc_msk_pc_det", width = 4, height = 4, type = "png")
plot_roc("mskimpact_pc_shuffled")
ggsave_publication("roc_msk_pc_shuffled", width = 4, height = 4, type = "png")




# MSK-IMPACT Prostate cancer ----------------------------------------------

plot_robustness(experiment = "mskimpact_pc")
ggsave_publication("X_msk_pc_robustness", width = 18, height = 4, type = "png")

plot_bias("mskimpact_pc")
ggsave_publication("X_msk_pc_scores", width = 14, height = 4, type = "png")

plot_bias_comparison("mskimpact_pc")
ggsave_publication("X_msk_pc_comparison", width = 14, height = 4, type = "png")


# Shuffled alternative ----------------------------------------------------

plot_bias("pnet_shuffled_each")
ggsave_publication("X_3b_each_scores", width = 14, height = 4, type = "png")

plot_roc("pnet_shuffled_each")
ggsave_publication("X_3c_each_roc", width = 4, height = 4, type = "png")

plot_bias_comparison("pnet_shuffled_each")
ggsave_publication("X_3d_each_comparison", width = 14, height = 4, type = "png")



# Differential node weights -----------------------------------------------

plot_importance_difference <- function(comparison, reference) {
  node_importance %>%
    filter(
      experiment == .env$comparison | experiment == .env$reference
    ) %>%
    pivot_wider(names_from = experiment, values_from = coef_combined) %>%
    mutate(delta_node_importance = .data[[comparison]] - .data[[reference]]) %>%
    summarise(
      .by = c(layer, reactome_id),
      delta_node_importance = median(delta_node_importance),
    ) %>%
    mutate(.by = layer, x = row_number(delta_node_importance)) %>%
    ggplot(aes(x, delta_node_importance)) +
    geom_hline(yintercept = 0) +
    geom_point(
      aes(color = delta_node_importance > 0),
      show.legend = FALSE,
      alpha = 1,
      size = .2
    ) +
    xlab("nodes ranked by median importance difference") +
    ylab("median node importance difference") +
    facet_wrap(vars(layer), nrow = 1, scales = "free") +
    ggtitle(str_glue("{comparison} vs {reference}")) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  # {.}
}

plot_importance_difference("pnet_shuffled", "pnet_original")
ggsave_publication("X_delta_shuffled_vs_original", width = 18, height = 4, type = "png")

plot_importance_difference("pnet_shuffled_each", "pnet_original")
ggsave_publication("X_delta_shuffled_each_vs_original", width = 18, height = 4, type = "png")

plot_importance_difference("pnet_deterministic", "pnet_original")
ggsave_publication("X_delta_determ_vs_original", width = 18, height = 4, type = "png")



plot_importance_difference("mskimpact_pc_shuffled", "mskimpact_pc")
ggsave_publication("X_delta_msk_shuffled_vs_original", width = 18, height = 4, type = "png")

plot_importance_difference("mskimpact_pc_deterministic", "mskimpact_pc")
ggsave_publication("X_delta_msk_determ_vs_original", width = 18, height = 4, type = "png")



# Figure 1D ---------------------------------------------------------------

get_mean_difference <- function(seed) {
  node_importance %>%
    filter(
      experiment == "pnet_original",
      seed %in% c("234_20080808", {{seed}})
    ) %>%
    pivot_wider(names_from = seed, values_from = coef_combined) %>%
    rename(original_seed = `234_20080808`, other_seed = .data[[seed]]) %>%
    mutate(delta = other_seed - original_seed) %>%
    summarise(mean_delta = sum(delta)) %>%
    pull(mean_delta) %>%
    {.}
}

# get_mean_difference("3_3")

df <-
  tibble(seed = unique(node_importance$seed) %>% setdiff(original_seed)) %>%
  rowwise() %>%
  mutate(mean_delta = get_mean_difference(seed))

df %>%
  write_csv("plots/mean_delta.csv")

read_csv("plots/mean_delta.csv")



# MSK-IMPACT Correlaction -------------------------------------------------

plot_msk_heatmap <- function(layer = NULL,
                             show_legends = TRUE,
                             heatmap_size = 50) {
  if (!is.null(layer)) {
    node_importance <-
      node_importance %>%
      filter(layer == {{layer}})
  }

  corr_mat <-
    node_importance %>%
    filter(
      str_starts(experiment, "mskimpact_"),
      !str_ends(experiment, "once"),
      # seed %in% c(
      #   "0_0", "1_1", "2_2", "3_3", "4_4"
      #   # "5_5", "6_6", "7_7", "8_8", "9_9"
      # )
    ) %>%
    unite(experiment, seed, col = "exp_seed", sep = "+") %>%
    select(exp_seed, reactome_id, coef_combined) %>%
    pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
    select(!reactome_id) %>%
    cor(use = "pairwise.complete.obs")

  col_metadata <-
    tibble(exp_seed = colnames(corr_mat)) %>%
    separate(exp_seed, into = c("experiment", "seed"), sep = "\\+") %>%
    mutate(type = case_when(
      str_ends(experiment, "deterministic") ~ "deterministic",
      str_ends(experiment, "shuffled") ~ "shuffled",
      .default = "original"
    ))

  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(min(corr_mat), max(corr_mat), length.out = 9),
      color("davos", reverse = TRUE)(9),
    ),

    show_heatmap_legend = show_legends,
    name = "correlation of\nnode importance\nscores",
    heatmap_legend_param = list(
      at = round(c(min(corr_mat), max(corr_mat)), 2),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),

    clustering_distance_rows = function(m) as.dist(1 - m[, rownames(m)]),
    clustering_distance_columns = function(m) as.dist(1 - m[, colnames(m)]),
    cluster_row_slices = FALSE,
    row_dend_gp = gpar(lwd = 0.5),
    row_title = "seeds",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

    row_split = col_metadata$type,

    width = unit(heatmap_size, "mm"),
    height = unit(heatmap_size, "mm"),
    border = FALSE,

    show_column_dend = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,

    left_annotation = rowAnnotation(
      experiment = col_metadata$experiment,
      type = col_metadata$type,
      # col = list(experiment = EXPERIMENT_COLORS),
      show_annotation_name = FALSE,
      show_legend = show_legends,
      annotation_legend_param = list(
        experiment = list(
          title = "experiment",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  )
}

(p <- plot_msk_heatmap(heatmap_size = 100))
ggsave_publication("X_msk_heatmap_split_nocluster", plot = p,
                   width = 15, height = 15, type = "png")

df <-
  node_importance %>%
  filter(
    str_starts(experiment, "mskimpact_"),
    !str_ends(experiment, "once"),
    # seed %in% c(
    #   "0_0", "1_1", "2_2", "3_3", "4_4"
    #   # "5_5", "6_6", "7_7", "8_8", "9_9"
    # )
  ) %>%
  # unite(experiment, seed, col = "exp_seed", sep = "+") %>%
  # select(exp_seed, reactome_id, coef_combined) %>%
  # pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
  # select(!reactome_id) %>%
  # cor(use = "pairwise.complete.obs") %>%
  {.}
df


# Corrected node weights --------------------------------------------------

library(limma)

imp_df <-
  node_importance %>%
  filter(experiment %in% c("pnet_original", "pnet_shuffled"), layer == 1) %>%
  pivot_wider(names_from = reactome_id, values_from = coef_combined)
imp_df

imp_mat <-
  imp_df %>%
  unite(experiment, seed, layer, col = "esl", sep = "|") %>%
  column_to_rownames("esl") %>%
  as.matrix()

imp_mat %>% dim()

imp_mat_scaled <-
  imp_mat %>%
  scale() %>%
  t()

imp_mat_scaled %>% dim()


model_mat <- model.matrix(
  ~experiment,
  data =
    imp_df %>%
    select(experiment, seed) %>%
    mutate(experiment = fct_relevel(experiment, "pnet_shuffled"))
)

efit <-
  imp_mat_scaled %>%
  lmFit(model_mat) %>%
  eBayes()

summary(decideTests(efit))

topTable(efit, number = Inf) %>%
  as_tibble(rownames = "node") %>%
  filter(adj.P.Val <= 0.05)





# DTox --------------------------------------------------------------------

dtox_relevance <-
  dir_ls(glob = "data/dtox/*/module_relevance.tsv", recurse = TRUE) %>%
  map(
    \(file) {
      df <- read_tsv(file)
      colnames(df)[1] <- "compound"
      df
    }
  ) %>%
  list_rbind(names_to = "file") %>%
  separate_wider_regex(
    file,
    c("data/dtox/", seed = ".+", "/.+"),,
  ) %>%
  {.}


plot_dtox_relevance <- function(compound = "CID_75146", by_median = TRUE) {
  data_long <-
    dtox_relevance %>%
    filter(compound == {{compound}}) %>%
    pivot_longer(
      !c(seed, compound),
      names_to = "node",
      values_to = "importance"
    )

  if (by_median) {
    top_nodes <-
      data_long %>%
      summarise(.by = node, importance = median(importance))
  } else {
    top_nodes <-
      data_long %>%
      summarise(.by = node, importance = max(importance) - min(importance))
  }

  top_nodes <-
    top_nodes %>%
    slice_max(importance, n = 30) %>%
    pull(node)

  data_long %>%
    filter(node %in% top_nodes) %>%
    mutate(node = fct_reorder(node, importance)) %>%
    ggplot(aes(node, importance)) +
    geom_boxplot() +
    coord_flip() +
    ggtitle(compound)
}

plot_dtox_relevance() + ylim(-2, 2)
ggsave_publication("X_dtox_relevance_example", width = 15, height = 12)
plot_dtox_relevance(by_median = FALSE)
ggsave_publication("X_dtox_relevance_example2", width = 15, height = 12)


dtox_relevance %>%
  filter(compound == "CID_75146") %>%
  select(seed, "P27986") %>%
  View() %>%
  {.}




# Dtox ROC ----------------------------------------------------------------

predictions_dtox <-
  dir_ls(glob = "data/dtox/*/test_labels.csv", recurse = TRUE) %>%
  map(\(file) read_csv(file)) %>%
  list_rbind(names_to = "file") %>%
  separate_wider_regex(
    file,
    c("data/dtox/", seed = ".+", "/.+")
  ) %>%
  mutate(
    truth =
      if_else(truth > 0.9, "one", "zero") %>%
      as_factor() %>%
      fct_relevel("one", "zero")
  )

predictions_dtox %>%
  distinct(truth)



plot_roc_dtox <- function() {
  pred <-
    predictions_dtox %>%
    group_by(seed)

  plot_data <-
    group_keys(pred) %>%
    mutate(
      metrics =
        pred %>%
        group_split() %>%
        map(roc_curve, truth, predicted)
    ) %>%
    unnest(metrics)

  ggplot(plot_data, aes(1 - specificity, sensitivity, group = seed)) +
    geom_path(alpha = .1, linewidth = BASE_LINEWIDTH) +
    geom_path(
      data = plot_data %>% filter(seed == "0"),
      linewidth = BASE_LINEWIDTH,
      color = ORIGINAL_SEED_COLOR
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_fixed() +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_roc_dtox()
ggsave_publication("X_dtox_roc", width = 4, height = 4)
