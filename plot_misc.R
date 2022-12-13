# More ideas --------------------------------------------------------------

## Importance vs measure ----

normalize <- function(x) {
  if (max(x) == min(x))
    rep.int(0.5, length(x))
  else
    (x - min(x)) / (max(x) - min(x))
}

plot_importance_vs_measure <- function(measure, method = "pearson") {
  plot_data <-
    node_importance %>%
    filter(!modified) %>%
    group_by(experiment, layer, reactome_id)

  if (method == "pearson") {
    plot_data <-
      plot_data %>%
      summarise(coef = mean(coef)) %>%
      left_join(graph_stats, by = c("reactome_id", "layer")) %>%
      mutate(
        measure = normalize({{measure}}),
        coef = normalize(coef)
      )
    y_label <- "relative node importance (mean)"
  } else if (method == "pearson_scaled") {
    plot_data <-
      plot_data %>%
      summarise(coef = mean(coef_scaled)) %>%
      left_join(graph_stats, by = c("reactome_id", "layer")) %>%
      mutate(
        measure = normalize({{measure}}),
        coef = normalize(coef)
      )
    y_label <- "relative node importance (mean scaled)"
  } else {
    plot_data <-
      plot_data %>%
      summarise(coef = median(coef_rank)) %>%
      left_join(graph_stats, by = c("reactome_id", "layer")) %>%
      mutate(
        measure = rank({{measure}}),
        measure = normalize({{measure}}),
        coef = normalize(coef)
      )
    y_label <- "relative node importance (median rank)"
  }

  x_label <- str_glue("{rlang::as_name(rlang::enquo(measure))} (relative)")

  ggplot(plot_data, aes(measure, coef)) +
    geom_point(size = 1, alpha = .25) +
    geom_smooth(method = "lm", size = BASE_LINE_SIZE) +
    scale_x_continuous(x_label, limits = 0:1, breaks = 0:1) +
    scale_y_continuous(y_label, limits = 0:1, breaks = 0:1) +
    facet_grid(vars(experiment), vars(layer)) +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_importance_vs_measure(reachability, "pearson")
ggsave_publication("z_importance_vs_reachability_pearson",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(reachability, "pearson_scaled")
ggsave_publication("z_importance_vs_reachability_pearson_scaled",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(reachability, "spearman")
ggsave_publication("z_importance_vs_reachability_spearman",
                   width = 18, height = 9, type = "png")


plot_importance_vs_measure(betweenness, "pearson")
ggsave_publication("z_importance_vs_betweenness_pearson",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(betweenness, "pearson_scaled")
ggsave_publication("z_importance_vs_betweenness_pearson_scaled",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(betweenness, "spearman")
ggsave_publication("z_importance_vs_betweenness_spearman",
                   width = 18, height = 9, type = "png")


plot_importance_vs_measure(degree, "pearson")
ggsave_publication("z_importance_vs_degree_pearson",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(degree, "pearson_scaled")
ggsave_publication("z_importance_vs_degree_pearson_scaled",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(degree, "spearman")
ggsave_publication("z_importance_vs_degree_spearman",
                   width = 18, height = 9, type = "png")


plot_importance_vs_measure(indegree, "pearson")
ggsave_publication("z_importance_vs_indegree_pearson",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(indegree, "pearson_scaled")
ggsave_publication("z_importance_vs_indegree_pearson_scaled",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(indegree, "spearman")
ggsave_publication("z_importance_vs_indegree_spearman",
                   width = 18, height = 9, type = "png")


plot_importance_vs_measure(outdegree, "pearson")
ggsave_publication("z_importance_vs_outdegree_pearson",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(outdegree, "pearson_scaled")
ggsave_publication("z_importance_vs_outdegree_pearson_scaled",
                   width = 18, height = 9, type = "png")

plot_importance_vs_measure(outdegree, "spearman")
ggsave_publication("z_importance_vs_outdegree_spearman",
                   width = 18, height = 9, type = "png")



## Correlations boxplot ----

correlate_importances <- function(method = "pearson", scaled = FALSE) {
  if (scaled)
    node_importance$coef <- node_importance$coef_scaled

  node_importance %>%
    filter(!modified) %>%
    left_join(graph_stats, by = c("reactome_id", "layer")) %>%
    group_split(experiment, seed, layer) %>%
    map_dfr(
      ~bind_cols(
        distinct(.x, experiment, seed, layer),
        cor(.x$coef, select(.x, indegree:betweenness), method = method) %>%
          as_tibble()
      )
    )
}

cor_data_pearson <- correlate_importances("pearson")
cor_data_pearson_scaled <- correlate_importances("pearson", TRUE)
cor_data_spearman <- correlate_importances("spearman")


plot_correlation_boxplots <- function(cor_data) {
  cor_data %>%
    pivot_longer(
      indegree:betweenness,
      names_to = "measure",
      values_to = "correlation"
    ) %>%
    mutate(
      layer = factor(layer),
      measure = fct_inorder(measure)
    ) %>%
    ggplot(aes(experiment, correlation, fill = experiment)) +
    geom_boxplot(
      outlier.size = .5,
      outlier.alpha = .25,
      key_glyph = "rect",
      size = BASE_LINE_SIZE,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = EXPERIMENT_COLORS,
      guide = guide_legend()
    ) +
    facet_grid(vars(layer), vars(measure)) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      axis.text.x = element_text(angle = 90, vjust = .5)
    )
}

plot_correlation_boxplots(cor_data_pearson)
ggsave_publication("z_correlations_pearson", type = "png", width = 8, height = 12)

plot_correlation_boxplots(cor_data_pearson_scaled)
ggsave_publication("z_correlations_pearson_scaled", type = "png", width = 8, height = 12)

plot_correlation_boxplots(cor_data_spearman)
ggsave_publication("z_correlations_spearman", type = "png", width = 8, height = 12)



## Correlation heatmap ----

plot_correlation_heatmap <- function(cor_data) {
  cor_data %>%

    # group_by(experiment) %>%
    # summarise(across(indegree:betweenness, mean, na.rm = TRUE)) %>%

    # group_by(experiment, layer) %>%
    # summarise(across(indegree:betweenness, mean, na.rm = TRUE)) %>%
    # unite(experiment, experiment, layer) %>%

    group_by(experiment, layer) %>%
    summarise(across(indegree:betweenness, mean, na.rm = TRUE)) %>%
    pivot_longer(indegree:betweenness, names_to = "measure") %>%
    unite(experiment, experiment, measure) %>%
    pivot_wider(names_from = layer) %>%

    column_to_rownames("experiment") %>%
    as.matrix() %>%
    t() %>%
    Heatmap(
      col = circlize::colorRamp2(
        seq(min(., na.rm = TRUE), max(., na.rm = TRUE), length.out = 9),
        color("davos", reverse = TRUE)(9),
      ),
      name = "correlation of\nnode importance\nscores vs\nnetwork measures",
      heatmap_legend_param = list(
        at = round(c(min(., na.rm = TRUE), max(., na.rm = TRUE)), 2),
        border = FALSE,
        grid_width = unit(2, "mm"),
        labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        legend_height = unit(15, "mm"),
        title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
      ),

      row_dend_gp = gpar(lwd = 0.5),
      row_title = "network measure",
      row_title_side = "right",
      row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

      column_dend_gp = gpar(lwd = 0.5),
      column_title = "experiment",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

      width = unit(5, "mm") * ncol(.),
      height = unit(5, "mm") * nrow(.),
      border = FALSE,
    )
}

plot_correlation_heatmap <- function(cor_data) {
  cor_data %>%
    group_by(experiment, layer) %>%
    summarise(across(indegree:betweenness, mean, na.rm = TRUE)) %>%
    pivot_longer(indegree:betweenness, names_to = "measure") %>%
    unite(experiment, experiment, measure) %>%
    pivot_wider(names_from = layer) %>%
    column_to_rownames("experiment") %>%
    as.matrix() %>%
    t() %>%
    Heatmap(
      col = circlize::colorRamp2(
        seq(min(., na.rm = TRUE), max(., na.rm = TRUE), length.out = 9),
        color("davos", reverse = TRUE)(9),
      ),
      name = "correlation of\nnode importance\nscores vs\nnetwork measures",
      heatmap_legend_param = list(
        at = round(c(min(., na.rm = TRUE), max(., na.rm = TRUE)), 2),
        border = FALSE,
        grid_width = unit(2, "mm"),
        labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
        legend_height = unit(15, "mm"),
        title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
      ),

      cluster_rows = FALSE,
      cluster_columns = FALSE,

      row_dend_gp = gpar(lwd = 0.5),
      row_title = "layer",
      row_title_side = "right",
      row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

      column_dend_gp = gpar(lwd = 0.5),
      column_title = "experiment_measure",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

      width = unit(5, "mm") * ncol(.),
      height = unit(5, "mm") * nrow(.),
      border = FALSE,
    )
}

(p <- plot_correlation_heatmap(cor_data_pearson))
ggsave_publication("z_corr_heatmap_pearson3", plot = p,
                   width = 10, height = 6, type = "png")

(p <- plot_correlation_heatmap(cor_data_pearson_scaled))
ggsave_publication("z_corr_heatmap_pearson_scaled", plot = p,
                   width = 6, height = 6, type = "png")

(p <- plot_correlation_heatmap(cor_data_spearman))
ggsave_publication("z_corr_heatmap_spearman", plot = p,
                   width = 6, height = 6, type = "png")

