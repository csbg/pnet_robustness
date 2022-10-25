library(tidyverse)
library(fs)
library(yardstick)
library(ComplexHeatmap)
library(khroma)
library(igraph)
source("styling.R")

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(1, "mm")
)



# Load data ---------------------------------------------------------------

original_seed <- "234_20080808"

reactome_names <-
  read_tsv(
    "pnet_prostate_paper/_database/pathways/Reactome/ReactomePathways.txt",
    col_names = c("reactome_id", "node", "species")) %>%
  filter(species == "Homo sapiens") %>%
  select(!species)

node_importance <-
  dir_ls(
    glob = str_glue("data/*/*/node_importance.csv"),
    recurse = TRUE
  ) %>%
  map_dfr(
    read_csv,
    .id = "file"
  ) %>%
  rename(node = ...1) %>%
  left_join(reactome_names, by = "node") %>%
  mutate(reactome_id = if_else(layer == 1, node, reactome_id)) %>%
  extract(
    file,
    into = c("experiment", "seed"),
    regex = str_glue("data/(.+)/(.+)/"),
    convert = TRUE
  ) %>%
  group_by(experiment, seed, layer) %>%
  mutate(
    modified = coef_graph > mean(coef_graph) + 5 * sd(coef_graph),
    coef_scaled = scale(coef)[,1],
    coef_rank = rank(coef),
    experiment = factor(experiment, levels = names(EXPERIMENT_COLORS))
  ) %>%
  ungroup() %>%
  select(experiment, seed, layer, coef, coef_scaled, coef_rank,
         reactome_id, modified)

predictions <-
  dir_ls(
    glob = "data/*/*/predictions_test.csv",
    recurse = TRUE
  ) %>%
  map_dfr(
    ~read_csv(., col_types = "cfdf"),
    .id = "file"
  ) %>%
  extract(
    file,
    into = c("experiment", "seed"),
    regex = str_glue("data/(.+?)/(.+?)/"),
    convert = TRUE
  ) %>%
  transmute(
    experiment,
    seed,
    obs = ...1,
    truth = fct_recode(y, class_1 = "1", class_0 = "0"),
    class_1 = pred_scores,
    predicted = fct_recode(pred, class_1 = "1.0", class_0 = "0.0")
  )

# only include nodes with connectivity; therefore, load link weight matrices of
# NN and search for nonzero entries
pnet_edges <-
  map_dfr(
    0:6,
    function(layer) {
      df <- read_csv(str_glue("pnet_prostate_paper/analysis/extracted/link_weights_{layer}.csv"))
      if (layer == 0) {
        colnames(df)[1] <- "to"
        colnames(df)[2] <- "from"
      } else {
        colnames(df)[1] <- "from"
        df <-
          df %>%
          pivot_longer(!from, names_to = "to", values_to = "has_edge") %>%
          filter(has_edge != 0) %>%
          select(!has_edge)
      }
      df %>%
        transmute(
          from = paste0(from, ":", layer - 1),
          to = paste0(to, ":", layer)
        )
    }
  )

pnet_graph <-
  pnet_edges %>%
  graph_from_data_frame()

# PNET sets indegrees for all genes (layer 0) to 1 and degree = outdegree
# here, indegree is (correctly) 3
graph_stats <-
  tibble(
    reactome_id = names(V(pnet_graph)),
    indegree = degree(pnet_graph, mode = "in"),
    outdegree = degree(pnet_graph, mode = "out"),
    degree = indegree + outdegree,
    reachability = map_int(
      reactome_id,
      ~subcomponent(pnet_graph, .x, "in") %>% length()
    ),
    betweenness = betweenness(pnet_graph)
  ) %>%
  separate(reactome_id, into = c("reactome_id", "layer"), sep = "\\:") %>%
  mutate(layer = as.integer(layer) + 1L)



# Figure 1 ----------------------------------------------------------------

## c ----

plot_roc <- function(experiment) {
  pred <-
    predictions %>%
    filter(experiment == {{experiment}}) %>%
    group_by(seed)

  plot_data <-
    group_keys(pred) %>%
    mutate(
      metrics =
        pred %>%
        group_split() %>%
        map(roc_curve, truth, class_1)
    ) %>%
    unnest(metrics)

  ggplot(plot_data, aes(1 - specificity, sensitivity, group = seed)) +
    geom_path(alpha = .1, size = BASE_LINE_SIZE) +
    geom_path(
      data = plot_data %>% filter(seed == original_seed),
      size = BASE_LINE_SIZE,
      color = ORIGINAL_SEED_COLOR
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_fixed() +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_roc("default")
ggsave_publication("1c_roc", width = 4, height = 4)


## d ----

plot_robustness <- function(top_nodes = 10) {
  node_importance <-
    node_importance %>%
    filter(experiment == "default")

  top_nodes_per_layer <-
    node_importance %>%
    filter(seed == original_seed) %>%
    group_by(layer) %>%
    slice_max(coef_combined, n = top_nodes, with_ties = FALSE) %>%
    arrange(layer, desc(coef_combined))

  node_importance <-
    node_importance %>%
    semi_join(top_nodes_per_layer, by = c("layer", "reactome_id")) %>%
    mutate(reactome_id = factor(reactome_id,
                                levels = top_nodes_per_layer$reactome_id)) %>%
    group_by(layer) %>%
    mutate(coef_combined = scale(coef_combined)) %>%
    ungroup()

  ggplot(node_importance, aes(reactome_id, coef_combined)) +
    geom_boxplot(
      aes(color = experiment),
      outlier.shape = NA,
      lwd = BASE_BOXPLOT_SIZE,
      show.legend = FALSE
    ) +
    geom_point(
      data = node_importance %>% filter(seed == original_seed),
      color = ORIGINAL_SEED_COLOR,
      size = .75,
    ) +
    scale_x_discrete(paste("ten most important nodes per layer",
                           "ordered by importance in the original network")) +
    scale_y_continuous("node importance (z-score)") +
    scale_color_manual(values = EXPERIMENT_COLORS) +
    facet_wrap(vars(layer), scales = "free_x", nrow = 1) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
}

plot_robustness()
ggsave_publication("1d_robustness", width = 18, height = 4)


## e ----

plot_changes <- function(seed = "1_1", top_nodes = c(4, 43)) {
  plot_data <-
    node_importance %>%
    filter(experiment == "default", layer == 1)

  top_nodes <-
    plot_data %>%
    filter(seed == original_seed) %>%
    slice_max(coef_combined, n = top_nodes[2], with_ties = FALSE) %>%
    slice_min(coef_combined, n = top_nodes[2] - top_nodes[1] + 1, with_ties = FALSE) %>%
    arrange(desc(coef_combined))

  plot_data <-
    plot_data %>%
    filter(seed %in% c(original_seed, {{seed}})) %>%
    semi_join(top_nodes, by = "reactome_id") %>%
    mutate(reactome_id = factor(reactome_id,
                                levels = top_nodes$reactome_id)) %>%
    select(seed, reactome_id, coef_combined)

  ggplot(plot_data, aes(reactome_id, coef_combined)) +
    geom_segment(
      data = plot_data %>%
        mutate(
          seed = fct_recode(seed,
                            original = original_seed,
                            replicate = {{seed}})
        ) %>%
        pivot_wider(names_from = seed, values_from = coef_combined),
      aes(xend = reactome_id, yend = replicate, y = original),
      arrow = arrow(length = unit(1, "mm"), ends = "last", type = "closed"),
      size = BASE_LINE_SIZE
    ) +
    geom_point(
      data =
        plot_data %>%
        filter(seed == original_seed),
      size = 1,
      color = ORIGINAL_SEED_COLOR,
      show.legend = FALSE
    ) +
    scale_x_discrete("node") +
    scale_y_continuous("node importance") +
    theme_pub(rotate_x_labels = TRUE) +
    theme(
      panel.grid = element_blank()
    )
}

plot_changes("28_28")
ggsave_publication("1e_changes", width = 8, height = 4)



# Figure 2 ----------------------------------------------------------------

## b ----

plot_bias <- function(experiment, top_nodes = 5) {
  node_importance <-
    node_importance %>%
    filter(experiment == {{experiment}})

  top_nodes_per_layer <-
    node_importance %>%
    group_by(layer, reactome_id) %>%
    summarise(coef_combined = median(coef_combined)) %>%
    {
      bind_rows(
        slice_max(., coef_combined, n = top_nodes, with_ties = FALSE),
        slice_min(., coef_combined, n = top_nodes, with_ties = FALSE)
      )
    } %>%
    arrange(layer, desc(coef_combined))

  node_importance <-
    node_importance %>%
    semi_join(top_nodes_per_layer, by = c("layer", "reactome_id")) %>%
    mutate(reactome_id = factor(reactome_id,
                                levels = top_nodes_per_layer$reactome_id)) %>%
    group_by(layer) %>%
    mutate(coef_combined = scale(coef_combined)) %>%
    ungroup()

  ggplot(node_importance, aes(reactome_id, coef_combined)) +
    geom_boxplot(
      aes(color = experiment),
      outlier.shape = NA,
      lwd = BASE_BOXPLOT_SIZE,
      show.legend = FALSE
    ) +
    scale_x_discrete(paste("five most and least important nodes per layer",
                           "ordered by median importance")) +
    scale_y_continuous("node importance (z-score)") +
    scale_color_manual(values = EXPERIMENT_COLORS) +
    facet_wrap(vars(layer), scales = "free_x", nrow = 1) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
}

plot_bias("correlated")
ggsave_publication("2b_correlated", width = 14, height = 4)


## c ----

plot_bias_comparison <- function(experiment, top_nodes = 5) {
  node_importance <-
    node_importance %>%
    filter(experiment %in% c("default", {{experiment}}))

  top_nodes_per_layer <-
    node_importance %>%
    group_by(experiment, layer, reactome_id) %>%
    summarise(coef_combined = median(coef_combined)) %>%
    {
      bind_rows(
        slice_max(., coef_combined, n = top_nodes, with_ties = FALSE),
        slice_min(., coef_combined, n = top_nodes, with_ties = FALSE)
      )
    } %>%
    arrange(experiment, layer, desc(coef_combined)) %>%
    mutate(x_order = row_number()) %>%
    select(!coef_combined)

  node_importance <-
    node_importance %>%
    inner_join(
      top_nodes_per_layer,
      by = c("experiment", "layer", "reactome_id")
    ) %>%
    mutate(reactome_id = factor(x_order)) %>%
    group_by(layer) %>%
    mutate(coef_combined = scale(coef_combined)) %>%
    ungroup()

  ggplot(node_importance, aes(reactome_id, coef_combined)) +
    geom_boxplot(
      aes(color = experiment),
      outlier.shape = NA,
      lwd = BASE_LINE_SIZE * 1.5,
      show.legend = FALSE
    ) +
    scale_x_discrete(paste("five most and least important nodes per layer",
                           "ordered by median importance")) +
    scale_y_continuous("node importance (z-score)") +
    scale_color_manual(values = EXPERIMENT_COLORS) +
    facet_wrap(vars(layer), scales = "free_x", nrow = 1) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
}

plot_bias_comparison("correlated")
ggsave_publication("2c_correlated_comparison", width = 14, height = 4)



## d ----

plot_roc("correlated")
ggsave_publication("2d_roc", width = 4, height = 4)



# Figure 3 ----------------------------------------------------------------

## b ----

plot_bias("scrambled")
ggsave_publication("3b_scrambled", width = 14, height = 4)


## c ----

plot_bias_comparison("scrambled")
ggsave_publication("3c_scrambled_comparison", width = 14, height = 4)


## d ----

plot_roc("scrambled")
ggsave_publication("3d_roc", width = 4, height = 4)



# Figure 4 ----------------------------------------------------------------

plot_corr_heatmap <- function(layer = NULL) {
  if (!is.null(layer)) {
    node_importance <-
      node_importance %>%
      filter(layer == {{layer}})
  }

  corr_mat <-
    node_importance %>%
    unite(experiment, seed, col = "exp_seed", sep = "+") %>%
    select(exp_seed, reactome_id, coef_combined) %>%
    pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
    select(!reactome_id) %>%
    cor()

  distance <- as.dist(1 - corr_mat)

  col_metadata <-
    tibble(exp_seed = colnames(corr_mat)) %>%
    separate(exp_seed, into = c("experiment", "seed"), sep = "\\+")
  col_metadata

  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(min(corr_mat), max(corr_mat), length.out = 9),
      color("davos", reverse = TRUE)(9),
    ),
    name = "correlation of\nnode importance\nscores",
    heatmap_legend_param = list(
      at = round(c(min(corr_mat), max(corr_mat)), 2),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),

    clustering_distance_rows = distance,
    clustering_distance_columns = distance,
    row_dend_gp = gpar(lwd = 0.5),
    row_title = "seeds",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

    width = unit(50, "mm"),
    height = unit(50, "mm"),
    border = FALSE,

    show_column_dend = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,

    left_annotation = rowAnnotation(
      experiment = col_metadata$experiment,
      col = list(experiment = EXPERIMENT_COLORS),
      show_annotation_name = FALSE,
      show_legend = TRUE,
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

(p <- plot_corr_heatmap())
ggsave_publication("4_corr_heatmap", plot = p,
                   width = 11, height = 6, type = "png")

walk(
  1:6,
  function(l) {
    p <- plot_corr_heatmap(l)
    ggsave_publication(str_glue("4_corr_heatmap_layer_{l}"), plot = p,
                       width = 11, height = 6, type = "png")
  }
)



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
ggsave_publication("z_correlations_pearson", width = 8, height = 12)

plot_correlation_boxplots(cor_data_pearson_scaled)
ggsave_publication("z_correlations_pearson_scaled", width = 8, height = 12)

plot_correlation_boxplots(cor_data_spearman)
ggsave_publication("z_correlations_spearman", width = 8, height = 12)



## Correlation heatmap ----

p <-
  cor_data %>%
  group_by(experiment) %>%
  summarise(across(indegree:betweenness, mean, na.rm = TRUE)) %>%
  column_to_rownames("experiment") %>%
  as.matrix() %>%
  t() %>%
  Heatmap(
    col = circlize::colorRamp2(
      seq(min(.), max(.), length.out = 9),
      color("davos", reverse = TRUE)(9),
    ),
    name = "correlation of\nnode importance\nscores vs\nnetwork measures",
    heatmap_legend_param = list(
      at = round(c(min(.), max(.)), 2),
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

    width = unit(15, "mm"),
    height = unit(25, "mm"),
    border = FALSE,
  )
ggsave_publication("z_corr_heatmap", plot = p,
                   width = 6, height = 6, type = "png")









# Unused ------------------------------------------------------------------

## Importance vs degree (hexbin) ----

plot_importance_vs_degree <- function(experiment, measure) {
  print(str_glue("{experiment}, {measure}"))
  plot_data <-
    node_importance %>%
    filter(experiment == {{experiment}}, !modified) %>%
    # filter(seed == original_seed) %>%
    left_join(
      graph_stats %>%
        separate(reactome_id, into = c("reactome_id", "layer"), sep = "\\:") %>%
        mutate(layer = as.integer(layer) + 1),
      by = c("reactome_id", "layer")
    )

  color_max <- 100
  # color_max <- quantile(plot_data %>% pull(.data[[measure]]), 0.999)

  ggplot(plot_data, aes(.data[[measure]], coef)) +
    # geom_point(alpha = .15, size = .25) +
    geom_hex(bins = 25) +
    geom_smooth(method = "lm", size = BASE_LINE_SIZE) +
    # scale_y_log10() +
    scale_fill_distiller(
      palette = "YlGnBu",
      direction = 1,
      limits = c(0, color_max),
      breaks = c(0, color_max),
      labels = c("0", str_glue("{color_max} and above")),
      oob = scales::oob_squish,
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      ),
    ) +
    xlab(measure) +
    ylab("node importance") +
    facet_wrap(vars(layer), nrow = 1, scales = "free") +
    theme_pub()
}


plot_params <- tribble(
  ~experiment, ~measure,
  "default",    "degree",
  "correlated", "degree",
  "scrambled",  "degree",

  "default",    "indegree",
  "correlated", "indegree",
  "scrambled",  "indegree",

  "default",    "outdegree",
  "correlated", "outdegree",
  "scrambled",  "outdegree",

  "default",    "reachability",
  "correlated", "reachability",
  "scrambled",  "reachability",

  "default",    "betweenness",
  "correlated", "betweenness",
  "scrambled",  "betweenness"
)

pwalk(
  plot_params,
  function(experiment, measure) {
    plot_importance_vs_degree(experiment, measure)
    ggsave_publication(
      str_glue("s1_importance_vs_{measure}_{experiment}"),
      width = 18,
      height = 3,
      type = "png"
    )
  }
)
