library(tidyverse)
library(fs)
library(yardstick)
library(ComplexHeatmap)
library(khroma)
library(igraph)
library(limma)
library(ggrepel)
library(patchwork)
source("scripts/styling.R")



# Load data ---------------------------------------------------------------


## P-NET ----

ORIGINAL_SEED_PNET <- "234_20080808"

reactome_names <-
  read_tsv(
    "pnet_data/original/ReactomePathways.txt",
    col_names = c("reactome_id", "node", "species")) %>%
  filter(species == "Homo sapiens") %>%
  select(!species) %>%
  distinct(node, .keep_all = TRUE)

node_importance <-
  dir_ls(
    glob = str_glue("data/*/*/node_importance.csv"),
    recurse = TRUE
  ) %>%
  map(
    \(file) {
      read_csv(
        file,
        show_col_types = FALSE,
      )
    }
  ) %>%
  list_rbind(names_to = "file") %>%
  separate_wider_regex(
    file,
    c("data/", experiment = ".+", "/", seed = ".+", "/.+"),
  ) %>%
  rename(node = ...1) %>%
  left_join(reactome_names, by = "node", relationship = "many-to-one") %>%
  mutate(reactome_id = if_else(layer == 1, node, reactome_id)) %>%
  select(experiment, seed, layer, coef_combined, reactome_id)

predictions <-
  dir_ls(
    glob = "data/*/*/predictions_test.csv",
    recurse = TRUE
  ) %>%
  map(
    \(file) {
      read_csv(
        file,
        skip = 1,
        col_names = c("sample_id", "predicted", "p_predicted", "truth"),
        col_types = "ccdc"
      )
    }
  ) %>%
  list_rbind(names_to = "file") %>%
  separate_wider_regex(
    file,
    c("data/", experiment = ".+", "/", seed = ".+", "/.+"),
  ) %>%
  mutate(
    truth =
      truth %>%
      fct_recode(metastatic = "1", primary = "0") %>%
      fct_relevel("metastatic", "primary"),
    predicted =
      predicted %>%
      fct_recode(metastatic = "1.0", primary = "0.0") %>%
      fct_relevel("metastatic", "primary")
  )

# only include nodes with connectivity; therefore, load link weight matrices of
# NN and search for nonzero entries
pnet_edges <-
  map_dfr(
    0:6,
    function(layer) {
      df <- read_csv(
        str_glue("data/pnet_original/234_20080808/link_weights_{layer}.csv")
      )
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


## DTox ----

ORIGINAL_SEED_DTOX <- "0"

node_importance_dtox <-
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
    c("data/dtox/", seed = ".+", "/.+"),
  )

predictions_dtox <-
  dir_ls(glob = "data/dtox/*/test_labels.csv", recurse = TRUE) %>%
  map(read_csv) %>%
  list_rbind(names_to = "file") %>%
  separate_wider_regex(
    file,
    c("data/dtox/", seed = ".+", "/.+")
  ) %>%
  mutate(
    experiment = "dtox",
    truth =
      if_else(truth > 0.9, "one", "zero") %>%
      as_factor() %>%
      fct_relevel("one", "zero")
  ) %>%
  rename(p_predicted = predicted) %>%
  relocate(experiment)



# Figure 1 ----------------------------------------------------------------

## b ----

plot_roc <- function(predictions,
                     experiment,
                     original_seed = ORIGINAL_SEED_PNET) {
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
        map(roc_curve, truth, p_predicted)
    ) %>%
    unnest(metrics)

  ggplot(plot_data, aes(1 - specificity, sensitivity, group = seed)) +
    geom_path(alpha = .1, linewidth = BASE_LINEWIDTH) +
    geom_path(
      data = plot_data %>% filter(seed == original_seed),
      linewidth = BASE_LINEWIDTH,
      color = ORIGINAL_SEED_COLOR
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_fixed() +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_roc(predictions, "pnet_original")
ggsave_publication("1b_roc", width = 4, height = 4)


## c ----

plot_robustness <- function(top_nodes = 10, experiment = "pnet_original") {
  node_importance <-
    node_importance %>%
    filter(experiment == {{experiment}})

  top_nodes_per_layer <-
    node_importance %>%
    filter(seed == ORIGINAL_SEED_PNET) %>%
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
      data = node_importance %>% filter(seed == ORIGINAL_SEED_PNET),
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
ggsave_publication("1c_robustness", width = 18, height = 4)


## d ----

plot_changes <- function(seed = "28_28", top_nodes = c(4, 43)) {
  plot_data <-
    node_importance %>%
    filter(experiment == "pnet_original", layer == 1)

  top_nodes <-
    plot_data %>%
    filter(seed == ORIGINAL_SEED_PNET) %>%
    slice_max(coef_combined, n = top_nodes[2], with_ties = FALSE) %>%
    slice_min(coef_combined, n = top_nodes[2] - top_nodes[1] + 1, with_ties = FALSE) %>%
    arrange(desc(coef_combined))

  plot_data <-
    plot_data %>%
    filter(seed %in% c(ORIGINAL_SEED_PNET, {{seed}})) %>%
    semi_join(top_nodes, by = "reactome_id") %>%
    mutate(reactome_id = factor(reactome_id,
                                levels = top_nodes$reactome_id)) %>%
    select(seed, reactome_id, coef_combined)

  ggplot(plot_data, aes(reactome_id, coef_combined)) +
    geom_segment(
      data = plot_data %>%
        mutate(
          seed = fct_recode(seed,
                            original = ORIGINAL_SEED_PNET,
                            replicate = {{seed}})
        ) %>%
        pivot_wider(names_from = seed, values_from = coef_combined),
      aes(xend = reactome_id, yend = replicate, y = original),
      arrow = arrow(length = unit(1, "mm"), ends = "last", type = "closed"),
      linewidth = BASE_LINEWIDTH
    ) +
    geom_point(
      data =
        plot_data %>%
        filter(seed == ORIGINAL_SEED_PNET),
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

plot_changes()
ggsave_publication("1d_changes", width = 8, height = 4)



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

plot_bias("pnet_deterministic")
ggsave_publication("2b_scores", width = 14, height = 4)


## c ----

plot_roc(predictions, "pnet_deterministic")
ggsave_publication("2c_roc", width = 4, height = 4)


## d ----

plot_bias_comparison <- function(experiment, top_nodes = 5) {
  node_importance <-
    node_importance %>%
    filter(experiment %in% c("pnet_original", {{experiment}}))

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
      linewidth = BASE_LINEWIDTH * 1.5,
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

plot_bias_comparison("pnet_deterministic")
ggsave_publication("2d_comparison", width = 14, height = 4)



# Figure 3 ----------------------------------------------------------------

## b ----

plot_bias("pnet_shuffled")
ggsave_publication("3b_scores", width = 14, height = 4)


## c ----

plot_roc(predictions, "pnet_shuffled")
ggsave_publication("3c_roc", width = 4, height = 4)


## d ----

plot_bias_comparison("pnet_shuffled")
ggsave_publication("3d_comparison", width = 14, height = 4)



# Figure 4 ----------------------------------------------------------------

plot_cor_heatmap <- function(layer = NULL,
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
      experiment %in% c("pnet_original", "pnet_deterministic", "pnet_shuffled")
    ) %>%
    unite(experiment, seed, col = "exp_seed", sep = "+") %>%
    select(exp_seed, reactome_id, coef_combined) %>%
    pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
    select(!reactome_id) %>%
    cor()

  distance <- as.dist(1 - corr_mat)

  col_metadata <-
    tibble(exp_seed = colnames(corr_mat)) %>%
    separate(exp_seed, into = c("experiment", "seed"), sep = "\\+")

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

    clustering_distance_rows = distance,
    clustering_distance_columns = distance,
    row_dend_gp = gpar(lwd = 0.5),
    row_title = "seeds",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

    width = unit(heatmap_size, "mm"),
    height = unit(heatmap_size, "mm"),
    border = FALSE,

    show_column_dend = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,

    left_annotation = rowAnnotation(
      experiment = col_metadata$experiment,
      col = list(experiment = EXPERIMENT_COLORS),
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

(p <- plot_cor_heatmap())
ggsave_publication("4_cor_heatmap", plot = p,
                   width = 9, height = 5.5, type = "png")



# Figure 5 ----------------------------------------------------------------

## a ----

cor_data_pearson <-
  node_importance %>%
  filter(
    experiment %in% c("pnet_original", "pnet_deterministic", "pnet_shuffled")
  ) %>%
  left_join(graph_stats, by = c("reactome_id", "layer")) %>%
  group_split(experiment, seed, layer) %>%
  map_dfr(
    ~bind_cols(
      distinct(.x, experiment, seed, layer),
      cor(
        .x$coef_combined,
        select(.x, indegree:betweenness),
        method = "pearson"
      ) %>%
        as_tibble()
    )
  )

plot_network_heatmap <- function(cor_data) {
  mat <-
    cor_data %>%
    group_by(experiment, layer) %>%
    summarise(across(indegree:betweenness, mean, na.rm = TRUE)) %>%
    pivot_longer(indegree:betweenness, names_to = "measure") %>%
    unite(experiment, experiment, measure, sep = ":") %>%
    pivot_wider(names_from = layer) %>%
    column_to_rownames("experiment") %>%
    as.matrix() %>%
    t()

  mat_min <- min(mat, na.rm = TRUE)
  mat_max <- max(mat, na.rm = TRUE)

  col_metadata <-
    tibble(colname = colnames(mat)) %>%
    separate(colname, into = c("experiment", "measure"), sep = ":")

  colnames(mat) <- col_metadata$measure

  Heatmap(
    mat,
    col = circlize::colorRamp2(
      seq(mat_min, mat_max, length.out = 9),
      color("tokyo", reverse = TRUE)(9),
    ),
    na_col = "gray90",
    name = "mean correlation of\nnode importance scores\nand centrality measure",
    heatmap_legend_param = list(
      at = round(c(mat_min, mat_max), 2),
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
    column_title = "{recode(x[1], !!!EXPERIMENT_NAMES)}",
    column_title_side = "top",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),

    column_split = col_metadata$experiment,
    column_gap = unit(1, "pt"),

    width = unit(3, "mm") * ncol(mat),
    height = unit(3, "mm") * nrow(mat),
    border = FALSE,

    top_annotation = HeatmapAnnotation(
      experiment = col_metadata$experiment,
      col = list(experiment = EXPERIMENT_COLORS),
      show_annotation_name = TRUE,
      show_legend = FALSE,
      annotation_legend_param = list(
        experiment = list(
          title = "experiment",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      ),
      annotation_name_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
  ) %>%
    draw(
      column_title = "network measure",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
}

(p <- plot_network_heatmap(cor_data_pearson))
ggsave_publication("5a_corr_network_measure", plot = p,
                   width = 8, height = 4, type = "png")


## b ----

normalize <- function(x) {
  if (max(x) == min(x))
    rep.int(0.5, length(x))
  else
    (x - min(x)) / (max(x) - min(x))
}

plot_selected_scatterplots <- function() {
  node_importance %>%
    filter(
      experiment %in% c("pnet_original", "pnet_deterministic", "pnet_shuffled"),
      layer == 5
    ) %>%
    group_by(experiment, layer, reactome_id) %>%
    summarise(coef_combined = mean(coef_combined)) %>%
    left_join(graph_stats, by = c("reactome_id", "layer")) %>%
    mutate(
      reachability = normalize(reachability),
      coef_combined = normalize(coef_combined)
    ) %>%
    ggplot(aes(reachability, coef_combined, color = experiment)) +
    geom_point(size = .75, alpha = .5, shape = 16, show.legend = FALSE) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      show.legend = FALSE,
      linewidth = BASE_LINEWIDTH
    ) +
    scale_x_continuous(
      "reachability in layer 5 (relative)",
      limits = 0:1,
      breaks = 0:1
    ) +
    scale_y_continuous(
      "mean node importance (relative)",
      limits = 0:1,
      breaks = 0:1
    ) +
    scale_color_manual(values = EXPERIMENT_COLORS) +
    coord_fixed() +
    facet_wrap(vars(experiment), nrow = 1) +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_selected_scatterplots()
ggsave_publication("5b_reachability_vs_importance", width = 8, height = 4)



# Figure 6 ----------------------------------------------------------------

plot_correction <- function(cutoff = 1) {
  # matrix, rows = nodes, columns = experiments/seeds,
  # values = quantile-normalized node importances
  normalized_mat <-
    node_importance %>%
    filter(experiment %in% c("pnet_original", "pnet_shuffled")) %>%
    unite(experiment, seed, col = "exp_seed", sep = "/") %>%
    select(!layer) %>%
    pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
    column_to_rownames("reactome_id") %>%
    as.matrix() %>%
    normalizeQuantiles()

  idx_original <- str_starts(colnames(normalized_mat), "pnet_original")

  # data frame, one row per node
  # original_seed: node importance in original setup with original seed
  # original: average node importance in original setup
  # control: average node importance with shuffled labels
  plot_data <-
    node_importance %>%
    filter(str_starts(experiment, "pnet")) %>%
    distinct(layer, reactome_id) %>%
    left_join(reactome_names, by = "reactome_id") %>%
    mutate(
      original_seed = normalized_mat[, "pnet_original/234_20080808"],
      original = rowMeans(normalized_mat[, idx_original]),
      control = rowMeans(normalized_mat[, !idx_original]),
      positive = original - control > 0,
      node = coalesce(node, reactome_id)
    ) %>%
    mutate(
      .by = layer,
      rank_original_seed = rank(-original_seed)
    )

  color_limits <-
    plot_data %>%
    filter(control < cutoff, original_seed < cutoff) %>%
    {hexbin::hexbin(.$control, .$original, xbins = 100)} %>%
    slot("count") %>%
    range()

  ggplot(plot_data, aes(control, original)) +
    geom_hex(
      data = plot_data %>% filter(control < cutoff, original_seed < cutoff),
      bins = 100
    ) +
    geom_segment(
      data = plot_data %>% filter(control > cutoff | original_seed > cutoff),
      aes(yend = control, xend = control, color = positive),
      show.legend = FALSE,
      linewidth = BASE_LINEWIDTH
    ) +
    geom_point(
      data = plot_data %>% filter(control > cutoff | original_seed > cutoff),
      aes(color = positive),
      show.legend = FALSE,
      size = 0.25
    ) +
    geom_abline(linewidth = BASE_LINEWIDTH) +
    geom_text_repel(
      data = plot_data %>% filter(rank_original_seed <= 5),
      aes(label = node),
      min.segment.length = 0,
      size = BASE_TEXT_SIZE_MM,
      segment.size = BASE_LINEWIDTH,
      segment.color = "gray70"
    ) +
    scale_color_manual(values = c("blue", "red")) +
    scale_fill_gradient(
      low = "grey",
      high = "black",
      limits = color_limits,
      breaks = color_limits,
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      ),
    ) +
    facet_wrap(vars(layer), ncol = 3) +
    xlab("average normalized importance score (control, i.e., shuffled labels)") +
    ylab("average normalized importance score (original setup)") +
    theme_pub() +
    theme(
      panel.grid = element_blank()
    )
}

plot_correction()
ggsave_publication("6_correction_approach", width = 18, height = 12)



# Figure S1 ---------------------------------------------------------------

left_join(
  node_importance %>%
    filter(experiment == "pnet_original", seed == ORIGINAL_SEED_PNET),
  node_importance %>%
    filter(experiment == "pnet_original", seed != ORIGINAL_SEED_PNET) %>%
    summarise(
      .by = c(layer, reactome_id),
      coef_combined_avg = mean(coef_combined)
    ),
  by = c("layer", "reactome_id")
) %>%
  ggplot(aes(coef_combined - coef_combined_avg, color = factor(layer))) +
  stat_ecdf() +
  xlab("difference between node importance (original seed)\n and average node importance (other seeds)") +
  ylab("cumulative density") +
  scale_color_brewer("layer", palette = "Dark2") +
  theme_pub() +
  theme(
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.position = c(.85, .35),
    panel.grid = element_blank()
  )
ggsave_publication("S1_delta_node_importance_distribution", width = 5, height = 4)



# Figure S2 ---------------------------------------------------------------

## a ----

plot_roc(predictions_dtox, "dtox", ORIGINAL_SEED_DTOX)
ggsave_publication("S2a_roc", width = 4, height = 4)


## b (boxblot) ----

corr_data_dtox <-
  map(
    unique(node_importance_dtox$compound),
    \(compound) {
      node_importance_dtox %>%
        filter(compound == {{compound}}) %>%
        select(!compound) %>%
        column_to_rownames("seed") %>%
        as.matrix() %>%
        t() %>%
        cor() %>%
        as_tibble(rownames = "seed_1") %>%
        pivot_longer(
          !seed_1,
          names_to = "seed_2",
          values_to = "correlation"
        ) %>%
        mutate(compound = {{compound}}, .before = 1)
    }
  ) %>%
  list_rbind()

corr_data_dtox %>%
  filter(seed_2 > seed_1) %>%
  mutate(compound = fct_reorder(compound, -correlation, median)) %>%
  ggplot(aes(compound, correlation)) +
  geom_boxplot(
    outlier.size = 0.1,
    outlier.alpha = 0.1,
    linewidth = BASE_LINEWIDTH
  ) +
  ylab("correlation of replicate networks") +
  theme_pub(TRUE) +
  theme(panel.grid = element_blank())
ggsave_publication("S2b_correlations", width = 18, height = 6)


## b (heatmaps) ----

plot_dtox_heatmap <- function(compound,
                              show_legends = TRUE,
                              heatmap_size = 30) {
  corr_mat <-
    corr_data_dtox %>%
    filter(compound == {{compound}}) %>%
    select(!compound) %>%
    pivot_wider(names_from = seed_2, values_from = correlation) %>%
    column_to_rownames("seed_1") %>%
    as.matrix()

  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(-1, 1, length.out = 9),
      color("roma", reverse = TRUE)(9),
    ),

    show_heatmap_legend = show_legends,
    name = "correlation of\nnode importance\nscores",
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
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

    show_row_names = FALSE,
    row_title = "seeds",
    row_title_side = "right",
    column_title = compound
  )
}

(p <- plot_dtox_heatmap("CID_150311"))
ggsave_publication("S2b_heatmap_150311",
                   plot = p, height = 4, width = 6)

(p <- plot_dtox_heatmap("CID_5281162"))
ggsave_publication("S2b_heatmap_5281162",
                   plot = p, height = 4, width = 6)



# Figure S3 ---------------------------------------------------------------

walk(
  1:6,
  function(l) {
    p <- plot_cor_heatmap(layer = l, show_legends = FALSE, heatmap_size = 40)
    ggsave_publication(str_glue("S3{letters[l]}_cor_heatmap_layer_{l}"),
                       plot = p, width = 6, height = 4.5, type = "png")
  }
)



# Figure S4 ---------------------------------------------------------------

label_middle <- function(labels) {
  middle_index <-
    tibble(labels = labels) %>%
    mutate(r = row_number()) %>%
    group_by(labels) %>%
    summarise(r = ceiling(mean(r))) %>%
    pull(r)
  replace(character(length(labels)), middle_index, labels[middle_index])
}

plot_importance_vs_measure <- function() {
  plot_data <-
    node_importance %>%
    filter(
      experiment %in% c("pnet_original", "pnet_deterministic", "pnet_shuffled")
    ) %>%
    group_by(experiment, layer, reactome_id) %>%
    summarise(coef_combined = mean(coef_combined)) %>%
    left_join(graph_stats, by = c("reactome_id", "layer")) %>%
    mutate(across(coef_combined:betweenness, normalize)) %>%
    pivot_longer(
      indegree:betweenness,
      names_to = "measure",
      names_transform =
        list(measure = fct_inorder)
    ) %>%
    mutate(experiment_long = recode(experiment, !!!EXPERIMENT_NAMES))

  ggplot(plot_data, aes(value, coef_combined, color = experiment)) +
    geom_point(size = .1, alpha = 1, shape = 16, show.legend = FALSE) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      show.legend = FALSE,
      linewidth = BASE_LINEWIDTH
    ) +
    scale_x_continuous(
      "value of network measure (relative)",
      limits = 0:1,
      breaks = 0:1,
      sec.axis = facet_title("experiment and network measure")
    ) +
    scale_y_continuous(
      "value of mean node importance (relative)",
      limits = 0:1,
      breaks = 0:1,
      sec.axis = facet_title("layer")
    ) +
    scale_color_manual(values = EXPERIMENT_COLORS) +
    coord_fixed() +
    facet_grid(
      vars(layer),
      vars(experiment_long, measure),
      labeller = labeller(experiment_long = label_middle)
    ) +
    theme_pub() +
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(1, "mm"),
      strip.clip = "off"
    )
}

plot_importance_vs_measure()
ggsave_publication("S4_importance_vs_measure",
                   width = 28, height = 12, type = "png")



# Figure S5 ---------------------------------------------------------------

graph_stats %>%
  filter(layer %>% between(2, 6)) %>%
  ggplot(aes(reachability, betweenness)) +
  geom_point(alpha = .25, size = .5, shape = 16) +
  scale_x_continuous(sec.axis = facet_title("layer")) +
  facet_wrap(vars(layer), nrow = 1, scales = "free") +
  theme_pub() +
  theme(panel.grid = element_blank())

ggsave_publication("S5_reachability_vs_betweenness",
                   height = 4, width = 18, type = "png")



# Figure S6 ---------------------------------------------------------------

## a ----

wrap_plots(
  plot_roc(predictions, "mskimpact_nsclc_original"),
  plot_roc(predictions, "mskimpact_bc_original"),
  plot_roc(predictions, "mskimpact_cc_original"),
  plot_roc(predictions, "mskimpact_pc_original"),
  plot_roc(predictions, "mskimpact_nsclc_shuffled"),
  plot_roc(predictions, "mskimpact_bc_shuffled"),
  plot_roc(predictions, "mskimpact_cc_shuffled"),
  plot_roc(predictions, "mskimpact_pc_shuffled"),
  nrow = 2
)
ggsave_publication("S6a_mskimpact_roc", width = 18, height = 9)


## b ----

plot_msk_heatmap <- function() {
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

  corr_data <-
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
    select(!reactome_id)

  corr_mat_biased <-
    corr_data %>%
    select(ends_with("biased")) %>%
    cor(use = "pairwise.complete.obs")

  corr_mat_corrected <-
    corr_data %>%
    select(ends_with("corrected")) %>%
    cor(use = "pairwise.complete.obs")

  color_limits <- range(corr_mat_biased, corr_mat_corrected)

  plot_heatmap <- function(corr_mat,
                           color_limits,
                           col_title = "",
                           show_legends = TRUE,
                           heatmap_size = 20) {
    cor_range <- range(corr_mat[lower.tri(corr_mat)]) %>% round(2)
    message("Correlation range (", col_title, "): ",
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
      row_labels = EXPERIMENT_NAMES[rownames(corr_mat)],
      row_title = "cancer",
      column_title = col_title
    )
  }

  (p <- plot_heatmap(corr_mat_biased, color_limits,
                     col_title = "original (biased)", show_legends = FALSE))
  ggsave_publication("S6b_msk_correlation_biased",
                     plot = p, height = 3, width = 6)

  (p <- plot_heatmap(corr_mat_corrected, color_limits, col_title = "corrected"))
  ggsave_publication("S6b_msk_correlation_corrected", plot = p, height = 3, width = 6)
}

plot_msk_heatmap()



# Values in main text -----------------------------------------------------

## P-NET importance score correlation ----

node_importance %>%
  filter(
    experiment %in% c("pnet_original", "pnet_deterministic", "pnet_shuffled")
  ) %>%
  unite(experiment, seed, col = "exp_seed", sep = "+") %>%
  select(exp_seed, reactome_id, coef_combined) %>%
  pivot_wider(names_from = exp_seed, values_from = coef_combined) %>%
  select(!reactome_id) %>%
  cor() %>%
  magrittr::inset(!lower.tri(.), NA) %>%
  as_tibble(rownames = "from") %>%
  pivot_longer(
    !from,
    names_to = c("to_exp", "to_seed"),
    names_sep = "\\+",
    values_to = "corr"
  ) %>%
  separate(from, into = c("from_exp", "from_seed"), sep = "\\+") %>%
  filter(!is.na(corr)) %>%
  group_by(from_exp, to_exp) %>%
  summarise(mean_corr = mean(corr) %>% round(2))



## DTox correlation extrema ----

corr_data_dtox %>%
  filter(seed_1 != seed_2) %>%
  pull(correlation) %>%
  range() %>%
