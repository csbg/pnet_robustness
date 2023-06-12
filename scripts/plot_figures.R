library(tidyverse)
library(fs)
library(yardstick)
library(ComplexHeatmap)
library(khroma)
library(igraph)
source("scripts/styling.R")



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
    show_col_types = FALSE,
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
  mutate(experiment = fct_recode(experiment, !!!EXPERIMENT_NAMES)) %>%
  ungroup() %>%
  select(experiment, seed, layer, coef_combined, reactome_id)

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
    experiment = fct_recode(experiment, !!!EXPERIMENT_NAMES),
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

## b ----

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

plot_roc("original setup")
ggsave_publication("1b_roc", width = 4, height = 4)


## c ----

plot_robustness <- function(top_nodes = 10) {
  node_importance <-
    node_importance %>%
    filter(experiment == "original setup")

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
ggsave_publication("1c_robustness", width = 18, height = 4)


## d ----

plot_changes <- function(seed = "28_28", top_nodes = c(4, 43)) {
  plot_data <-
    node_importance %>%
    filter(experiment == "original setup", layer == 1)

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
      linewidth = BASE_LINEWIDTH
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

plot_bias("deterministic inputs")
ggsave_publication("2b_scores", width = 14, height = 4)


## c ----

plot_roc("deterministic inputs")
ggsave_publication("2c_roc", width = 4, height = 4)


## d ----

plot_bias_comparison <- function(experiment, top_nodes = 5) {
  node_importance <-
    node_importance %>%
    filter(experiment %in% c("original setup", {{experiment}}))

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

plot_bias_comparison("deterministic inputs")
ggsave_publication("2d_comparison", width = 14, height = 4)



# Figure 3 ----------------------------------------------------------------

## b ----

plot_bias("shuffled labels")
ggsave_publication("3b_scores", width = 14, height = 4)


## c ----

plot_roc("shuffled labels")
ggsave_publication("3c_roc", width = 4, height = 4)


## d ----

plot_bias_comparison("shuffled labels")
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



# Figure S1 ---------------------------------------------------------------

walk(
  1:6,
  function(l) {
    p <- plot_cor_heatmap(layer = l, show_legends = FALSE, heatmap_size = 40)
    ggsave_publication(str_glue("S1{letters[l]}_cor_heatmap_layer_{l}"),
                       plot = p, width = 6, height = 4.5, type = "png")
  }
)



# Figure 5 ----------------------------------------------------------------

## a ----

cor_data_pearson <-
  node_importance %>%
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
    unite(experiment, experiment, measure) %>%
    pivot_wider(names_from = layer) %>%
    column_to_rownames("experiment") %>%
    as.matrix() %>%
    t()

  mat_min <- min(mat, na.rm = TRUE)
  mat_max <- max(mat, na.rm = TRUE)

  col_metadata <-
    tibble(colname = colnames(mat)) %>%
    separate(colname, into = c("experiment", "measure"), sep = "_")

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
    filter(layer == 5) %>%
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



# Figure S2 ---------------------------------------------------------------

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
    group_by(experiment, layer, reactome_id) %>%
    summarise(coef_combined = mean(coef_combined)) %>%
    left_join(graph_stats, by = c("reactome_id", "layer")) %>%
    mutate(across(coef_combined:betweenness, normalize)) %>%
    pivot_longer(
      indegree:betweenness,
      names_to = "measure",
      names_transform =
        list(measure = fct_inorder)
    )

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
      vars(experiment, measure),
      labeller = labeller(experiment = label_middle)
    ) +
    theme_pub() +
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(1, "mm"),
      strip.clip = "off"
    )
}

plot_importance_vs_measure()
ggsave_publication("S2_importance_vs_measure",
                   width = 28, height = 12, type = "png")



# Figure S3 ---------------------------------------------------------------

graph_stats %>%
  filter(layer %>% between(2, 6)) %>%
  ggplot(aes(reachability, betweenness)) +
  geom_point(alpha = .25, size = .5, shape = 16) +
  scale_x_continuous(sec.axis = facet_title("layer")) +
  facet_wrap(vars(layer), nrow = 1, scales = "free") +
  theme_pub() +
  theme(panel.grid = element_blank())

ggsave_publication("S3_reachability_vs_betweenness",
                   height = 4, width = 18, type = "png")



# Values in main text -----------------------------------------------------

## Importance score correlation ----

node_importance %>%
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
  summarise(mean_corr = mean(corr))

