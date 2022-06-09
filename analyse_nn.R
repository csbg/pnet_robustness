library(yardstick)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(fs)


ggsave_default <- function(filename,
                           type = "png",
                           plot = NULL,
                           width = 297,
                           height = 210,
                           crop = TRUE,
                           ...) {
  if (is.null(filename))
    return()

  filename <- stringr::str_glue("plots/{filename}.{type}")
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()

  if (is.null(plot)) {
    ggplot2::ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
                    width = width, height = height, ...)
  } else {
    rlang::exec(type, filename, res = 300, units = "mm",
                width = width, height = height,  ...)
    print(plot)
    dev.off()
  }

  if (crop)
    knitr::plot_crop(filename)

  invisible(filename)
}



# Load data ---------------------------------------------------------------

loaded_experiments <- c(
  "default",
  "correlated",
  "dropout_none",
  "dropout_high",
  "scrambled_labels",
  "scrambled_labels_balanced",
  "scrambled_labels_seeds",
  "scrambled_features_0.5_seed_0",
  "scrambled_features_0.5_seed_1",
  "scrambled_features_0.05_seed_0",
  "scrambled_features_0.05_seed_1",
  "scrambled_features_0.001_seed_0",
  NULL
)

loaded_seeds <- c(
  "234_20080808",
  str_c(0:9, "_", 0:9),
  -1:9
)

node_importance <-
  dir_ls(
    glob = str_glue("data/*/*/node_importance_graph_adjusted.csv"),
    recurse = TRUE
  ) %>%
  map_dfr(
    read_csv,
    .id = "file"
  ) %>%
  rename(node = X1) %>%
  extract(
    file,
    into = c("experiment", "seed"),
    regex = str_glue("data/(.+)/(.+)/"),
    convert = TRUE
  ) %>%
  filter(experiment %in% loaded_experiments, seed %in% loaded_seeds)


predictions <-
  dir_ls(
    regex = "P-net_ALL_testing.csv",
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
    obs = X1,
    truth = fct_recode(y, class_1 = "1", class_0 = "0"),
    class_1 = pred_scores,
    predicted = fct_recode(pred, class_1 = "1.0", class_0 = "0.0")
  ) %>%
  filter(experiment %in% loaded_experiments, seed %in% loaded_seeds)



# Node importance ---------------------------------------------------------

plot_node_importance <- function(layer,
                                 top_nodes = Inf,
                                 experiments = loaded_experiments,
                                 show_raw_data = FALSE,
                                 show_original_seed = FALSE) {
  node_importance <-
    node_importance %>%
    filter(experiment %in% experiments) %>%
    mutate(experiment = factor(experiment, levels = experiments))

  top_nodes <-
    node_importance %>%
    filter(layer == {{layer}}, experiment == experiments[1]) %>%
    group_by(node) %>%
    summarise(coef = median(coef)) %>%
    slice_max(coef, n = top_nodes, with_ties = FALSE) %>%
    arrange(coef) %>%
    pull(node)

  if (show_raw_data) {
    raw_data <- geom_point(
      alpha = .25,
      size = 1,
      position = position_jitter(width = .25, seed = 1)
    )
  } else {
    layer_raw_data <- NULL
  }

  if (show_original_seed) {
    layer_original_seed <- geom_point(
      data =
        node_importance %>%
        filter(seed == "234_20080808", layer == {{layer}}, node %in% top_nodes),
      aes(fill = experiment),
      shape = 21,
      position = position_dodge(width = .75),
      show.legend = FALSE
    )
  } else {
    layer_original_seed <- NULL
  }

  node_importance %>%
    filter(layer == {{layer}}, node %in% top_nodes) %>%
    mutate(node = factor(node, levels = top_nodes)) %>%
    ggplot(aes(node, coef)) +
    geom_boxplot(
      aes(fill = experiment),
      outlier.shape = 21,
      outlier.alpha = .25,
      key_glyph = "rect"
    ) +
    layer_raw_data +
    layer_original_seed +
    ylab("node importance score") +
    coord_flip() +
    ggtitle(str_glue("Layer {layer}")) +
    theme_bw() +
    theme()
}


plot_node_importance_all_layers <- function(experiments = loaded_experiments,
                                            top_nodes = 20) {
  node_importance <-
    node_importance %>%
    filter(experiment %in% experiments) %>%
    mutate(experiment = factor(experiment, levels = experiments))

  top_nodes_per_layer <-
    node_importance %>%
    filter(experiment == experiments[1]) %>%
    group_by(layer, node) %>%
    summarise(coef = median(coef)) %>%
    slice_max(coef, n = top_nodes, with_ties = FALSE) %>%
    arrange(layer, desc(coef))

  node_importance %>%
    semi_join(top_nodes_per_layer, by = c("layer", "node")) %>%
    mutate(node = factor(node, levels = top_nodes_per_layer$node)) %>%
    ggplot(aes(node, coef)) +
    geom_boxplot(
      aes(fill = experiment),
      outlier.shape = 21,
      outlier.alpha = .25,
      key_glyph = "rect"
    ) +
    ylab("node importance score") +
    facet_wrap(vars(layer), scales = "free", nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}


## Dropout ----

EXP_DROPOUT <- c(
  "default",
  "dropout_high",
  "dropout_none"
)

plot_node_importance_all_layers(EXP_DROPOUT)
ggsave_default("node_importance_dropout", width = 600, height = 200)

plot_node_importance(6, experiments = EXP_DROPOUT)
ggsave_default("node_importance_6_dropout", width = 250, height = 200)

plot_node_importance(5, 20, experiments = EXP_DROPOUT)
ggsave_default("node_importance_5_dropout", width = 200)

plot_node_importance(1, 20, experiments = EXP_DROPOUT)
ggsave_default("node_importance_1_dropout", width = 200)


## Correlated ----

EXP_CORRELATED <- c(
  "default",
  "correlated"
)

plot_node_importance_all_layers(EXP_CORRELATED)
ggsave_default("node_importance_correlated", width = 500, height = 200)

plot_node_importance(6, experiments = EXP_CORRELATED)
ggsave_default("node_importance_6_correlated", width = 250, height = 200)

plot_node_importance(5, 20, experiments = EXP_CORRELATED)
ggsave_default("node_importance_5_correlated", width = 200)

plot_node_importance(1, 20, experiments = EXP_CORRELATED)
ggsave_default("node_importance_1_correlated", width = 200)


## Scrambled labels ----

EXP_LABELS <- c(
  "default",
  "scrambled_labels",
  "scrambled_labels_balanced",
  "scrambled_labels_seeds"
)

plot_node_importance_all_layers(EXP_LABELS[c(1, 2)])
ggsave_default("node_importance_scrambled_labels", width = 500, height = 200)

plot_node_importance_all_layers(EXP_LABELS[c(1, 3)])
ggsave_default("node_importance_scrambled_labels_balanced", width = 500, height = 200)

plot_node_importance_all_layers(EXP_LABELS[c(1, 4)])
ggsave_default("node_importance_scrambled_labels_seeds", width = 500, height = 200)


plot_node_importance(6, experiments = EXP_LABELS)
ggsave_default("node_importance_6_labels", width = 250, height = 300)

plot_node_importance(5, 20, experiments = EXP_LABELS)
ggsave_default("node_importance_5_labels", width = 200)

plot_node_importance(1, 20, experiments = EXP_LABELS)
ggsave_default("node_importance_1_labels", width = 200)


## Scrambled features ----

EXP_FEATURES_A <- c(
  "default",
  "scrambled_features_0.5_seed_0",
  "scrambled_features_0.5_seed_1"
)

plot_node_importance_all_layers(EXP_FEATURES_A[c(1, 2)])
ggsave_default("node_importance_scrambled_features_0.5_seed_0", width = 500, height = 200)

plot_node_importance_all_layers(EXP_FEATURES_A[c(1, 3)])
ggsave_default("node_importance_scrambled_features_0.5_seed_1", width = 500, height = 200)

plot_node_importance(6, experiments = EXP_FEATURES_A)
ggsave_default("node_importance_6_features_0.5", width = 250, height = 200)

plot_node_importance(5, 20, experiments = EXP_FEATURES_A)
ggsave_default("node_importance_5_features_0.5", width = 200)

plot_node_importance(1, 20, experiments = EXP_FEATURES_A)
ggsave_default("node_importance_1_features_0.5", width = 200)


EXP_FEATURES_B <- c(
  "default",
  "scrambled_features_0.05_seed_0",
  "scrambled_features_0.05_seed_1"
)

plot_node_importance_all_layers(EXP_FEATURES_B[c(1, 2)])
ggsave_default("node_importance_scrambled_features_0.05_seed_0", width = 500, height = 200)

plot_node_importance_all_layers(EXP_FEATURES_B[c(1, 3)])
ggsave_default("node_importance_scrambled_features_0.05_seed_1", width = 500, height = 200)

plot_node_importance(6, experiments = EXP_FEATURES_B)
ggsave_default("node_importance_6_features_0.05", width = 250, height = 200)

plot_node_importance(5, 20, experiments = EXP_FEATURES_B)
ggsave_default("node_importance_5_features_0.05", width = 200)

plot_node_importance(1, 20, experiments = EXP_FEATURES_B)
ggsave_default("node_importance_1_features_0.05", width = 200)


EXP_FEATURES_C <- c(
  "default",
  "scrambled_features_0.001_seed_0"
)

plot_node_importance_all_layers(EXP_FEATURES_C)
ggsave_default("node_importance_scrambled_features_0.001", width = 500, height = 200)

plot_node_importance(6, experiments = EXP_FEATURES_C)
ggsave_default("node_importance_6_features_0.001", width = 250, height = 200)

plot_node_importance(5, 20, experiments = EXP_FEATURES_C)
ggsave_default("node_importance_5_features_0.001", width = 200)

plot_node_importance(1, 20, experiments = EXP_FEATURES_C)
ggsave_default("node_importance_1_features_0.001", width = 200)



# Metrics -----------------------------------------------------------------

multi_metric <- metric_set(accuracy, roc_auc, pr_auc, f_meas, precision, recall)

pred <-
  predictions %>%
  group_by(experiment, seed)

metrics <-
  group_keys(pred) %>%
  mutate(
    metrics =
      pred %>%
    group_split() %>%
    map(multi_metric, truth, class_1, estimate = predicted)
  ) %>%
  unnest(metrics) %>%
  select(!.estimator) %>%
  rename(metric = .metric, estimate = .estimate)


plot_metrics <- function(experiments = loaded_experiments,
                         show_original_seed = FALSE) {
  metrics <-
    metrics %>%
    filter(experiment %in% experiments) %>%
    mutate(experiment = factor(experiment, levels = experiments))

  if (show_original_seed) {
    layer_original_seed <- geom_point(
      data = metrics %>% filter(seed == "234_20080808"),
      aes(fill = experiment),
      shape = 21,
      position = position_dodge(width = .75)
    )
  } else {
    layer_original_seed <- NULL
  }

  ggplot(metrics, aes(metric, estimate)) +
    geom_boxplot(
      aes(fill = experiment),
      outlier.size = 1,
      outlier.alpha = 0.25
    ) +
    layer_original_seed +
    coord_flip() +
    theme_bw()
}

plot_metrics(EXP_DROPOUT)
ggsave_default("metrics_dropout", width = 200, height = 120)

plot_metrics(EXP_CORRELATED)
ggsave_default("metrics_correlated", width = 200, height = 120)

plot_metrics(EXP_LABELS)
ggsave_default("metrics_labels", width = 200, height = 120)

plot_metrics(unique(c(EXP_FEATURES_A, EXP_FEATURES_B, EXP_FEATURES_C)))
ggsave_default("metrics_labels", width = 200, height = 120)



# Curves ------------------------------------------------------------------

plot_curve <- function(curve_fun, x, y) {
  group_keys(pred) %>%
    mutate(
      metrics =
        pred %>%
        group_split() %>%
        map(curve_fun, truth, class_1)
    ) %>%
    unnest(metrics) %>%
    ggplot(aes({{x}}, {{y}}, group = seed)) +
    geom_path(color = "blue", alpha = .1) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    coord_fixed() +
    facet_wrap(vars(experiment)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}

plot_curve(roc_curve, 1 - specificity, sensitivity)
ggsave_default("curve_roc", width = 250)

plot_curve(pr_curve, recall, precision)
ggsave_default("curve_pr", width = 250)




# Unused analyses ---------------------------------------------------------

## Load weights

# # only load weights for layers 3–6 due to memory limits
# # data for layer 1 is too big even if loaded separately :(
# link_weights <-
#   dir_ls(
#     regex = str_glue("data/{DATA_DIR}/.+/link_weights_[3-6]\\.csv"),
#     recurse = TRUE
#   ) %>%
#   map_dfr(
#     ~read_csv(.) %>%
#       rename(source = X1) %>%
#       pivot_longer(!source, names_to = "sink", values_to = "weight"),
#     .id = "file"
#   ) %>%
#   extract(
#     file,
#     into = c("seed", "layer"),
#     regex = str_glue("{DATA_DIR}/(.+)/link_weights_(.)"),
#     convert = TRUE
#   ) %>%
#   unite(source, sink, col = name, remove = FALSE)



## Clustering

# plot_corr_mat <- function(layer,
#                           type = c("node", "link"),
#                           abs_weights = FALSE) {
#   type <- match.arg(type)
#
#   if (type == "node")
#     corr_mat <-
#       node_importance %>%
#       filter(layer %in% {{layer}}, !near(coef, 0)) %>%
#       select(seed, node, coef) %>%
#       pivot_wider(names_from = seed, values_from = coef) %>%
#       column_to_rownames("node")
#   else
#     corr_mat <-
#       link_weights %>%
#       filter(layer == {{layer}}, !near(weight, 0)) %>%
#       select(seed, name, weight) %>%
#       pivot_wider(names_from = seed, values_from = weight) %>%
#       column_to_rownames("name")
#
#   corr_mat <-
#     corr_mat %>%
#     as.matrix() %>%
#     replace_na(0) %>%
#     cor(use = "pairwise.complete.obs")
#
#   if (abs_weights) {
#     corr_mat <- abs(corr_mat)
#     type <- "link_abs"
#   }
#
#   distance <- as.dist(1 - corr_mat)
#
#   p <- Heatmap(
#     corr_mat,
#     col = circlize::colorRamp2(
#       seq(min(corr_mat), max(corr_mat), length.out = 9),
#       brewer.pal(9, "YlGnBu"),
#     ),
#     name = "correlation",
#     heatmap_legend_param = list(
#       at = round(c(min(corr_mat), max(corr_mat)), 2)
#     ),
#
#     show_row_names = FALSE,
#     show_column_dend = FALSE,
#     show_column_names = FALSE,
#
#     clustering_distance_rows = distance,
#     clustering_distance_columns = distance,
#
#     width = unit(100, "mm"),
#     height = unit(100, "mm"),
#   )
#   ggsave_default(
#     str_glue("corr_{type}_{str_c(layer, collapse = '+')}"),
#     plot = p,
#     height = 300
#   )
#   p
# }
#
# 1:6 %>% walk(plot_corr_mat, "node")
# plot_corr_mat(1:6, "node")
# 3:6 %>% walk(plot_corr_mat, "link", abs_weights = FALSE)
# 3:6 %>% walk(plot_corr_mat, "link", abs_weights = TRUE)



## Weights

# plot_link_weights <- function(layer, top_links = Inf, abs_weights = FALSE) {
#   vis_data <-
#     link_weights %>%
#     filter(layer == {{layer}})
#
#   selected_edges <-
#     vis_data %>%
#     group_by(name) %>%
#     summarise(keep = sum(weight) > 0) %>%
#     filter(keep) %>%
#     pull(name)
#
#   vis_data <-
#     vis_data %>%
#     filter(name %in% selected_edges)
#
#   if (abs_weights) {
#     vis_data <-
#       vis_data %>%
#       mutate(weight = abs(weight))
#     sort_fun <- median
#   } else {
#     sort_fun <- sd
#   }
#
#   top_links <-
#     vis_data %>%
#     filter() %>%
#     group_by(name) %>%
#     summarise(weight = sort_fun(weight)) %>%
#     slice_max(weight, n = top_links, with_ties = FALSE) %>%
#     pull(name)
#
#   vis_data <-
#     vis_data %>%
#     filter(name %in% top_links) %>%
#     mutate(name = factor(name) %>% fct_reorder(weight, .fun = sort_fun))
#
#   ggplot(vis_data, aes(name, weight)) +
#     geom_violin(scale = "width") +
#     geom_point(
#       alpha = .25,
#       size = 1,
#       position = position_jitter(width = .25, seed = 1)
#     ) +
#     geom_point(
#       data =
#         vis_data %>% filter(seed == "234_20080808"),
#       color = "red"
#     ) +
#     xlab("edge") +
#     ylab(if_else(abs_weights, "absolute weight", "weight")) +
#     coord_flip() +
#     theme_bw()
# }
#
# plot_link_weights(6)
# ggsave_default("link_weights_6", height = 150, width = 200)
# plot_link_weights(6, 200, abs_weights = TRUE)
# ggsave_default("link_weights_6_abs", height = 150, width = 200)
#
# plot_link_weights(5, 100)
# ggsave_default("link_weights_5", height = 300, width = 200)
# plot_link_weights(5, 100, abs_weights = TRUE)
# ggsave_default("link_weights_5_abs", height = 300, width = 200)
#
# plot_link_weights(4, 100)
# ggsave_default("link_weights_4", height = 600, width = 200)
# plot_link_weights(4, 100, abs_weights = TRUE)
# ggsave_default("link_weights_4_abs", height = 600, width = 200)
#
# plot_link_weights(3, 100)
# ggsave_default("link_weights_3", height = 600, width = 200)
# plot_link_weights(3, 100, abs_weights = TRUE)
# ggsave_default("link_weights_3_abs", height = 600, width = 200)
