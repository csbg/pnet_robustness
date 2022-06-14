library(yardstick)
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

reactome_names <-
  read_tsv(
    "pnet_prostate_paper/_database/pathways/Reactome/ReactomePathways.txt",
    col_names = c("reactome_id", "node", "species")) %>%
  filter(species == "Homo sapiens") %>%
  select(!species)

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
  left_join(reactome_names, by = "node") %>%
  mutate(reactome_id = if_else(layer == 1, node, reactome_id)) %>%
  extract(
    file,
    into = c("experiment", "seed"),
    regex = str_glue("data/(.+)/(.+)/"),
    convert = TRUE
  ) %>%
  filter(experiment %in% loaded_experiments, seed %in% loaded_seeds) %>%
  group_by(experiment, seed, layer) %>%
  mutate(
    modified = coef_graph > mean(coef_graph) + 5 * sd(coef_graph),
  ) %>%
  ungroup()


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
    group_by(reactome_id) %>%
    summarise(coef_combined = median(coef_combined)) %>%
    slice_max(coef_combined, n = top_nodes, with_ties = FALSE) %>%
    arrange(coef_combined) %>%
    pull(reactome_id)

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
        filter(
          seed == "234_20080808",
          layer == {{layer}},
          reactome_id %in% top_nodes
        ),
      aes(fill = experiment),
      shape = 21,
      position = position_dodge(width = .75),
      show.legend = FALSE
    )
  } else {
    layer_original_seed <- NULL
  }

  node_importance %>%
    filter(layer == {{layer}}, reactome_id %in% top_nodes) %>%
    mutate(reactome_id = factor(reactome_id, levels = top_nodes)) %>%
    ggplot(aes(reactome_id, coef_combined)) +
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
    group_by(layer, reactome_id) %>%
    summarise(coef_combined = median(coef_combined)) %>%
    slice_max(coef_combined, n = top_nodes, with_ties = FALSE) %>%
    arrange(layer, desc(coef_combined))

  node_importance %>%
    semi_join(top_nodes_per_layer, by = c("layer", "reactome_id")) %>%
    mutate(reactome_id = factor(reactome_id,
                                levels = top_nodes_per_layer$reactome_id)) %>%
    ggplot(aes(reactome_id, coef_combined)) +
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
ggsave_default("node_importance_dropout", width = 600, height = 100)

plot_node_importance(6, experiments = EXP_DROPOUT)
ggsave_default("node_importance_6_dropout", width = 250)

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
ggsave_default("node_importance_correlated", width = 500, height = 100)

plot_node_importance(6, experiments = EXP_CORRELATED)
ggsave_default("node_importance_6_correlated", width = 250)

plot_node_importance(5, 20, experiments = EXP_CORRELATED)
ggsave_default("node_importance_5_correlated", width = 200)

plot_node_importance(1, 20, experiments = EXP_CORRELATED)
ggsave_default("node_importance_1_correlated", width = 200)


## Scrambled labels ----

EXP_LABELS <- c(
  "default",
  "scrambled_labels",
  "scrambled_labels_balanced"
)

plot_node_importance_all_layers(EXP_LABELS[c(1, 2)])
ggsave_default("node_importance_scrambled_labels", width = 500, height = 100)

plot_node_importance_all_layers(EXP_LABELS[c(1, 3)])
ggsave_default("node_importance_scrambled_labels_balanced",
               width = 500, height = 100)


plot_node_importance(6, experiments = EXP_LABELS)
ggsave_default("node_importance_6_labels", width = 250)

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
ggsave_default("node_importance_scrambled_features_0.5_seed_0",
               width = 500, height = 100)

plot_node_importance_all_layers(EXP_FEATURES_A[c(1, 3)])
ggsave_default("node_importance_scrambled_features_0.5_seed_1",
               width = 500, height = 100)

plot_node_importance(6, experiments = EXP_FEATURES_A)
ggsave_default("node_importance_6_features_0.5", width = 250)

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
ggsave_default("node_importance_scrambled_features_0.05_seed_0",
               width = 500, height = 100)

plot_node_importance_all_layers(EXP_FEATURES_B[c(1, 3)])
ggsave_default("node_importance_scrambled_features_0.05_seed_1",
               width = 500, height = 100)

plot_node_importance(6, experiments = EXP_FEATURES_B)
ggsave_default("node_importance_6_features_0.05", width = 250)

plot_node_importance(5, 20, experiments = EXP_FEATURES_B)
ggsave_default("node_importance_5_features_0.05", width = 200)

plot_node_importance(1, 20, experiments = EXP_FEATURES_B)
ggsave_default("node_importance_1_features_0.05", width = 200)


EXP_FEATURES_C <- c(
  "default",
  "scrambled_features_0.001_seed_0"
)

plot_node_importance_all_layers(EXP_FEATURES_C)
ggsave_default("node_importance_scrambled_features_0.001",
               width = 500, height = 100)

plot_node_importance(6, experiments = EXP_FEATURES_C)
ggsave_default("node_importance_6_features_0.001", width = 250)

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



# Node degrees ------------------------------------------------------------

node_importance %>%
  filter(experiment == "default", seed == "234_20080808") %>%
  ggplot(aes(coef_graph, coef)) +
  geom_point(aes(color = modified), alpha = .25) +
  facet_wrap(vars(layer), scales = "free_x")

# scrambled_labels_balanced
# correlated
# default
node_importance %>%
  # filter(experiment == "default") %>%
  filter(experiment %in% c("default",
                           "correlated",
                           "scrambled_labels_balanced",
                           "scrambled_features_0.5_seed_0")) %>%
  # filter(seed == "234_20080808") %>%
  filter(!modified) %>%
  ggplot(aes(coef_graph, coef)) +
  geom_point(alpha = .15) +
  geom_smooth(method = "lm") +
  # facet_wrap(vars(layer), scales = "free_x")
  facet_grid(vars(experiment), vars(layer), scales = "free_x")
ggsave_default("importance_vs_degree", height = 250)
