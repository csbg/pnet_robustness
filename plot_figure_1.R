library(tidyverse)
library(fs)
source("styling.R")


# Load data ---------------------------------------------------------------

loaded_experiments <- c(
  "default",
  "correlated",
  "scrambled_labels_balanced"
)

original_seed <- "234_20080808"

loaded_seeds <- c(
  original_seed,
  str_c(0:49, "_", 0:49)
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
  filter(experiment %in% loaded_experiments, seed %in% loaded_seeds)



# Figure 1a ---------------------------------------------------------------

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
    ggtitle(paste("Training of replicate networks",
                  "reveals variability of interpretations")) +
    facet_wrap(vars(layer), scales = "free_x", nrow = 1) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
}

plot_robustness()
ggsave_publication("1a_robustness", width = 18, height = 4)



# Figure 1b ---------------------------------------------------------------

plot_bias <- function(experiment, title = "", top_nodes = 5) {
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
    ggtitle(title) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
}

plot_bias("correlated", title = "Artificial inputs reveal network biases")
ggsave_publication("1b_correlated", width = 18, height = 4)



# Figure 1c ---------------------------------------------------------------

plot_bias_comparison <- function(experiment, title = "", top_nodes = 5) {
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
    ggtitle(title) +
    theme_pub() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
}

plot_bias_comparison(
  "correlated",
  title = paste("Comparison of importance measures from control",
                "and real inputs quantifies network biases")
)
ggsave_publication("1c_correlated_comparison", width = 18, height = 4)



# Figure 1d ---------------------------------------------------------------

plot_bias(
  "scrambled_labels_balanced",
  title = "Randomized target labels provide a second control for network biases"
)
ggsave_publication("1d_scrambled", width = 18, height = 4)



# Figure 1e ---------------------------------------------------------------

plot_bias_comparison(
  "scrambled_labels_balanced",
  title = paste("Comparison of importance measures from randomized",
                "and real targets quantifies network biases")
)
ggsave_publication("1e_scrambled_comparison", width = 18, height = 4)
