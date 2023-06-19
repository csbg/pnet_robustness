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
