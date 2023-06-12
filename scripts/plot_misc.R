plot_roc("MSK-IMPACT")
ggsave_publication("X_msk_ROC", width = 4, height = 4)

plot_robustness(experiment = "MSK-IMPACT")
ggsave_publication("X_msk_robustness", width = 18, height = 4)

plot_bias("MSK-IMPACT")
ggsave_publication("X_msk_scores", width = 14, height = 4)

plot_bias_comparison("MSK-IMPACT")
ggsave_publication("X_msk_comparison", width = 14, height = 4)
