
# Do we construct the correct graph? --------------------------------------

gs_old <-
  map_dfr(
    0:5,
    function(layer) {
      df <- read_csv(str_glue("pnet_prostate_paper/analysis/extracted/graph_stats_h{layer}.csv"))
      colnames(df)[1] <- "reactome_id"
      df %>%
        mutate(layer = {{layer}})
    }
  )

graph_stats %>%
  separate(reactome_id, into = c("reactome_id", "layer"), sep = "\\:", convert = TRUE) %>%
  left_join(gs_old, by = c("reactome_id", "layer")) %>%
  filter(degree.x != degree.y) %>%
  # filter(outdegree != fanout) %>%
  # View() %>%
  {.}
