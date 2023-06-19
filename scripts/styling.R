library(tidyverse)
library(fs)
library(ComplexHeatmap)



# Common definitions ------------------------------------------------------

EXPERIMENT_NAMES <- c(
  "pnet_original" = "original setup",
  "pnet_deterministic" = "deterministic inputs",
  "pnet_shuffled" = "shuffled labels"
)

EXPERIMENT_COLORS <- c(
  "pnet_original" = "gray50",
  "pnet_deterministic" = "#3182bd",
  "pnet_shuffled" = "#e69f00",
  "mskimpact" = "#66a61e",
  "mskimpact_pc" = "#66a61e"
)

ORIGINAL_SEED_COLOR <- "red"



# ggplot functions --------------------------------------------------------

BASE_TEXT_SIZE_MM = 1.76  # mm, corresponds to 5 pt, use in geom_text()
BASE_TEXT_SIZE_PT = 5 # pt, use in theme()
BASE_LINEWIDTH = 0.25 # pt
BASE_BOXPLOT_SIZE = 0.5

#' Common theme for figures in the publication.
#'
#' This theme bases upon `theme_bw()` and ensures
#' - common line widths of `BASE_LINEWIDTH`
#' - common text sizes of `BASE_TEXT_SIZE_PT`
#' - a uniform plot margin of 1 mm
#' - a medium strip text, an empty strip background, and
#'   1 mm padding between strip text and panel
#'
#' @param rotate_x_labels If `TRUE`, rotate x-axis tick labels by 90°.
#' @param ... Other parameters passed to `theme_bw()`.
#'
#' @return A theme object.
theme_pub <- function(rotate_x_labels = FALSE, ...){
  res <-
    theme_bw(...) +
    theme(
      line = element_line(size = BASE_LINEWIDTH),
      axis.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      axis.title = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.background = element_blank(),
      legend.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.title = element_text(size = BASE_TEXT_SIZE_PT),
      panel.border = element_rect(size = BASE_LINEWIDTH * 2),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      strip.background = element_blank(),
      strip.text = element_text(
        color = "black",
        size = BASE_TEXT_SIZE_PT
      ),
      strip.text.x = element_text(margin = margin(b = 1, unit = "mm")),
      strip.text.y = element_text(margin = margin(l = 1, unit = "mm")),
      plot.title = element_text(size = BASE_TEXT_SIZE_PT, face = "bold")
    )

  if (rotate_x_labels)
    res <-
      res +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  res
}


#' Save a publication-quality plot.
#'
#' @param filename Filename, will be saved in subfolder `plots/final`.
#' @param plot Plot object.
#' @param width Width in cm.
#' @param height Height in cm.
#' @param type Type of image file.
#' @param dpi Resolution.
#' @param ... Other parameters passed to the plotting function.
ggsave_publication <- function(filename,
                               plot = NULL,
                               width = 4,
                               height = 4,
                               type = "pdf",
                               dpi = 1200,
                               ...) {
  filename <- str_glue("plots/{filename}.{type}")
  filename %>%
    path_dir() %>%
    dir_create()

  if (is.null(plot)) {
    # if last_plot() is available, use ggsave()
    ggsave(
      filename,
      dpi = dpi,
      units = "cm",
      limitsize = FALSE,
      width = width,
      height = height,
      ...
    )
  } else {
    # for non-ggplot objects, use the base R functions directly;
    # only png and pdf backends are supported
    if (type == "png") {
      png(
        filename,
        res = dpi,
        units = "cm",
        width = width,
        height = height,
        ...
      )
    } else if (type == "pdf") {
      pdf(
        filename,
        width = width / 2.54,  # dimensions for pdf() must be inches
        height = height / 2.54,
        ...
      )
    } else {
      stop("Type", type, "cannot be saved.")
    }

    print(plot)
    dev.off()
  }
}

# shorthand for adding a facet title as secondary axis
# use like scale_x_continuous (sec.axis = facet_title("..."))
facet_title <- function(name) {
  dup_axis(name = name, breaks = NULL, labels = NULL)
}



# Heatmap appearance ------------------------------------------------------

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(2, "mm"),
  heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_border = FALSE
)
