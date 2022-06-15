library(tidyverse)
library(fs)



# Common definitions ------------------------------------------------------

EXPERIMENT_COLORS <- c(
  "default" = "black",
  "correlated" = "#3182bd",
  "scrambled_labels_balanced" = "#e69f00"
)

ORIGINAL_SEED_COLOR <- "red"


# ggplot functions --------------------------------------------------------

BASE_TEXT_SIZE_MM = 1.76  # mm, corresponds to 5 pt, use in geom_text()
BASE_TEXT_SIZE_PT = 5 # pt, use in theme()
BASE_LINE_SIZE = 0.25 # pt
BASE_BOXPLOT_SIZE = 0.5

#' Common theme for figures in the publication.
#'
#' This theme bases upon `theme_bw()` and ensures
#' - common line widths of `BASE_LINE_SIZE`
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
      line = element_line(size = BASE_LINE_SIZE),
      axis.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      axis.title = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.background = element_blank(),
      legend.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.title = element_text(size = BASE_TEXT_SIZE_PT),
      panel.border = element_rect(size = BASE_LINE_SIZE * 2),
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
#' @param type Type of image file.
#' @param dpi Resolution.
#' @param ... Other parameters passed to the plotting function.
ggsave_publication <- function(filename,
                               type = "pdf",
                               dpi = 1200,
                               ...) {
  filename <- str_glue("plots/final/{filename}.{type}")
  filename %>%
    path_dir() %>%
    dir_create()

  ggsave(
    filename,
    dpi = dpi,
    units = "cm",
    limitsize = FALSE,
    ...
  )
}
