# Three-axis alluvial plot with one cluster highlighted (Figure 10f).

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggalluvial)
})

#' Alluvial plot of cluster -> organ -> perturbation with one cluster picked out
#'
#' Every flow is drawn, but only those belonging to the highlighted cluster are
#' coloured and fully opaque; the rest are grey and translucent.
#'
#' @param df3_plot Data frame with `count`, `annotation_level2`,
#'   `organ_simplified`, `condition_detailed_simplified`, plus the `fill_key`
#'   and `alpha_val` columns that select the highlight.
#' @param mypal_alluv Named colours for `fill_key`.
#' @param alpha_range Opacity range mapped from `alpha_val`.
#' @param x_limits Axis order, left to right before the coordinate flip.
PlotAlluvialHighlight <- function(df3_plot,
                                  mypal_alluv,
                                  title = "",
                                  alpha_range = c(0.2, 1),
                                  x_limits = c("condition_detailed_simplified",
                                               "organ_simplified",
                                               "annotation_level2")) {
    ggplot(
        df3_plot,
        aes(
            y     = count,
            axis2 = organ_simplified,
            axis1 = annotation_level2,
            axis3 = condition_detailed_simplified
        )
    ) +
        geom_alluvium(aes(fill = fill_key, alpha = alpha_val),
                      width = 1 / 5, colour = NA) +
        geom_stratum(width = 1 / 5, fill = "white", color = "grey50") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)),
                  size = 3, angle = 90) +
        scale_x_discrete(limits = x_limits, expand = c(.05, .05)) +
        coord_flip() +
        scale_fill_manual(values = mypal_alluv) +
        scale_alpha(range = alpha_range, guide = "none") +
        ggtitle(title) +
        theme_minimal() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title       = element_blank(),
            axis.text        = element_blank(),
            axis.ticks       = element_blank()
        )
}
