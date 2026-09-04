# Small shared utilities.

#' Order experiment identifiers numerically
#'
#' Experiment identifiers are strings of the form "IGT5", "IGT27", "IGT96".
#' Sorting them as characters puts IGT10 before IGT5, which scrambles every
#' axis and legend ordered by experiment. This returns a factor ordered by the
#' numeric part instead.
#'
#' @param igt Character or factor vector of IGT identifiers.
#' @return A factor with levels in numeric order.
ReorderIGT <- function(igt) {
    extract_numeric <- function(x) as.numeric(gsub("\\D", "", x))
    ordered_vector <- unique(igt)[order(extract_numeric(unique(igt)))] %>%
        as.character()
    factor(igt, levels = ordered_vector)
}
