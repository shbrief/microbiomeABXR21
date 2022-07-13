#' Plot summary of multiple ABX exposure status
#'
#' @param x A table. The number of participants for the different number of
#' ABX exposures.
#'
plotMultiABX <- function(x) {
    tbl <- stack(x)
    p <- ggplot(tbl, aes(x = ind, y = values)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_text(aes(label = values), vjust = -0.3, size = 3) +
        labs(x = "# of ABX exposures", y = "# of Participants") +
        theme_minimal()
    p
}
