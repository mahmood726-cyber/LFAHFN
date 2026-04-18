# lfahfn-plots.R
# Visualization Functions for LFAHFN

#' Network Structure Plot
#' @param x LFAHFN object or data frame
#' @param ... Additional arguments
#' @export
lfahfn_network_plot <- function(x, ...) {
  fit <- if (inherits(x, "lfahfn")) x$main_nma else x
  netmeta::netgraph(fit, ...)
}

#' Forest Plot for LFAHFN
#' @param x LFAHFN object
#' @param ... Additional arguments
#' @export
lfahfn_forest <- function(x, ...) {
  if (!inherits(x, "lfahfn")) stop("x must be a 'lfahfn' object")
  meta::forest(x$main_nma, ...)
}

#' Plot Rankogram
#' @param x LFAHFN object
#' @export
lfahfn_rankogram <- function(x) {
  if (!inherits(x, "lfahfn")) stop("x must be a 'lfahfn' object")
  if (is.null(x$diagnostics$ranking)) stop("Ranking data not available")
  
  treatment <- pscore <- NULL 
  df <- data.frame(
    treatment = names(x$diagnostics$ranking$Pscore.random),
    pscore = as.numeric(x$diagnostics$ranking$Pscore.random)
  )
  df <- df[order(df$pscore, decreasing = TRUE), ]
  df$treatment <- factor(df$treatment, levels = df$treatment)
  
  ggplot2::ggplot(df, ggplot2::aes(x = treatment, y = pscore, fill = pscore)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_gradient(low = "#ffeda0", high = "#feb24c") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Treatment Ranking (P-scores)", y = "P-score", x = "Treatment") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot LFAHFN Results
#' @param x LFAHFN object
#' @param type Type of plot
#' @param ... Additional arguments
#' @export
plot.lfahfn <- function(x, type = c("forest", "network", "funnel", "rank", "transitivity"), ...) {
  type <- match.arg(type)
  switch(type,
    forest = lfahfn_forest(x, ...),
    network = lfahfn_network_plot(x, ...),
    funnel = meta::funnel(x$main_nma, ...),
    rank = lfahfn_rankogram(x),
    transitivity = {
      args <- list(...)
      if (is.null(args$covariate)) stop("Must specify 'covariate' for transitivity plot")
      lfahfn_transitivity_plot(x$data, args$covariate)
    }
  )
}

#' Transitivity Check Plot
#' @param data NMA data frame
#' @param covariate Name of the covariate
#' @export
lfahfn_transitivity_plot <- function(data, covariate) {
  comparison <- NULL
  data$comparison <- paste(data$treat1, "vs", data$treat2)
  ggplot2::ggplot(data, ggplot2::aes(x = comparison, y = rlang::.data[[covariate]], fill = comparison)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
    ggplot2::labs(title = paste("Transitivity Check:", covariate), x = "Comparison", y = covariate)
}
