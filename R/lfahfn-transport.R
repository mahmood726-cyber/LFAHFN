# lfahfn-transport.R
# Transportability Analysis Functions

#' Calculate Transportability Weights
#'
#' @param data Study data with characteristics
#' @param target_population Target population characteristics
#' @param metric Distance metric
#' @param kernel Kernel function
#' @param truncation Minimum weight threshold
#' @export
calculate_transport_weights <- function(data,
                                       target_population,
                                       metric = c("mahalanobis", "euclidean"),
                                       kernel = c("gaussian", "tricube"),
                                       truncation = 0.01) {

  metric <- match.arg(metric)
  kernel <- match.arg(kernel)
  if (is.null(target_population)) return(rep(1, nrow(data)))
  vars <- intersect(names(target_population), names(data))
  if (length(vars) == 0) return(rep(1, nrow(data)))

  study_chars <- as.matrix(data[, vars, drop = FALSE])
  target_chars <- unlist(target_population[vars])
  complete_rows <- complete.cases(study_chars)
  weights <- rep(1, nrow(data))

  if (sum(complete_rows) > 0) {
    study_chars_complete <- study_chars[complete_rows, , drop = FALSE]
    if (metric == "mahalanobis") {
      cov_mat <- stats::cov(study_chars_complete)
      distances <- if (det(cov_mat) > 1e-10) stats::mahalanobis(study_chars_complete, target_chars, cov_mat) else sqrt(rowSums((sweep(study_chars_complete, 2, target_chars, "-"))^2))
    } else {
      distances <- sqrt(rowSums((sweep(study_chars_complete, 2, target_chars, "-"))^2))
    }

    if (kernel == "gaussian") {
      h <- stats::median(distances[distances > 0])
      if (is.na(h) || h == 0) h <- 1
      kernel_weights <- exp(-0.5 * (distances / h)^2)
    } else {
      h <- stats::quantile(distances, 0.75)
      if (is.na(h) || h == 0) h <- max(distances)
      kernel_weights <- ifelse(distances < h, (1 - (distances/h)^3)^3, 0)
    }

    kernel_weights[kernel_weights < truncation] <- truncation
    kernel_weights <- kernel_weights / mean(kernel_weights)
    weights[complete_rows] <- kernel_weights
  }
  weights
}

#' Apply Transportability Weights to NMA
#' @param data NMA data
#' @param weights Transportability weights
#' @param config Configuration
#' @export
apply_transport_weights <- function(data, weights, config) {
  data$seTE_weighted <- data$seTE / sqrt(weights)
  netmeta::netmeta(TE = data$TE, seTE = data$seTE_weighted, treat1 = data$treat1, treat2 = data$treat2, 
                   studlab = data$studlab, data = data, sm = config$sm, common = FALSE, random = TRUE,
                   reference.group = config$ref_treatment %||% names(sort(table(c(data$treat1, data$treat2)), decreasing = TRUE))[1])
}
