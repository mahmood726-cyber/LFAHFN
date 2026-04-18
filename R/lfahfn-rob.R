# lfahfn-rob.R
# Risk of Bias Network Propagation (CINeMA-style)

#' Propagate Risk of Bias Across the Network
#'
#' @param nma_fit A standard NMA fit object from netmeta()
#' @param rob_data A data frame with columns `studlab` and `rob_score`
#' @return A data frame containing the expected network risk of bias for each comparison
#' @export
lfahfn_rob_flow <- function(nma_fit, rob_data) {
  
  if (!inherits(nma_fit, "netmeta")) stop("nma_fit must be a 'netmeta' object")
  if (!all(c("studlab", "rob_score") %in% names(rob_data))) stop("rob_data must contain 'studlab' and 'rob_score' columns")
  
  msg("Propagating Risk of Bias through network evidence flow...")
  
  contrib <- lfahfn_contribution_matrix(nma_fit)
  if (is.null(contrib)) return(NULL)
  
  c_matrix <- if (!is.null(contrib$random)) contrib$random else contrib$fixed
  direct_comps <- colnames(c_matrix)
  nma_data <- nma_fit$data
  
  merged_data <- merge(nma_data, rob_data, by = "studlab", all.x = TRUE)
  if (any(is.na(merged_data$rob_score))) {
    warning("Some studies in the network are missing RoB scores. Imputing median.")
    merged_data$rob_score[is.na(merged_data$rob_score)] <- median(merged_data$rob_score, na.rm=TRUE)
  }
  
  merged_data$comp_str <- paste(
    pmin(as.character(merged_data$treat1), as.character(merged_data$treat2)),
    pmax(as.character(merged_data$treat1), as.character(merged_data$treat2)),
    sep = ":"
  )
  
  direct_rob <- aggregate(rob_score ~ comp_str, data = merged_data, FUN = mean)
  direct_rob_vector <- numeric(length(direct_comps))
  names(direct_rob_vector) <- direct_comps
  
  for (i in seq_along(direct_comps)) {
    comp <- direct_comps[i]
    if (comp %in% direct_rob$comp_str) {
      direct_rob_vector[i] <- direct_rob$rob_score[direct_rob$comp_str == comp]
    } else {
      direct_rob_vector[i] <- mean(merged_data$rob_score)
    }
  }
  
  network_rob <- c_matrix %*% direct_rob_vector
  res <- data.frame(
    comparison = rownames(network_rob),
    expected_rob_score = as.numeric(network_rob),
    stringsAsFactors = FALSE
  )
  
  res$rob_category <- cut(res$expected_rob_score, 
                          breaks = c(-Inf, 1.5, 2.5, Inf), 
                          labels = c("Low", "Moderate", "High"))
  
  msg("Risk of Bias network propagation complete.")
  res
}
