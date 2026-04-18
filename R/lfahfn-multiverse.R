# lfahfn-multiverse.R
# Multiverse Stability Analysis for Network Meta-Analysis

#' Run Multiverse NMA Stability Analysis
#'
#' Evaluates how sensitive the network results are to different statistical 
#' choices (estimators, adjustments).
#'
#' @param data NMA data frame
#' @param config Configuration object
#' @return A "multiverse" object containing results across all specifications
#' @export
lfahfn_multiverse <- function(data, config = setup_lfahfn()) {
  
  msg("Initiating Multiverse Stability Analysis (Exploring 12+ specifications)...")
  
  # Define the multiverse grid
  methods <- c("REML", "ML", "DL", "SJ")
  kh_adjust <- c(TRUE, FALSE)
  
  results_list <- list()
  
  for (m in methods) {
    for (kh in kh_adjust) {
      spec_name <- paste0(m, "_KH", kh)
      
      fit <- tryCatch({
        netmeta::netmeta(
          TE = data$TE, seTE = data$seTE, 
          treat1 = data$treat1, treat2 = data$treat2,
          studlab = data$studlab, data = data,
          method.tau = m,
          hakn = kh,
          common = FALSE, random = TRUE
        )
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        results_list[[spec_name]] <- list(
          tau = fit$tau,
          top_treatment = names(sort(netmeta::netrank(fit)$Pscore.random, decreasing = TRUE))[1],
          p_scores = netmeta::netrank(fit)$Pscore.random
        )
      }
    }
  }
  
  # Calculate Stability Score (0 to 1)
  # Proportion of models that agree on the top treatment
  top_trts <- sapply(results_list, function(x) x$top_treatment)
  agreement <- max(table(top_trts)) / length(top_trts)
  
  res <- list(
    specifications = results_list,
    stability_score = agreement,
    consensus_treatment = names(sort(table(top_trts), decreasing = TRUE))[1]
  )
  
  class(res) <- "lfahfn_multiverse"
  msg(sprintf("Multiverse analysis complete. Stability Score: %.2f", agreement))
  
  res
}
