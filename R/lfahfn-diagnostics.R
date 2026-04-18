# lfahfn-diagnostics.R
# Diagnostic and Sensitivity Analysis Functions

#' Run Comprehensive Diagnostics
#'
#' @param nma_fit Network meta-analysis object
#' @param data Original data
#' @param config Configuration
#' @return List of diagnostic results
#' @export
lfahfn_diagnostics <- function(nma_fit, data, config = setup_lfahfn()) {

  diag <- list()
  msg("Running diagnostic analyses...")

  # Inconsistency assessment
  if (config$enable_ume %||% TRUE) {
    diag$inconsistency <- tryCatch({
      netmeta::netsplit(nma_fit)
    }, error = function(e) NULL)
  }

  # Global Inconsistency assessment (Design-by-Treatment Interaction)
  diag$global_inconsistency <- tryCatch({
    netmeta::netmeasures(nma_fit)
  }, error = function(e) NULL)
  
  if (!is.null(nma_fit$Q.inconsistency)) {
    diag$inconsistency_stats <- list(
      Q = nma_fit$Q.inconsistency,
      df = nma_fit$df.Q.inconsistency,
      p_value = 1 - stats::pchisq(nma_fit$Q.inconsistency, nma_fit$df.Q.inconsistency)
    )
  }

  # Heterogeneity
  diag$heterogeneity <- list(
    tau = nma_fit$tau,
    I2 = nma_fit$I2.random,
    Q = nma_fit$Q
  )

  # Ranking
  diag$ranking <- tryCatch({
    netmeta::netrank(nma_fit)
  }, error = function(e) NULL)

  # Publication Bias
  diag$bias <- lfahfn_publication_bias(nma_fit, data)
  
  # Formal network-wide test for small-study effects (Egger's test for NMA)
  diag$small_study_effects <- tryCatch({
    meta::metabias(nma_fit, method.bias = "Egger")
  }, error = function(e) { NULL })

  class(diag) <- "lfahfn_diagnostics"
  diag
}

#' Publication Bias Assessment
#' @param nma_fit NMA object
#' @param data Data
#' @return Publication bias results
#' @export
lfahfn_publication_bias <- function(nma_fit, data) {
  if (!inherits(nma_fit, "netmeta")) stop("nma_fit must be a 'netmeta' object")
  bias <- list()
  msg("Assessing publication bias...")
  
  bias$funnel <- suppressWarnings(tryCatch({
    meta::funnel(nma_fit, order = NULL, method.bias = "Egger")
  }, error = function(e) {
    NULL
  }))
  
  if (!is.null(bias$funnel)) {
    bias$egger <- list(
      intercept = bias$funnel$intercept,
      se = bias$funnel$se.intercept,
      p_value = bias$funnel$p.value
    )
  }
  bias
}
