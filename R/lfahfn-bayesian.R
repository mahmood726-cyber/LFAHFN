# lfahfn-bayesian.R
# Bayesian Network Meta-Analysis Functions

#' Run Bayesian Network Meta-Analysis
#'
#' Performs Bayesian NMA using JAGS through the gemtc package.
#' Includes convergence diagnostics (R-hat).
#'
#' @param data NMA data frame
#' @param config Configuration object
#' @return List with Bayesian NMA results
#' @export
lfahfn_bayesian <- function(data, config) {
  if (!requireNamespace("gemtc", quietly = TRUE)) {
    message("Bayesian analysis skipped: 'gemtc' package not installed")
    return(NULL)
  }

  msg("Preparing data for Bayesian analysis...")
  studies <- unique(data$studlab)
  treatments <- unique(c(data$treat1, data$treat2))

  network <- tryCatch({
    gemtc::mtc.network(data.re = data.frame(
      study = data$studlab,
      treatment = data$treat1,
      diff = data$TE,
      std.err = data$seTE,
      stringsAsFactors = FALSE
    ), treatments = data.frame(
      id = treatments,
      description = treatments,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) { NULL })

  if (is.null(network)) return(NULL)
  msg("Running Bayesian MCMC...")

  model <- tryCatch({
    gemtc::mtc.model(
      network = network,
      type = "consistency",
      likelihood = config$bayes_likelihood %||% "normal",
      link = config$bayes_link %||% "identity",
      n.chain = config$bayes_chains %||% 3,
      linearModel = "random"
    )
  }, error = function(e) { NULL })

  if (is.null(model)) return(NULL)

  results <- tryCatch({
    gemtc::mtc.run(
      model,
      n.adapt = config$bayes_warmup %||% 5000,
      n.iter = config$bayes_iter %||% 20000,
      thin = max(1, (config$bayes_iter %||% 20000) / 1000)
    )
  }, error = function(e) { NULL })

  if (is.null(results)) return(NULL)
  msg("Bayesian analysis completed")

  # Convergence check (Potential Scale Reduction Factor)
  psrf <- tryCatch({
    coda::gelman.diag(results$samples, multivariate = FALSE)
  }, error = function(e) { NULL })

  if (!is.null(psrf)) {
    max_rhat <- max(psrf$psrf[,1], na.rm = TRUE)
    if (max_rhat > 1.1) {
      warning(sprintf("MCMC convergence warning: Max R-hat is %.2f (recommended < 1.1)", max_rhat))
    }
  }

  list(
    model = model,
    results = results,
    summary = summary(results),
    dic = tryCatch(gemtc::mtc.deviance(results), error = function(e) NULL),
    convergence = psrf
  )
}

#' Bayesian Node-Splitting Analysis
#' @param data NMA data
#' @param config Configuration
#' @export
lfahfn_nodesplit_bayesian <- function(data, config) {
  if (!requireNamespace("gemtc", quietly = TRUE)) return(NULL)
  msg("Preparing node-splitting Bayesian analysis...")
  
  treatments <- unique(c(data$treat1, data$treat2))
  network <- gemtc::mtc.network(
    data.re = data.frame(
      study = data$studlab,
      treatment = data$treat1,
      diff = data$TE,
      std.err = data$seTE,
      stringsAsFactors = FALSE
    ),
    treatments = data.frame(
      id = treatments,
      description = treatments,
      stringsAsFactors = FALSE
    )
  )
  
  msg("Running node-splitting sampling...")
  results <- tryCatch({
    gemtc::mtc.nodesplit(
      network,
      type = "consistency",
      likelihood = config$bayes_likelihood %||% "normal",
      link = config$bayes_link %||% "identity",
      n.chain = config$bayes_chains %||% 3,
      n.adapt = config$bayes_warmup %||% 5000,
      n.iter = config$bayes_iter %||% 20000,
      thin = max(1, (config$bayes_iter %||% 20000) / 1000)
    )
  }, error = function(e) { NULL })
  
  if (!is.null(results)) msg("Node-splitting analysis completed")
  results
}
