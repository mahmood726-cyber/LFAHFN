# lfahfn-priors.R
# Data-Driven Empirical Priors for NMA Heterogeneity

#' Generate Empirical Prior for Heterogeneity (Tau)
#' @param outcome_type Outcome type: "objective", "semi-objective", or "subjective"
#' @param intervention_type Comparison type: "pharma_vs_pharma", "pharma_vs_placebo", etc.
#' @export
lfahfn_empirical_prior <- function(outcome_type = c("objective", "semi-objective", "subjective"),
                                 intervention_type = c("pharma_vs_pharma", "pharma_vs_placebo", "non_pharma", "any")) {

  outcome_type <- match.arg(outcome_type)
  intervention_type <- match.arg(intervention_type)

  # Default fallback
  mu <- -2.56; sigma <- 1.74
  
  if (outcome_type == "objective") {
    if (intervention_type == "pharma_vs_placebo") { mu <- -3.44; sigma <- 1.52 }
    else if (intervention_type == "pharma_vs_pharma") { mu <- -4.13; sigma <- 1.48 }
    else { mu <- -3.11; sigma <- 1.60 }
  } else if (outcome_type == "semi-objective") {
    if (intervention_type == "pharma_vs_placebo") { mu <- -2.71; sigma <- 1.54 }
    else if (intervention_type == "pharma_vs_pharma") { mu <- -3.40; sigma <- 1.51 }
    else { mu <- -2.38; sigma <- 1.62 }
  } else if (outcome_type == "subjective") {
    if (intervention_type == "pharma_vs_placebo") { mu <- -2.23; sigma <- 1.55 }
    else if (intervention_type == "pharma_vs_pharma") { mu <- -2.92; sigma <- 1.52 }
    else { mu <- -1.90; sigma <- 1.63 }
  }
  
  precision <- 1 / (sigma^2)
  list(mu = mu, sigma = sigma, precision = precision, gemtc_prior = sprintf("dlnorm(%.2f, %.4f)", mu, precision))
}

#' Apply Informative Priors to Bayesian NMA
#' @param data NMA data
#' @param config Configuration
#' @export
lfahfn_apply_prior <- function(data, config) {
  if (is.null(config$empirical_prior_outcome) || is.null(config$empirical_prior_intervention)) return(NULL)
  lfahfn_empirical_prior(outcome_type = config$empirical_prior_outcome, intervention_type = config$empirical_prior_intervention)
}
