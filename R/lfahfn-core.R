# lfahfn-core.R
# Main LFAHFN Package Functions

#' @keywords internal
"_PACKAGE"

#' @import netmeta
#' @importFrom meta forest funnel metabias
#' @import coda
#' @importFrom stats complete.cases qnorm aggregate median quantile rnorm runif sd na.omit pchisq qlogis plogis cov AIC coef predict mahalanobis
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import tibble
#' @importFrom rlang .data %||%
#' @importFrom utils capture.output write.csv installed.packages sessionInfo combn
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom grDevices dev.off pdf png svg
#' @importFrom graphics plot
NULL

#' Setup LFAHFN Configuration
#'
#' @param sm Summary measure (default: "HR")
#' @param use_transport Use transportability weighting
#' @param use_bayesian Run Bayesian analysis
#' @param run_metareg Run meta-regression
#' @param metareg_covariates Covariates for meta-regression
#' @param run_diagnostics Run diagnostic analyses
#' @param run_ntsa Run Network Trial Sequential Analysis (NTSA)
#' @param run_component Run Component Network Meta-Analysis
#' @param component_sep Separator for component NMA (default: "+")
#' @param run_ml_covariates Use Machine Learning for covariate selection
#' @param use_robust_bayes Use robust t-distribution for Bayesian NMA
#' @param empirical_prior_outcome Outcome type for empirical prior
#' @param empirical_prior_intervention Intervention type for empirical prior
#' @param run_geometry Calculate network geometry and contribution matrix
#' @param run_rob_flow Propagate Risk of Bias through network evidence flow
#' @param rob_data Data frame containing Risk of Bias scores
#' @param run_multiverse Run Multiverse Stability Analysis
#' @param run_fragility Calculate Network Fragility Index
#' @param rwe_data Optional Registry/Observational data for cross-design synthesis
#' @param ref_treatment Reference treatment
#' @param enable_ume Enable UME inconsistency model
#' @param bayes_chains Number of MCMC chains
#' @param bayes_iter MCMC iterations
#' @param bayes_warmup MCMC warmup iterations
#' @param bayes_likelihood Likelihood for Bayesian model
#' @param bayes_link Link function for Bayesian model
#' @param seed Random seed
#' @param quiet Suppress progress messages
#' @return Configuration list with class "lfahfn_config"
#' @export
setup_lfahfn <- function(sm = "HR",
                        use_transport = FALSE,
                        use_bayesian = FALSE,
                        run_metareg = FALSE,
                        metareg_covariates = NULL,
                        run_diagnostics = TRUE,
                        run_ntsa = FALSE,
                        run_component = FALSE,
                        component_sep = "+",
                        run_ml_covariates = FALSE,
                        use_robust_bayes = FALSE,
                        empirical_prior_outcome = NULL,
                        empirical_prior_intervention = NULL,
                        run_geometry = FALSE,
                        run_rob_flow = FALSE,
                        rob_data = NULL,
                        run_multiverse = FALSE,
                        run_fragility = FALSE,
                        rwe_data = NULL,
                        ref_treatment = NULL,
                        enable_ume = TRUE,
                        bayes_chains = 3,
                        bayes_iter = 20000,
                        bayes_warmup = 5000,
                        bayes_likelihood = "normal",
                        bayes_link = "identity",
                        seed = 42,
                        quiet = FALSE) {

  config <- as.list(environment())
  class(config) <- "lfahfn_config"
  options(lfahfn.quiet = quiet)
  config
}

#' Run LFAHFN Analysis
#'
#' @param data Data frame with NMA data
#' @param target_population Target population for transportability
#' @param config Configuration from setup_lfahfn()
#' @return LFAHFN results object with class "lfahfn"
#' @export
run_lfahfn_analysis <- function(data,
                               target_population = NULL,
                               config = setup_lfahfn()) {

  set.seed(config$seed)
  validate_lfahfn_input(data)
  data <- lfahfn_clean_data(data)

  # Auto-select reference
  ref_treatment <- names(sort(table(c(data$treat1, data$treat2)), decreasing = TRUE))[1]
  msg("Reference treatment:", ref_treatment)

  results <- list(data = data, config = config, ref_treatment = ref_treatment)

  msg("Running primary network meta-analysis...")
  results$main_nma <- netmeta::netmeta(
    TE = data$TE, seTE = data$seTE, treat1 = data$treat1, treat2 = data$treat2,
    studlab = data$studlab, data = data, sm = config$sm,
    common = FALSE, random = TRUE, reference.group = ref_treatment
  )

  # Feature Integration
  if (config$use_transport && !is.null(target_population)) {
    results$transport_weights <- calculate_transport_weights(data, target_population)
    results$transport_nma <- apply_transport_weights(data, results$transport_weights, config)
  }

  if (config$use_bayesian) results$bayesian <- lfahfn_bayesian(data, config)
  if (config$run_metareg) results$metareg <- lfahfn_metareg(data, config$metareg_covariates, config)
  if (config$run_diagnostics) results$diagnostics <- lfahfn_diagnostics(results$main_nma, data, config)
  if (config$run_ntsa) results$ntsa <- lfahfn_ntsa(results$main_nma)
  if (config$run_component) results$component <- lfahfn_component_model(results$main_nma, sep = config$component_sep)
  if (config$run_ml_covariates) results$ml_covariates <- lfahfn_ml_covariate_selection(data)
  if (config$use_robust_bayes) results$robust_bayesian <- lfahfn_robust_bayesian(data, config)
  if (!is.null(config$empirical_prior_outcome)) results$empirical_prior <- lfahfn_apply_prior(data, config)
  if (config$run_geometry) results$geometry <- lfahfn_contribution_matrix(results$main_nma)
  if (config$run_rob_flow && !is.null(config$rob_data)) {
    results$rob_flow <- lfahfn_rob_flow(results$main_nma, config$rob_data)
  }

  # Novel Methodology Integration
  if (config$run_multiverse) results$multiverse <- lfahfn_multiverse(data, config)
  if (config$run_fragility) results$fragility <- lfahfn_fragility(results$main_nma)
  if (!is.null(config$rwe_data)) {
    results$cross_design <- lfahfn_cross_design(data, config$rwe_data, config = config)
  }

  class(results) <- "lfahfn"
  msg("Analysis complete")
  results
}

#' Validate LFAHFN Input
#' @param data Input data
#' @export
validate_lfahfn_input <- function(data) {
  need <- c("studlab", "treat1", "treat2", "TE", "seTE")
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  if (any(data$seTE <= 0, na.rm = TRUE)) stop("seTE must be strictly positive")
  invisible(TRUE)
}

#' Clean LFAHFN Data
#' @param data Input data frame
#' @export
lfahfn_clean_data <- function(data) {
  data[!is.na(data$TE) & !is.na(data$seTE) & data$seTE > 0, ]
}

#' Null coalescing operator
#' @name null-coalesce
#' @rdname null-coalesce
#' @keywords internal
#' @export
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

#' @importFrom rlang .data
NULL

#' Print progress message
#' @keywords internal
msg <- function(...) {
  if (!getOption("lfahfn.quiet", FALSE)) {
    message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
  }
}
