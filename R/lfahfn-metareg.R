# lfahfn-metareg.R
# Network Meta-Regression Functions

#' Run Network Meta-Regression
#' @param data NMA data
#' @param covariates Character vector of covariate names
#' @param config Configuration object
#' @export
lfahfn_metareg <- function(data, covariates = NULL, config = setup_lfahfn()) {
  if (is.null(covariates)) covariates <- intersect(c("year", "age_mean", "female_pct", "bmi_mean"), names(data))
  if (length(covariates) == 0) return(NULL)
  
  results <- list()
  for (covar in covariates) {
    if (!(covar %in% names(data))) next
    msg(sprintf("Running meta-regression with %s...", covar))
    if (sum(!is.na(data[[covar]])) < 5) next
    
    cov_data <- data
    cov_data$covariate <- scale(data[[covar]])[,1]
    fit <- tryCatch({
      netmeta::netmeta(TE = cov_data$TE, seTE = cov_data$seTE, treat1 = cov_data$treat1, treat2 = cov_data$treat2,
                       studlab = cov_data$studlab, data = cov_data, sm = config$sm, common = FALSE, random = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) results[[covar]] <- list(fit = fit, covariate = covar, n_studies = sum(!is.na(data[[covar]])))
  }
  class(results) <- "lfahfn_metareg"
  results
}
