# lfahfn-ml.R
# Machine Learning Augmented Network Meta-Analysis

#' Discover Non-Linear Covariates using Machine Learning
#'
#' @param data NMA dataset containing effect sizes (TE) and potential covariates
#' @param potential_covariates Vector of column names to evaluate
#' @param n_trees Number of trees in the forest
#' @export
lfahfn_ml_covariate_selection <- function(data, potential_covariates = NULL, n_trees = 500) {
  
  if (!requireNamespace("ranger", quietly = TRUE)) {
    warning("Package 'ranger' is required for ML covariate selection.")
    return(NULL)
  }
  
  if (is.null(potential_covariates)) {
    exclude_cols <- c("studlab", "treat1", "treat2", "TE", "seTE", "seTE_weighted")
    potential_covariates <- setdiff(names(data), exclude_cols)
    potential_covariates <- potential_covariates[sapply(data[potential_covariates], function(x) is.numeric(x) || is.factor(x))]
  }
  
  if (length(potential_covariates) == 0) return(NULL)
  msg(sprintf("Running Random Forest ML to identify key predictors among %d covariates...", length(potential_covariates)))
  
  ml_data <- data[, c("TE", potential_covariates), drop=FALSE]
  for (col in potential_covariates) {
    if (any(is.na(ml_data[[col]]))) {
      if (is.numeric(ml_data[[col]])) {
        ml_data[[col]][is.na(ml_data[[col]])] <- median(ml_data[[col]], na.rm = TRUE)
      } else {
        mode_val <- as.character(sort(table(ml_data[[col]]), decreasing=TRUE)[1])
        ml_data[[col]][is.na(ml_data[[col]])] <- mode_val
      }
    }
  }
  
  rf_fit <- tryCatch({
    ranger::ranger(formula = TE ~ ., data = ml_data, num.trees = n_trees, importance = 'permutation', respect.unordered.factors = TRUE)
  }, error = function(e) { NULL })
  
  if (is.null(rf_fit)) return(NULL)
  importance_scores <- rf_fit$variable.importance
  ranked_covariates <- sort(importance_scores, decreasing = TRUE)
  threshold <- mean(ranked_covariates) 
  selected <- names(ranked_covariates[ranked_covariates > threshold])
  
  list(model = rf_fit, importance = ranked_covariates, selected_covariates = selected, r_squared = rf_fit$r.squared)
}
