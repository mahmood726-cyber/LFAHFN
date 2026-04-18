# lfahfn-component.R
# Component Network Meta-Analysis (CNMA)

#' Run Component Network Meta-Analysis
#' @param nma_fit A standard NMA fit object from netmeta()
#' @param sep Separator used to denote combinations in treatment names
#' @export
lfahfn_component_model <- function(nma_fit, sep = "+") {
  if (!inherits(nma_fit, "netmeta")) stop("nma_fit must be a 'netmeta' object")
  msg("Running Component Network Meta-Analysis...")
  
  if (!any(grepl(sep, nma_fit$trts, fixed = TRUE))) {
    msg("No combination treatments detected using separator: ", sep)
    return(NULL)
  }
  
  comp_model <- tryCatch({
    netmeta::netcomb(nma_fit, sep = sep)
  }, error = function(e) { NULL })
  
  if (!is.null(comp_model)) msg("Component Analysis complete.")
  comp_model
}
