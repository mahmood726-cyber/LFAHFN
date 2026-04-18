# lfahfn-geometry.R
# Network Geometry and Contribution Matrix

#' Calculate Network Contribution Matrix
#' @param nma_fit A standard NMA fit object from netmeta()
#' @export
lfahfn_contribution_matrix <- function(nma_fit) {
  if (!inherits(nma_fit, "netmeta")) stop("nma_fit must be a 'netmeta' object")
  msg("Calculating network contribution matrix...")
  contrib <- tryCatch({
    netmeta::netcontrib(nma_fit)
  }, error = function(e) {
    warning("Failed to calculate contribution matrix: ", e$message)
    NULL
  })
  contrib
}
