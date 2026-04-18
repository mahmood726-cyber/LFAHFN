# lfahfn-ntsa.R
# Network Trial Sequential Analysis (NTSA)

#' Network Trial Sequential Analysis (NTSA)
#' @param nma_fit A standard NMA fit object
#' @param alpha Overall Type I error rate
#' @param beta Overall Type II error rate (default: 0.20 for 80\% power)
#' @param mcid Minimum Clinically Important Difference
#' @export
lfahfn_ntsa <- function(nma_fit, alpha = 0.05, beta = 0.20, mcid = 0.2) {
  
  msg("Performing Network Trial Sequential Analysis (NTSA)...")
  if (!inherits(nma_fit, "netmeta")) stop("nma_fit must be a 'netmeta' object")
  
  z_alpha <- abs(qnorm(alpha / 2))
  z_beta <- abs(qnorm(beta))
  i2 <- nma_fit$I2.random %||% 0
  d2 <- max(0, min(0.99, i2))
  
  variance_base <- 4 
  ris_unadjusted <- variance_base * ((z_alpha + z_beta) / mcid)^2
  daris <- ris_unadjusted / (1 - d2)
  
  current_info_approx <- nma_fit$k 
  target_studies_approx <- max(current_info_approx, round(daris / 100))
  info_fraction <- min(1, current_info_approx / target_studies_approx)
  
  adjusted_z_alpha <- if (info_fraction > 0 && info_fraction < 1) z_alpha / sqrt(info_fraction) else z_alpha
  
  res <- list(daris = daris, info_fraction = info_fraction, adjusted_z_alpha = adjusted_z_alpha, 
              traditional_z_alpha = z_alpha, sufficient_evidence = info_fraction >= 1, mcid = mcid)
  
  class(res) <- "lfahfn_ntsa"
  msg(sprintf("NTSA: Information fraction is %.1f%%. Adjusted Z-boundary: %.2f", info_fraction * 100, adjusted_z_alpha))
  res
}
