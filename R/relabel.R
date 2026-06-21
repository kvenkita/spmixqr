## Label handling. order_components is vendored from the mixqr/mixqrgate internal
## (not exported). label_stability is new: it measures whether a single global label
## order is coherent across space, which a roughness penalty does NOT guarantee
## (Reviewer A; Green & Richardson 2002 note coherence needs an allocation prior,
## deferred to v2). The diagnostic detects crossing slope surfaces and warns.

#' Global component order from a constant-coefficient matrix.
#' @param beta_const `p x G` constant coefficients (intercept + constant slopes).
#' @param label_order `"slope"` or `"intercept"`.
#' @param order_var Slope index used when ordering by slope (1 = first slope).
#' @return An integer permutation (ascending key).
#' @keywords internal
order_components <- function(beta_const, label_order = "slope", order_var = 1L) {
  if (ncol(beta_const) < 2L) return(seq_len(ncol(beta_const)))
  key <- if (label_order == "intercept" || nrow(beta_const) < 2L)
    beta_const[1L, ] else beta_const[1L + order_var, ]
  order(key)
}

#' Spatial label-coherence diagnostic.
#'
#' Given each regime's ordering-covariate slope evaluated at every location
#' (`slope_surface`, an `n x G` matrix), returns the fraction of locations where the
#' pointwise ranking of the regimes disagrees with the global (mean-slope) ranking. A
#' high value means slope surfaces cross in space and a single global label is not
#' pointwise coherent (v1 limitation; per-region alignment is v2).
#' @keywords internal
label_stability <- function(slope_surface) {
  if (ncol(slope_surface) < 2L) return(0)
  global <- order(colMeans(slope_surface))
  disagree <- apply(slope_surface, 1L,
                    function(v) !identical(order(v), global))
  mean(disagree)
}
