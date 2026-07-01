## Observation (case/sampling) weights layer. Distinct from the SPATIAL weights
## matrix W (see weights.R / spq_weights): these weight each observation's
## contribution to the penalised EM, the log-likelihood, and inference.

#' Resolve, validate and normalize observation weights.
#'
#' @param weights A length-n numeric vector, a column name in `data`, a one-sided
#'   formula `~w`, or `NULL` (unweighted).
#' @param weights_type `"sampling"`, `"frequency"`, or `"precision"`.
#' @param data The model data frame (for column-name / formula resolution).
#' @param n Number of observations.
#' @return list(w, raw, type, sum_raw, weighted). `w` is normalized to mean 1
#'   (Σ = n); `raw` is the pre-normalization vector (or `NULL` when unweighted).
#' @keywords internal
resolve_weights <- function(weights, weights_type = c("sampling", "frequency", "precision"),
                            data = NULL, n) {
  weights_type <- match.arg(weights_type)
  if (is.null(weights))
    return(list(w = rep(1, n), raw = NULL, type = weights_type,
                sum_raw = as.numeric(n), weighted = FALSE))
  raw <- if (inherits(weights, "formula")) {
    v <- all.vars(weights)
    if (length(v) != 1L) stop("`weights` formula must name exactly one variable.", call. = FALSE)
    eval(as.name(v), envir = data)
  } else if (is.character(weights) && length(weights) == 1L) {
    if (is.null(data) || !weights %in% names(data))
      stop("`weights` column '", weights, "' not found in `data`.", call. = FALSE)
    data[[weights]]
  } else weights
  raw <- as.numeric(raw)
  if (length(raw) != n)
    stop("`weights` has length ", length(raw), " but the model has ", n,
         " observations (length mismatch).", call. = FALSE)
  if (anyNA(raw) || any(!is.finite(raw)))
    stop("`weights` must be finite (no NA/Inf).", call. = FALSE)
  if (any(raw < 0)) stop("`weights` must be non-negative.", call. = FALSE)
  s <- sum(raw)
  if (s <= 0) stop("`weights` must not be all zero.", call. = FALSE)
  list(w = raw * n / s, raw = raw, type = weights_type,
       sum_raw = s, weighted = TRUE)
}
