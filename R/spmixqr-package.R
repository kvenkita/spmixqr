#' spmixqr: Spatial Finite Mixtures of Quantile Regressions
#'
#' Fits spatial finite mixtures of quantile regressions. The mixing (gating)
#' probabilities vary over space through a low-rank spatial basis, and each
#' regime's covariate-effect slopes are spatially varying surfaces (with a scalar
#' per-regime intercept). Estimation is a penalised expectation-maximisation
#' algorithm that reuses the \pkg{mixqr} component machinery and a vendored,
#' spatially penalised multinomial-logit gate. The single-regime case is a
#' penalised spatially-varying-coefficient quantile regression in the lineage of
#' Reich, Fuentes and Dunson (2011).
#'
#' @section Entry points:
#' \itemize{
#'   \item [spmixqr()] --- fit a spatial mixture of quantile regressions.
#'   \item [spmixqr_select()] --- choose the number of regimes and/or the
#'         smoothing parameters.
#'   \item [sim_spmixqr()] --- simulate from the model (for validation).
#'   \item [coef_surface()], [gate_surface()] --- spatial-surface accessors.
#' }
#'
#' @references
#' Reich, B. J., Fuentes, M. and Dunson, D. B. (2011). Bayesian spatial quantile
#' regression. \emph{Journal of the American Statistical Association} 106, 6--20.
#'
#' Wu, C. and Yao, W. (2016). Mixtures of quantile regressions.
#' \emph{Computational Statistics & Data Analysis} 93, 162--176.
#'
#' Fernandes, M., Guerre, E. and Horta, E. (2021). Smoothing quantile
#' regressions. \emph{Journal of Business & Economic Statistics} 39, 338--357.
#'
#' @keywords internal
#' @importFrom stats AIC BIC logLik coef confint vcov fitted residuals nobs predict
#' @importFrom utils globalVariables
"_PACKAGE"

## Quiet the no-visible-binding NOTE for the mgcv smooth-constructor column names
## referenced via non-standard evaluation inside s().
globalVariables(c(".s1", ".s2", ".region"))
