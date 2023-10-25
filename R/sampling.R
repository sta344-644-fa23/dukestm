#' @title Sample from MVN normal
#'
#' @description Generate draws of the given MVN distribution.
#'
#' @returns A matrix of size nrow(Sigma) x n
#'
#' @param n Number of draws to sample
#' @param mu Mean vector
#' @param Sigma Covariance matrix
#' @param diag_adj Numeric. Small value that is added covariance's diagonal to avoid
#' numerical singularity issues.
#'
#' @export
#'
rmvnorm = function(n, mu=rep(0, nrow(Sigma)), Sigma, diag_adj = 1e-6) {
  diag(Sigma) = diag(Sigma) + diag_adj
  mu %*% matrix(1, ncol=n) + t(chol(Sigma)) %*% matrix(stats::rnorm(n*nrow(Sigma)), ncol=n)
}

#' @title GP conditional predictions
#'
#' @description Given a covariance function and parameters this samples from
#' the conditional MVN of y_pred | y
#'
#' @returns A matrix of size `n_pred` x `reps`.
#'
#' @param x_pred Vector or matrix of prediction locations (values).
#' @param x Vector or matrix of fitted locations (values).
#' @param mu Vector of mean values for the observed locations.
#' @param mu_pred Vector of mean values for the prediction locations.
#' @param cov Function. Covariance function that takes a distance matrix as input.
#' @param ... Parameter values used by `cov`.
#' @param reps Numeric. Number of samples to draw from conditional distribution.
#' @param diag_adj Numeric. Small value that is added covariance diagonal to avoid
#' numerical singularity issues.
#' @param cov_f_o Function. Calculates the covariance matrix for the observed locations.
#' @param cov_f_p Function. Calculates the covariance matrix for the prediction locations.
#' @param cov_f_o Function. Calculates the cross covariance matrix between the prediction and observed locations.
#'
#' @export
#'
cond_predict = function(
    y, x, x_pred,
    mu=rep(0,length(x)), mu_pred=rep(0,length(x_pred)),
    cov, ..., reps=1000, diag_adj = 1e-6,
    cov_f_o = cov, cov_f_p = cov, cov_f_po = cov
) {
  y = as.matrix(y)
  x = as.matrix(x)
  x_pred = as.matrix(x_pred)

  dist_o  = fields::rdist(x)
  dist_p  = fields::rdist(x_pred)
  dist_po = fields::rdist(x_pred, x)

  cov_o  = cov_f_o(dist_o, ...)
  cov_p  = cov_f_p(dist_p, ...)
  cov_po = cov_f_po(dist_po, ...)

  # Quick fix for singularity issues
  diag(cov_o) = diag(cov_o) + diag_adj
  diag(cov_p) = diag(cov_p) + diag_adj

  cond_cov = cov_p - cov_po %*% solve(cov_o, t(cov_po))
  cond_mu  = mu_pred + cov_po %*% solve(cov_o) %*% (y-mu)

  cond_mu %*% matrix(1, ncol=reps) + t(chol(cond_cov)) %*% matrix(stats::rnorm(nrow(x_pred)*reps), ncol=reps)
}
