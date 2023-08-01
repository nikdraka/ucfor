#' Forecast reconciliation
#'
#' @param forecast A matrix with a dimension of p (row) and h (column), where p is the number of time series and h is the forecast horizon
#' @param residual A matrix with a dimension of p (row) and T (column), where p is the number of time series and T is the number of observations
#' @param weights A reconciliation weight matrix.
#' @param linearComb A matrix of a linear combination. In hierarchical forecasting, this is the summation matrix or S.
#'
#' @return
#' Function returns the following variables:
#' \itemize{
#' \item{\code{reconciledForecast} - the reconciled forecasts}
#' \item{\code{baseForecast} - the base forecasts}
#' \item{\code{reconMatrix} - the reconciliation weights matrix}
#' \item{\code{smatrix} - the linear combination. Usually, it is the summation matrix}
#' \item{\code{invcovMatrix} - the inverse of the covariance matrix estimation approximation}
#' }
#'
#' @examples
#'
#' @importFrom MASS ginv
#' @importFrom corpcor invcov.shrink
#'
#' @export combForecast
combForecast <- function(forecast, residual, weights = c("wls", "cwls", "bshr", "pshr", "mint"), linearComb = S) {

	# p: number of series
	# t: number of obs
	# forecast: p x t
	# residuals: p x t

	S <- linearComb
	nSeries <- nrow(forecast)
	nObs <- ncol(forecast)
	nBottom <- ncol(S)

	if (!is.matrix(forecast) || !is.matrix(residual) || !is.matrix(S)) {
		stop("forecast and residual should be a matrix!")
	}

	if (nrow(S) != nrow(forecast)) {
		stop("Change the dimension of forecast. The number of rows is the number of time series in the hierarchy.")
	}

	if (nrow(S) != nrow(residual)) {
		stop("Change the dimension of residual. The number of rows is the number of time series in the hierarchy.")
	}

	yhat <- forecast
	residual <- residual

	if (weights == "ols" || is.null(weights)) {
		invSigma <- diag(nSeries)
	} else if (weights == "scl") {
		invSigma <- diag(1/rowSums(S))
	} else if (weights == "wls") {
		invSigma <- diag(1/rowMeans(residual^2))
	} else if (weights == "cwls") {
		invSigma <- diag(1/c(S %*% rowMeans(residual[nBottom:nSeries,]^2)))
	} else if (weights == "bshr") {
		mintSigma <- invcov.shrink(t(residual), verbose = FALSE)
		invSigma <- ginv(S %*% mintSigma[nBottom:nSeries, nBottom:nSeries] %*% t(S))
	} else if (weights == "mint") {
		invSigma <- invcov.shrink(t(residual), verbose = FALSE)[1:nSeries, 1:nSeries]
	} else if (weights == "pshr") {
		mintSigma <- invcov.shrink(t(residual), verbose = FALSE)
		mintSigma[1:(nBottom-1), nBottom:nSeries] <- 0
		mintSigma[nBottom:nSeries, 1:(nBottom-1)] <- 0
		invSigma <- mintSigma
	} else {
		invSigma <- diag(nSeries)
		warning("Weights is undefined. Return to 'ols'!")
	}

	G <- ginv(t(S) %*% invSigma %*% S) %*% t(S) %*% invSigma

	# Reconciliation
	ytilde <- S %*% G %*% yhat

	if (is.null(rownames(ytilde))) {
		rownames(ytilde) <- paste("Series", 1:nSeries)
	}

	listReturned <- list(reconciledForecast = ytilde,
											 baseForecast = forecast,
											 reconMatrix = G,
											 smatrix = S,
											 invcovMatrix = invSigma)

	return(structure(listReturned,class="reconciliation"))

}

