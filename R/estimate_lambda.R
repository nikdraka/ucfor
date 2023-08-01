#' Estimate a hyper-parameter in the shrinkage estimators for exponential smoothing models
#'
#' @param object An adam object
#' @param lambda An initialisation for lambda. Should be between zero and one.
#' @param origins A number of forecast origins. The default is the frequency of the time series.
#' @param ... Other parameters passed to adam.
#'
#' @return Function returns the following variables:
#' \itemize{
#' \item{\code{lambda_{min}} - the optimal hyper-parameter of ETS}
#' \item{\code{model} - the model}
#' \item{\code{loss} - the loss function, either "RIDGE" or "LASSO"}
#' \item{\code{data} - the training set}
#' }
#'
#' @author Kandrika Pritularga
#'
#' @references \itemize{
#' \item Pritularga, K., Svetunkov, I., Kourentzes, N., (2022) Shrinkage estimator for exponential smoothing models. International Journal of Forecasting, XX, pp. XX.
#' }
#'
#' @examples
#'
#' library(smooth)
#' library(greybox)
#'
#' # Estimate an exponential smoothing model with RIDGE first
#' fit <- adam(AirPassengers, loss = "RIDGE", lambda = 0.5)
#'
#' # Implement estimate_lambda() to find the optimal hyper-parameter
#' choosingLambda <- estimate_lambda(fit, lambda = 0.5, origins = 5)
#'
#' # Inspect the result
#' choosingLambda
#'
#' # Re-estimate adam() with the result of estimate_lambda()
#' fitShrinkage <- adam(choosingLambda$data, model = choosingLambda$model,
#'                      loss = choosingLambda$loss,
#'                      lambda = choosingLambda$lambda_min)
#'
#' # Inspect the effect of smoothing parameter shrinkage
#' fit$persistence
#' fitShrinkage$persistence
#'
#' @importFrom smooth is.adam
#' @importFrom stats frequency
#'
#' @export estimate_lambda
estimate_lambda <- function(object, lambda = 0.1, origins = 5, ...) {
  # I removed loss because it makes sense to align it with the loss used in the initial object.
  # I renamed x0 to lambda for consistency with adam()

  # Input checking
  if (!smooth::is.adam(object)) {
    stop("An adam object is needed to run this function. Use adam() as the object!")
  }

  if (lambda<0 || lambda >= 1 || !is.numeric(lambda)) {
    lambda <- 0.1
  } else {
    lambda <- lambda[1]
  }

  if (any(loss != c("RIDGE", "LASSO"))) {
    loss <- "RIDGE"
  } else {
    loss <- loss[1]
  }

  if (origins < 0 || is.integer(origins) || is.null(origins)) {
    origins <- stats::frequency(object$data)
  } else {
    origins <- origins[1]
  }

  # collect arguments from 'object'
  modelNames <- modelType(object)
  data <- object$data

  # A function to calculate error
  roLambda <- function(data = data, model = modelNames, loss = loss, lambda = x0, origins = origins) {
# The stuff that was here is not needed, you already define model, loss etc in the call of the function

    stringLambda <- paste0("lambda=",lambda)
    stringModel <- paste0("model=","'",as.character(model),"'")
    stringLoss <- paste0("loss=","'",as.character(loss), "'")

    ourCall <- paste("adam(data", stringModel, stringLoss, stringLambda, "h = 1,holdout=TRUE)", sep=",")
    ro.fit <- greybox::ro(data, h = 1, origins = origins, ourCall, "forecast")

    ro.error <- ro.fit$holdout - ro.fit$forecast
    yDenominator <- mean(abs(diff(ro.fit$actuals)))
    scaled.ro.error <- ro.error/yDenominator

    return(mean((scaled.ro.error)^2))

  }

  # optimisation
  opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel" = 1e-08, "maxeval" = 1000)
  lb <- 0
  ub <- 0.9999

  regCF <- nloptr::nloptr(x0, function(x) roLambda(data=data, model = modelNames, loss = loss, lambda = x, origins = origins),
                          lb = lb, ub = ub, opts = opts)

  listReturned <- list(lambda_min = regCF$solution,
                       model = modelNames,
                       loss = loss,
                       data = data)

  return(structure(listReturned,class="shrink"))

}

#' estimate_lambda.default <- function(object,x0 = 0.1,
#'                                     loss=c("RIDGE","LASSO"),
#'                                     origins = 5,...){
#'   return("The default method is not available")
#' }
