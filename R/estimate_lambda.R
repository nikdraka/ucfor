#' Estimate a hyper-parameter in the shrinkage estimators for exponential smoothing models
#'
#' @param object An adam object
#' @param lambda An initialisation for lambda. Should be between zero and one.
#' @param origins A number of forecast origins. The default is the frequency of the time series.
#' @param ci The parameter defines if the in-sample window size should be constant. If TRUE,
#' then with each origin one observation is added at the end of series and another one
#' is removed from the beginning.
#' @param co The parameter defines whether the holdout sample window size should be constant.
#' If TRUE, the rolling origin will stop when less than h observations are left in the holdout.
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
#' @seealso \code{\link[smooth]{adam}, \link[greybox]{ro}}

#' @author Kandrika Pritularga
#'
#' @references \itemize{
#' \item Pritularga, K., Svetunkov, I., Kourentzes, N., (2022) Shrinkage estimator for exponential smoothing models. International Journal of Forecasting, XX, pp. XX.
#' }
#'
#' @examples
#'
#' library(smooth)
#'
#' # Estimate an exponential smoothing model with RIDGE
#' fit <- adam(AirPassengers, model="MAM", loss = "RIDGE", lambda = 0.5)
#'
#' # Implement estimate_lambda() to find the optimal hyper-parameter
#' choosingLambda <- estimate_lambda(fit, lambda = 0.5, origins = 5)
#'
#' # Inspect the result
#' choosingLambda
#'
#' # Re-estimate adam() with the result of estimate_lambda()
#' fitShrinkage <- adam(AirPassengers, model = modelType(choosingLambda$model),
#'                      loss = choosingLambda$model$loss,
#'                      lambda = choosingLambda$lambda_min)
#'
#' # Inspect the effect of smoothing parameter shrinkage
#' fit$persistence
#' fitShrinkage$persistence
#'
#' @importFrom smooth is.adam modelType
#' @importFrom greybox ro
#' @importFrom stats frequency
#' @importFrom nloptr nloptr
#'
#' @export estimate_lambda
estimate_lambda <- function(object, lambda = 0.1, origins = 5,
                            ci=FALSE, co=TRUE, ...) {
  #!!! I removed loss because it makes sense to align it with the loss used in the initial object.
  #!!! I renamed x0 to lambda for consistency with adam()

  # Input checking
	#!!! Also, no need to do smooth::function - use @importFrom to use specific function, drop ::
  if (!is.adam(object)) {
    stop("An adam object is needed to run this function. Use adam() as the object!")
  }

  if (lambda<0 || lambda >= 1 || !is.numeric(lambda)) {
    warning("The initial lambda is not between 0 and 1. Setting it to 0.1.")
    lambda <- 0.1
  }

  # Grab loss from adam()
  loss <- object$loss
  if(all(loss!=c("LASSO","RIDGE"))){
    stop("adam() needs to be estimated with LASSO/RIDGE in order for this function to work")
  }

  if (origins < 0 || !is.numeric(origins) || is.null(origins)) {
    warning("The number of origins is not a positive number. Setting it to 5.")
    #!!! Why frequency?
    # origins <- stats::frequency(object$data)
    origins <- 5
  }

  # collect arguments from 'object'
  #!!! modelName - singular, right?
  modelName <- modelType(object)
  data <- object$data

  # A function to calculate error
  roLambda <- function(lambda = lambda, data = data, model = modelName, loss = loss, origins = origins) {
#!!! The stuff that was here is not needed, you already define model, loss etc in the call of the function

    stringLambda <- paste0("lambda=",lambda)
    stringModel <- paste0("model=","'",as.character(model),"'")
    stringLoss <- paste0("loss=","'",as.character(loss), "'")

    ourCall <- paste("adam(data", stringModel, stringLoss, stringLambda, "h = 1,holdout=TRUE,...)", sep=",")
    ro.fit <- ro(data, h = 1, origins = origins, call=ourCall, value="forecast")

    ro.error <- ro.fit$holdout - ro.fit$forecast
    yDenominator <- mean(abs(diff(ro.fit$actuals)))
    scaled.ro.error <- ro.error/yDenominator

    return(mean((scaled.ro.error)^2))

  }

  # optimisation
  opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel" = 1e-08, "maxeval" = 1000)
  lb <- 0
  ub <- 0.9999

  regCF <- nloptr(lambda, roLambda, lb = lb, ub = ub, opts = opts,
  								data=data, model = modelName, loss = loss, origins = origins,)

  objectUpdated <- adam(data, model = modelName, loss=loss, lambda = regCF$solution)

  # loss and data are now saved in model (which is now the provided object)
  listReturned <- list(lambda_min = regCF$solution,
                       model = objectUpdated)

  return(structure(listReturned,class="shrink"))

}

# A function to print outputs of the shrink class
# print.shrink <- function(x, ...){
# }

#' estimate_lambda.default <- function(object,lambda = 0.1,
#'                                     loss=c("RIDGE","LASSO"),
#'                                     origins = 5,...){
#'   return("The default method is not available")
#' }
