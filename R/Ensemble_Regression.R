Ensemble_Function_Continuous <- function(x,y){
  x <- as.matrix(x[!is.na(y),])
  y <- y[!is.na(y)]

  lasso_train <- glmnet::glmnet(x,y,family = "gaussian",alpha = 1)
  ridge_train <- glmnet::glmnet(x,y,family = "gaussian",alpha = 0)

  lasso_prs_tune <- predict(lasso_train,x)
  ridge_prs_tune <- predict(ridge_train,x)

  all <- cbind(lasso_prs_tune,ridge_prs_tune)

  R2_Vector <- vector()
  for(i in 1:ncol(all)){
    tmp <- data.frame(y = y, x_try = all[,i])
    R2_Vector[i] <- summary(lm(y~x_try,data = tmp))$r.square
  }

  coefficients_x <- coef(lm(y~.,data.frame(y = all[,which.max(R2_Vector)],x)))
  return(list(Coefficients = coefficients_x))
}
Ensemble_Function_Binary <- function(x,y){
  x <- as.matrix(x[!is.na(y),])
  y <- y[!is.na(y)]

  lasso_train <- glmnet::glmnet(x,y,family = "binomial",alpha = 1)
  ridge_train <- glmnet::glmnet(x,y,family = "binomial",alpha = 0)

  lasso_prs_tune <- predict(lasso_train,x)
  ridge_prs_tune <- predict(ridge_train,x)

  all <- cbind(lasso_prs_tune,ridge_prs_tune)

  AUC_Vector <- vector()
  for(i in 1:ncol(all)){
    tmp <- data.frame(y = y, x_try = all[,i])
    roc_obj <- RISCA::roc.binary(status = "y",
                          variable = "x_try",
                          confounders = "~1",
                          data = tmp,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_Vector[i] <- roc_obj$auc
  }

  coefficients_x <- coef(lm(y~.,data.frame(y = all[,which.max(AUC_Vector)],x)))
  return(list(Coefficients = coefficients_x))
}

#' @title Perform Rare Variant PRS Ensemble Regression
#' @description Performs ensemble regression of PRSs and continuous or binary response Y. Produces a set of coefficients for the PRSs.
#'
#' @param PRSs A data.frame or matrix containing p PRSs and n observations. p must be greater than 1 and n must be equal the number of observations in Y. PRSs and Y must be matched before performing ensemble regression.
#' @param Y A data.frame, matrix, or vector containing continuous or binary outcome Y. If a data.frame or matrix, first column will be used as response. Number of observations must equal the number of observations in PRSs. PRSs and Y must be matched before performing ensemble regression.
#' @param family Either continuous or binary depending on outcome Y. Defaults to continuous.
#' @returns A named list with one item, Coefficients. Coefficients contain coefficients of the PRSs columns to generate the predicted ensemble PRS. Coefficients include the intercept estimate.
#' @export
Ensemble_Function <- function(PRSs,Y,family = c("continuous","binary")){

  family <- match.arg(family)

  if(!is.matrix(PRSs) & !is.data.frame(PRSs)){
    stop("PRSs is not a data.frame or a matrix")
  }

  if(dim(PRSs)[2] < 2){
    stop("Number of columns of PRSs must be greater than 1")
  }

  if(!is.matrix(Y) & !is.data.frame(Y) & !is.vector(Y)){
    stop("Y is not a data.frame, a matrix, or a vector")
  }

  if(is.matrix(Y) | is.data.frame(Y)){
    Y <- Y[,1]
  }

  if(dim(PRSs)[1] != length(Y)){
    stop("Number of rows of PRSs must match the dimension of Y")
  }

  if(family == "continuous"){
    return(Ensemble_Function_Continuous(PRSs,Y))
  }else{
    return(Ensemble_Function_Binary(PRSs,Y))
  }
}
