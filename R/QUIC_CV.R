# creat folds, from package caret
creatFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE) {
    if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2) 
            cuts <- 2
        if (cuts > 5) 
            cuts <- 5
        breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
        y <- cut(y, breaks, include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            min_reps <- numInClass[i]%/%k
            if (min_reps > 0) {
                spares <- numInClass[i]%%k
                seqVector <- rep(1:k, min_reps)
                if (spares > 0) 
                  seqVector <- c(seqVector, sample(1:k, spares))
                foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
            }
            else {
                foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                  size = numInClass[i])
            }
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
            sep = "")
        if (returnTrain) 
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
}



#TODO: likelihood, cross validation function, push then merge to simulations
# log likelihood for cross validation
#' @name negLLrobOmega
#' @title -log Likelihood on test set
#' @description The default evaluation function in corss validation, -log liekihood on test set
#' @param Omega_hat the estimated *precision* matrix of training set
#' @param Sigma the *covariance* matrix of test sets
#' @return -log likelihood 
#' @export
negLLrobOmega <- function(Omega_hat, Sigma){
    -determinant(Omega_hat)$modulus + sum(diag(Omega_hat %*% Sigma))
}


## several helper functions

# this fit the QUIC on training and test with a given lambda
QUIC_lambda <- function(lambda, trainset, testset,covest,evaluation,...){
    S <- covest(trainset)
    res <- QUIC::QUIC(S, lambda, ...)
    S_test <- covest(testset)
    evaluation(res$X, S_test)
}

# this loop over lambdas in a fold
QUIC_fold <- function(fold, data, covest, lambdas, evaluation, ...){
    trainset <- data[-fold, ]
    testset <- data[fold, ]
    sapply(lambdas, QUIC_lambda, trainset, testset, covest, evaluation, ...)
}

#' @name cvQUIC
#' @title Cross validation to chose tuning parameter of QUIC
#' @description This routine use k fold cross validation to chose tuning parameter
#' @param data The full dataset, should be a matrix
#' @param k number of folds
#' @param covest a *function* that takes a matrix to estimate covariance
#' @param lambdas a vector of tuning parameter to be tested
#' @param evaluation a *function* that takes only two arguments, the estimated *precision* and the test *covariace*, when NULL, we use negative log likelihood on test sets
#' @param ... extra arguments send to QUIC
#' @return a matrix with k rows, each row is the evaluation loss of that fold
#' @examples cvQUIC(matrix(rnorm(100),20,5), msg = 0)
#' @export
cvQUIC <- function(data, k = 10, covest = cov, 
                    lambdas = seq(0.1, 1, 0.1), 
                    evaluation = negLLrobOmega, ...){
    data <- as.matrix(data)
    folds <- creatFolds(1:nrow(data), k = k)
    #browser()
    res <- lapply(folds, QUIC_fold, data, covest, lambdas, evaluation, ...)
    res <- Reduce(rbind, res)
    rownames(res) <- paste0("folds", 1:k)
    return(res)
}


#' @name robQUIC
#' @title QUIC with robust covariance estimations
#' @description This routine fits QUIC using a robust covariance matrix
#' @param data raw data, shoule be a matrix
#' @param covest a *function* that takes a matrix to estimate covariance
#' @param lambda a scalar or vector of tuning parameters, if CV=FALSE, shoule be a scalar, if CV=TRUE scalar input will be override and tuning parameter will be chosed based on CV
#' @param CV bool, whether doing corss validation for tuning parameter, if lambda is a scalar, the candidate will be chosen automatically by log spacing between 0.01 max covariance and max covariance with number of grids
#' @param k fold for corss validation if applicable
#' @param grids number of candidate tuning parameters in cross validation
#' @param evaluation a *function* that takes only two arguments, the estimated *precision* and the test *covariace*, when NULL, we use negative log likelihood on test sets
#' @param ... extra argument sent to QUIC::QUIC
#' @return a QUIC return (see ?QUIC::QUIC), most important one is $X the estimated sparse precision,with an extra entry of tuning parameter lambda
#' @examples robQUIC(matrix(rnorm(100),20,5))
#' @export
robQUIC <- function(data, covest = cov, lambda = 0.1, 
                    CV = FALSE, k = 10, grids = 15, evaluation = negLLrobOmega, ...){
    data <- as.matrix(data)
    S <- covest(data)
    if(length(lambda)!=1 & !CV){
        stop("Provide more than one tuning parameter while not doing cross validation\n")
    }
    if(CV){
        if(length(lambda)<=1){
            lambdas <- seq(log(0.01*max(S[upper.tri(S)])), log(max(S[upper.tri(S)])), length.out = grids)
            lambdas <- exp(lambdas)
        }
        else {
           lambdas <- lambda
        }
        cv_res <- cvQUIC(data, k = k, covest = covest, 
                    lambdas = lambdas, 
                    evaluation = evaluation, ...)
        mean_cv <- colMeans(cv_res)
        lambda <- lambdas[which(mean_cv==min(mean_cv))[1]]
    }
    res <- QUIC::QUIC(S, lambda, ...)
    res$lambda <- lambda
    return(res)
}

