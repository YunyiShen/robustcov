conta_normal_cell_5 <- function(n, Omega,
                        cont_rate = 0.05, 
                        mu = 10, sd = sqrt(0.2)){

    res <- rmvnorm(n, Omega)
    nn <- length(res)
    contind <- sample.int(nn , size = floor(cont_rate * nn))
    res[contind] <- rnorm(length(contind), mu, sd)
    return(res)
}


conta_normal_row_5 <- function(n, Omega,
                        cont_rate = 0.05, 
                        mu = 10, sd = sqrt(0.2)){

    res <- rmvnorm(n, Omega)
    p <- nrow(Omega)
    contind <- sample.int(n , size = floor(cont_rate * n))
    res[contind,] <- rnorm(length(contind)*p, mu, sd)
    return(res)
}

conta_normal_cell_10 <- function(n, Omega,
                        cont_rate = 0.1, 
                        mu = 10, sd = sqrt(0.2)){

    res <- rmvnorm(n, Omega)
    nn <- length(res)
    contind <- sample.int(nn , size = floor(cont_rate * nn))
    res[contind] <- rnorm(length(contind), mu, sd)
    return(res)
}


conta_normal_row_10 <- function(n, Omega,
                        cont_rate = 0.1, 
                        mu = 10, sd = sqrt(0.2)){

    res <- rmvnorm(n, Omega)
    p <- nrow(Omega)
    contind <- sample.int(n , size = floor(cont_rate * n))
    res[contind,] <- rnorm(length(contind)*p, mu, sd)
    return(res)
}
