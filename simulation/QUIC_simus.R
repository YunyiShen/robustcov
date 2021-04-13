#devtools::install_github("cran/QUIC")
#devtools::install_github("YunyiShen/RobustOmega")
#install.packages("flare")
library(flare) # for Tiger
library(QUIC) # for Glasso
library(RobustOmega) # robust covariances implemented and wrapper for QUIC
# for the performace evaluation function
source("./R/misc.R")
source("./R/graph_generator.R")
source("./simulation/cont_norm.R") # contaimanted norm

results_file <- "./simulation/Res/res1.csv"

graphs <- c("omega_sparse", "omega_banded", "omega_dense", "diag")
covmats <- c(cov,corSpearmanmat,corKendallmat, 
            #corQuadrantmat,covGKmat, 
            #covNPDmat, 
            covOGKmat, covSpearmanUmat)
covmats_name <- c("cov","corSpearmanmat","corKendallmat", 
            #"corQuadrantmat","covGKmat", 
            #"covNPDmat", 
            "covOGKmat", "covSpearmanUmat")
distributions <- c("rmvnorm","conta_normal_cell_5","conta_normal_cell_10", 
                    "conta_normal_row_5","conta_normal_row_10", "rmvt", 
                    "raltert")

ps <- c(120, 400)
n <- 200
rep <- 100

tot_num <- rep * length(ps) * length(graphs) * 
    length(covmats) * length(distributions)

res <- data.frame(matrix(NA, tot_num, 9))
colnames(res) <- c("cov_measure","p",
                    "graph", "distribution", "rep", 
                    "norm_Omega", "norm_Sigma",
                    "FP", "FN")

i_data <- 1
# since the algorithms are already parallelized, will do for here in R

for(p in ps){
    for(g in graphs){
        for(d in distributions){
            for(i in 1:length(covmats)){
                cat("covmat:" ,covmats_name[i], " p =",p," g =", g, " d =", d, " rep = CV", "\n")
                # will do CV here only ones
                if(g=="diag") omega <- diag(p)
                else{
                    omega <- do.call(g, list(p=p))
                }
                dist_para <- list(n = n, Omega = omega)
                if(d %in% c("rmvt","raltert")) dist_para$nu = 3
                data1 <- do.call(d, dist_para)
                QUICres <- tryCatch( robQUIC(data=data1, covest = covmats[[i]], CV = TRUE, msg = 0), 
                                    error =function(e){list()})
                if(length(QUICres)>0){
                    omega_hat <- QUICres$X
                    lambda <- QUICres$lambda # get the tunning
                    losses <- perform_metrics(omega_hat, omega)
                }
                else {
                   lambda = 0.1
                   losses <- c(NA,NA,NA,NA)
                }
                

                res$cov_measure[i_data] <- covmats_name[i]
                res$p[i_data] <- p
                res$graph[i_data] <- g
                res$distribution[i_data] <- d
                res$rep[i_data] <- 1
                res[i_data, 6:9] <- losses
                i_data <- i_data + 1
                write.csv(res,results_file, row.names = FALSE)
                for(j in 2:rep){
                    cat("covmat:" ,covmats_name[i], " p =",p," g =", g, " d =", d, " rep =",j, "\n")
                    data1 <- do.call(d, dist_para)
                    QUICres <- tryCatch( robQUIC(data=data1, covest = covmats[[i]], lambda = lambda,CV = FALSE, msg = 0), 
                                    error =function(e){list()})
                    if(length(QUICres)>0){
                        omega_hat <- QUICres$X
                        losses <- perform_metrics(omega_hat, omega)
                    }
                    else {
                        losses <- c(NA,NA,NA,NA)
                    }
                    res$cov_measure[i_data] <- covmats_name[i]
                    res$p[i_data] <- p
                    res$graph[i_data] <- g
                    res$distribution[i_data] <- d
                    res$rep[i_data] <- j
                    res[i_data, 6:9] <- losses
                    i_data <- i_data + 1
                    write.csv(res,results_file, row.names = FALSE)
                }
            }
        }
    }
}
