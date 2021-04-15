perform_metrics <- function(graph_hat,graph, tol = 1e-10){
    normOmega <- norm(graph_hat-graph, type = "I")
    normSigma <- norm(solve(graph_hat)-solve(graph), type = "I")
    graph <- graph != 0
    graph_hat <- abs(graph_hat) >= tol
    FP <- sum(graph_hat & !graph)/sum(!graph)
    FN <- sum(!graph_hat & graph)/sum(graph)
    return(c(normOmega,normSigma,FP,FN))
}
