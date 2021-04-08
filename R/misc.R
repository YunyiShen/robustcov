perform_metrics <- function(graph_hat,graph, tol = 1e-8){
    normOmega <- norm(graph_hat-graph, type = "I")
    normSigma <- norm(solve(graph_hat)-solve(graph), type = "I")
    graph <- graph != 0
    graph_hat <- abs(graph_hat) <= tol
    FP <- sum(graph_hat[upper.tri(graph_hat)] & !graph[upper.tri(!graph)])
    FN <- sum(!graph_hat[upper.tri(!graph_hat)] & graph[upper.tri(graph)])
    return(c(normOmega,normSigma,FP,FN))
}
