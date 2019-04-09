library(igraph)

evaluate.components <- function(rho, diss.matrix) {
  print(paste("Testing rho value: ", as.character(rho)))
  # filter the adjacency matrix and create a new graph
  gene.network <- replace(diss.matrix, abs(1-diss.matrix) <= rho, 0)
  gene.graph <- graph_from_adjacency_matrix(gene.network, weighted = T)
  
  # generate degree distributions
  distrib.deg <- degree_distribution(gene.graph)
  degree.dist <- data.frame(freq=distrib.deg, x=0:(length(distrib.deg)-1))
  degree.dist <- subset(degree.dist, freq != 0)
  
  # non-igraph alternative
  # library(plyr)
  # degree.dist <- count(rowSums(diss.matrix != 0))
  # degree.dist$freq <- degree.dist$freq/sum(degree.dist$freq)
  
  # linear model
  fitting <- lm(formula=log10(freq) ~ I(log10(x+1)), data=degree.dist)
  #  plots for each linear model (at growing rho)
  plot(log10(freq) ~ I(log10(x+1)), data=degree.dist,
       main=paste("rho:", rho, ", gamma:",round(fitting$coefficients[2], 2)))
  abline(fitting$coefficients, col="red", lwd=2)
  
  # signed R-squared
  r.squared <- ifelse(fitting$coefficients[2] <= 0, summary(fitting)$r.squared, 
                      -summary(fitting)$r.squared)
  
  # returning biggest component fraction and signed R^2
  return(c(max(components(gene.graph)$csize)/length(V(gene.graph)), r.squared))
}

connectivity.check <- function(diss.matrix, min.rho, max.rho, step.rho, actual.rho, parallel.exec=F) {
  # checks the connectivy of the network at varying pearson correlation coefficient threshold
  rhos <- seq(min.rho, max.rho, step.rho)
  par(mfrow=c(5,4), mar=c(2,2,2,2), cex=0.6)
  if (!parallel.exec) {
    evaluation <- lapply(rhos, evaluate.components, diss.matrix=diss.matrix)
    big.comp.frac <- unlist(lapply(evaluation, function(x) x[1]))
    scalefreeness <- unlist(lapply(evaluation, function(x) x[2]))
  } else {
    evaluation <- unlist(mclapply(rhos, evaluate.components, diss.matrix=diss.matrix,
                                  mc.cores = detectCores()-1))
    big.comp.frac <- unlist(lapply(evaluation, function(x) x[1]))
    scalefreeness <- unlist(lapply(evaluation, function(x) x[2]))
  }
  
  # eventually plots everything with a line to indicate the actual rho used
  # par(mfrow=c(1,2))
  plot(big.comp.frac~rhos, type='l', pch=19, col="blue",
       main="Biggest connected component", ylim=c(0, 1.2))
  abline(v=actual.rho, col="red")
  axis(1, at=seq(min.rho, max.rho, 0.1))
  plot(scalefreeness~rhos, type='l', pch=19, col="red",
       main="R^2 of fit to scale-free model",ylim=c(-0.7,1),xlim=c(0,1))
  abline(h=0.8, col="green" )
  return(evaluation)
}

# this is the code to run to get the same image I sent you by email
# diss.matrix <- 1 - t(cor(filtered.genes))
evaluation <- connectivity.check(diss.matrix, 0.1, 0.95, 0.05, 0.62, parralel.exec=F)