################ 1 - PREPROCESSING ########

clean.zeros <- function(dataset, threshold) {
  # if the number of zeros in a gene is higher than the threshold, the gene is removed from the dataset
  # assuming dataset is shaped like (genes, samples) 
  return(dataset[rowSums(dataset == 0) <= threshold*dim(dataset)[2],])
}

clean.iqr <- function(dataset, threshold) {
  # the gene is removed from the analysis if the IQR size does not pass the threshold level
  # assuming the dataset is shaped like (gene, samples)
  iqrs <- apply(dataset, 1, IQR)
  return(dataset[iqrs >= quantile(iqrs, threshold),])
}

draw.iqr.dist <- function(dataset, threshold) {
  # IQR for every gene
  iqrs <- apply(dataset, 1, IQR)
  hist(iqrs,col = "mistyrose", main = "IQRs distribution", breaks = 80)
  abline(v=quantile(iqrs, threshold), lty=2, lwd=2, col="red")
}

draw.iqr.zeros <- function(dataset, threshold.iqr, threshold.zeros) {
  # % of zero plotted against IQR values
  iqrs <- apply(dataset, 1, IQR)
  zeros <- apply(dataset == 0, 1, sum)* 100/dim(dataset)[2]
  plot(iqrs, 100 - zeros, xlab = "IQR", ylab = "100 - % of zeros", cex=0.6)
  abline(v=quantile(iqrs, threshold.iqr), h=100-threshold.zeros*100, lty=2, col="red")
}

################## 2 - FILTERING #########

t.tester <- function(dataset, A, B) {
  # applies a t-test to every gene to check if the mean expression between Tumor and Control is significantly different
  t.results <- apply(dataset, 1, function(x) {t.test(x[A[[1]]], x[B[[1]]], var.equal = T)})
  print("T-tests completed")
  # creating a dataframe that returns all necessary data to keep going
  p.value <- unlist(lapply(t.results, function(x) {x$p.value}))
  t.counts<- data.frame(p.value=p.value, fdr= p.adjust(p.value, method = 'fdr'))
  print("calculating FDR and log-fold change...")
  
  #Log-fold change
  means.condA <- apply(dataset[A[[1]]], 1, mean)
  means.condB <- apply(dataset[B[[1]]], 1, mean)
  t.counts$lfc <- means.condA - means.condB
  
  return(list('matrix'=dataset, 'stats'=t.counts))
}

t.filter <- function(dataset, lfc, FDR) {
  # filters genes according to some log-fold change and adjusted p-value thresholds
  t.dataset <- dataset$stats
  idx <- t.dataset$fdr <= FDR & abs(t.dataset$lfc) >= log2(lfc)
  dataset$matrix <- dataset$matrix[idx,]
  dataset$stats <- dataset$stats[idx,]
  return(dataset)
}

draw.volcanoplot <- function (dataset, lfcthresh=2, fdrthresh=0.05, main="Volcano Plot", legendpos="topright", cex=0.7, ...) {
  require(calibrate)
  # for values less than 1 it'll be negative, but it's a meaningless measure
  lfcthresh <- log2(lfcthresh)
  with(dataset, plot(lfc, -log10(fdr), pch=20, main=main, ...))
  with(subset(dataset, fdr<fdrthresh ), points(lfc, -log10(fdr), pch=20, col="red", ...))
  with(subset(dataset, abs(lfc)>lfcthresh), points(lfc, -log10(fdr), pch=20, col="orange", ...))
  with(subset(dataset, fdr<fdrthresh & abs(lfc)>lfcthresh), points(lfc, -log10(fdr), 
                                                                   pch=20, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",fdrthresh ,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"),
         pch=20, col=c("red","orange","green"), text.width = strwidth("|LogFC|>")[1]/8)
}

draw.lfc.hist <- function(dataset, lfcthresh) {
  lfcthresh <- log2(lfcthresh)
  base.hist <- hist(dataset$lfc, breaks = 300)
  hist.colors <- ifelse(base.hist$breaks <= -lfcthresh |
                          base.hist$breaks >= lfcthresh, "red", "gray")
  plot(base.hist, col=hist.colors, sub = "Distribution of log-fold change")
}

# CAMBIA COLORI
biclustering <- function(dataset) {
  heatmap.matrix <- as.matrix(dataset)
  mode(heatmap.matrix) <- "numeric"
  
  is.na(heatmap.matrix) <- sapply(heatmap.matrix, is.infinite)
  heatmap.matrix[is.na(heatmap.matrix)] <- 0
  heatmap.matrix[is.nan(heatmap.matrix)] <- 0
  heatmap.2(heatmap.matrix,
            trace = 'none',
            scale = 'row',
            col = colorRampPalette(c("blue", "blue3", "black", "yellow3", "yellow"))(100),
            hclustfun = function(x) hclust(x, method = "complete"),
            distfun = function(x) as.dist(1-cor(t(x))),
            labRow = F, labCol = F
            )
}

################ 3 - CORRELATION NETWORK ##########

evaluate.components <- function(rho, diss.matrix) {
  print(paste("Testing rho value: ", as.character(rho)))
  # filter the adjacency matrix and create a new graph
  gene.network <- replace(diss.matrix, abs(1-diss.matrix) <= rho, 0)
  
  # removing nodes under the threshold and transforming everything into an adjacency list
  gene.network[upper.tri(gene.network)] <- 0 
  gene.adj.list <- melt(gene.network)
  gene.adj.list <- filter(gene.adj.list, abs(1-value) >= rho) %>% filter(value != 0) %>% filter(Var1 != Var2) 
  gene.graph <- graph_from_data_frame(gene.adj.list, directed=F, vertices = NULL)
  
  # gene.network <- subset(diss.matrix, abs(1-diss.matrix) >= rho)
  # nonzero.genes <- row.names(subset(gene.network, rowSums(gene.network) != 0))
  
  # gene.network <- gene.network[nonzero.genes, nonzero.genes]
  
  gene.graph.connect <- graph_from_data_frame(gene.adj.list, directed=F, vertices = NULL)
  
  distrib.deg <- degree_distribution(gene.graph, cumulative = T)
  degree.dist <- data.frame(freq=distrib.deg[-1], x=1:(length(distrib.deg)-1))
  degree.dist <- subset(degree.dist, freq != 0)
  # degree.dist$freq <- degree.dist$freq/sum(degree.dist$freq)
  fitting <- lm(formula=log10(freq) ~ I(log10(x)), data=degree.dist)
  
  # plot(freq ~ x, log = "yx", data=degree.dist,
  #      main=paste("rho:", rho, ", gamma:",round(fitting$coefficients[2] - 1, 2)))
  # abline(fitting$coefficients, col="red", lwd=2)
  
  alpha <- fitting$coefficients[2]
  r.squared <- - sign(alpha) * summary(fitting)$r.squared
  
  # nnls.fitting <- nnls(as.matrix(-log10(degree.dist$x+1)), 
  #                      as.matrix(-log10(degree.dist$freq)))
  # 
  # plot(-log10(freq) ~ -log10(x + 1), data=degree.dist,
  #       main=paste("rho:", rho, ", gamma:",round(nnls.fitting$x, 2)))
  # abline(b=nnls.fitting$x, a=-log10(mean(degree.dist$freq)), col="red", lwd=3)
  # 
  # r.squared <- nnls.fitting$x
  return(c(max(components(gene.graph.connect)$csize)/length(V(gene.graph.connect)), r.squared))
}

connectivity.check <- function(diss.matrix, min.rho, max.rho, step.rho, actual.rho, parallel.exec=F) {
  # checks the connectivy of the network at varying pearson correlation coefficient threshold
  rhos <- seq(min.rho, max.rho, step.rho)
  # par(mfrow=c(5,4), mar=c(2,2,2,2), cex=0.6)
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
  par(mfrow=c(1,2))
  plot(big.comp.frac~rhos, type='l', pch=19, col="blue",
       main="Biggest connected component", ylim=c(0, 1.2))
  abline(v=actual.rho, col="red")
  axis(1, at=seq(min.rho, max.rho, 0.1))
  plot(scalefreeness~rhos, type='l', pch=19, col="red",
       main="R^2 of fit to scale-free model",ylim=c(-0.7,1),xlim=c(0,1))
  abline(h=0.8, col="green" )
  return(evaluation)
}


# HERE LIES EVALUATE.PEARS.THRESH FUNCTION

draw.pears.dist <- function(diss.matrix, threshold) {
  # plots a histogram of pearson coefficients distribution and sets a threshold for the |cor(x,y)| value of each pair
  base.hist <- hist(diss.matrix, breaks = 100)
  quantile.breaks <- quantile(base.hist$breaks, c(1-threshold, threshold))
  hist.colors <- ifelse(base.hist$breaks <= quantile.breaks[1] |
                          base.hist$breaks >= quantile.breaks[2], "red", "gray")
  plot(base.hist, col=hist.colors)
}

filter.pears <- function(diss.matrix, threshold) {
  #filters connections based on a threshold, if |cor(x,y)| < threshold, then (x,y) = 0
  gene.net.raw <- replace(diss.matrix, abs(1-diss.matrix) <= quantile(1-diss.matrix, threshold), NA)
  # removal of disconnected nodes (which lead to error downward)
  gene.network <- gene.net.raw[!is.nan(sum(gene.net.raw, na.rm = T)),]
  return(gene.network)
}

mean.edge <- function(vect) {
  # helper function to calculate mean of edges (not considering 0s)
  mean(vect[which(!is.na(vect) & vect != 0)])
}

apcc <- function(gene.network){
  gene.network.clean <- gene.network
  gene.network.clean[is.na(gene.network.clean)] <- 0
  gene.network.clean <- gene.network.clean[colSums(gene.network.clean > 0) > 5,]
  apcc <- apply(gene.network, 1, mean.edge)
  apcc.clean <- apply(gene.network.clean, 1, mean.edge)
  plot(density(apcc.clean, n = 100), main = "APCC density distribution")
  return(apcc)
}

kmeans.scree <- function(gene.network, max.iter, n.clusters, n.repl) {
  scree <- list()
  gene.network[is.na(gene.network)] <- 0
  for (clusters in 1:n.clusters) {
    scree[[as.character(clusters)]] <- kmeans(gene.network, clusters, iter.max = max.iter, nstart = n.repl, trace = T)
  }
  scree.plot.data <- unlist(lapply(scree, function(X){X$tot.withinss}))
  
  plot(scree.plot.data, 
       main = "Scree plot of k-means clustering", 
       ylab = "Within group sum squared error", 
       xlab = "# of clusters",
       type="b")
  axis(1, at = 1:n.clusters)
  return(scree)
}

############### 5 - switch mining #########

normalize <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}

cRamp <- function(x){
  cols <- colorRamp(matlab.like(100))(normalize(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}

# Cartographic map
# PROBLEMA DI MERGING DI APCC E CLUSTEROPHOBICITY/WMD (risolto <3)
cartographic.map <- function(gene.network, communities, g.apcc, ...) {
  # plots a cartographic map AND returns a cartography list 
  clust.coef <- clusterophobic.coefficient(gene.network = gene.network, communities = communities)
  within.mod.deg <- within.module.degree(gene.network = gene.network, communities = communities)
  degree <- colSums(!is.na(gene.network), na.rm = T)
  g.apcc <- g.apcc[names(clust.coef)]
  cartography.df <- data.frame("clusterophobicity"=clust.coef, "wmdegree"=within.mod.deg, "gene.apcc"=g.apcc, "degree"=degree)
  palette(matlab.like(100))
  plot(wmdegree~clusterophobicity, 
       pch=19, 
       col=cRamp(g.apcc), 
       ylim = c(-5, 2), 
       data = cartography.df, 
       ylab="Within module degree",
       xlab="Clusterophobic coefficient",
       ...)
  return(cartography.df)
}

.clusterophobic.coeff <- function(genes, clust.names) {
  in.clust <- sum(genes[clust.names] != 0, na.rm =T)
  tot.deg <- sum(genes != 0, na.rm= T)
  return(1 - (in.clust/tot.deg)^2)
}

clusterophobic.coefficient <- function(gene.network, communities) {
  # returns a list of clust. coeffs
  coefs <- c()
  for (clust.num in 1:max(communities$cluster)) {
  nodes <- names(which(communities$cluster == clust.num))
  clust.network <- gene.network[nodes,]
  cluster.coef <- apply(clust.network, 1, .clusterophobic.coeff, clust.names=nodes)
  coefs <- c(coefs, cluster.coef)
  }
  coefs[is.nan(coefs)] <- 0
  return(coefs)
}

within.module.degree <- function(gene.network, communities) {
  # calculate mean and std dev for each cluster
  network <- c()
  
  for (clust.num in 1:max(communities$cluster)) {
    print(paste("calculating cluster number: ", clust.num))
    nodes <- names(which(communities$cluster == clust.num))
    clust.network.within <- gene.network[nodes, nodes]
    clust.network <- gene.network[nodes,]
    clust.degrees <- apply(clust.network != 0, 1, sum, na.rm=T)
    within.degrees <- apply(clust.network.within != 0, 1, sum, na.rm=T)
    # list of within module degrees
    clust.mean <- mean(clust.degrees)
    clust.std <- sd(clust.degrees)
    within.mod.deg <- (within.degrees - clust.mean)/clust.std
    network <- c(network, within.mod.deg)
  }
  return(network)
}

# TODO NOT THE SAME, CLOSER BUT WITH ONE YELLOW STRIPE MORE...
switch.mine <- function(cartography, dataset="t.dataset") {
  # plots switch biclustering and APCC
  # returns a vector of genes
  switch.genes <- subset(cartography, (clusterophobicity >= 0.8) & (wmdegree <= 2.5) & (degree >= 5) & (gene.apcc < 0))
  switch.names <- row.names(switch.genes)
  biclustering(dataset$matrix[switch.names,])
  return(switch.names)
}

############### 6 - evaluation ##########

average.s.path <- function(x, g.graph) {
  asp.tot <- c() 
  for (gene in 1:length(x)) {
    graph=delete.vertices(graph = g.graph, v = x[1:gene])
    asp <- average.path.length(graph)
    asp.tot <- c(asp.tot, asp)
    cat("Deleting gene: ", gene, " - Avg path length: ", asp, '\n')
  }
  return(asp.tot)
}

robustness.check <- function(gene.graph, switches, cartography) {
  # performs a robustness analysis on:
  # - random removal
  # - hubs removal
  # - switch removals
  #
  # plots everything gracefully (# of nodes removed vs shortest path)
  
  # order by hubs, order by fight-club, order by date, order by party and order by random
  # for each removal calculate avg shortest path and add it to a vector
  cartography.hubs <- subset(cartography, degree >= 5)
  setorder(cartography, -degree)
  setorder(cartography.hubs, -degree)
  genes.hubs <- row.names(cartography.hubs)[1:length(switches)]
  genes.fight.club <- row.names(subset(cartography.hubs, gene.apcc < 0))[1:length(switches)]
  genes.date <- row.names(subset(cartography.hubs, gene.apcc >= 0 & gene.apcc < 0.5))[1:length(switches)]
  genes.party <- row.names(subset(cartography.hubs, gene.apcc >= 0.5))[1:length(switches)]
  
  genes.random <- sample(x=row.names(cartography), 
                         size=length(switches))
  
  # order by non-switch hubs, order by switches, order by random (already done)
  # for each removal calculate the average shortest path and add it to a vector
  non.switch.hubs <- row.names(cartography.hubs[!row.names(cartography.hubs) %in% switches,])[1:length(switches)]
  switch.hubs <- row.names(cartography.hubs[switches,])
  
  targets <- list("hubs"=genes.hubs, "fight-club"=genes.fight.club, "date"=genes.date, "party"=genes.party,
                  "random"=genes.random, "non-switch"=non.switch.hubs, "switch"=switch.hubs)
  
  results <- mclapply(targets, average.s.path, g.graph=gene.graph,
                      mc.cores = detectCores()-1, mc.silent = F, mc.preschedule = T)

  max.frac <- max(sapply(results, length))
  
  tot.nodes <- length(V(gene.graph))
  
  par(mfrow=c(1,2))
  # plot average sp against fraction of removed nodes
  plot(results$hubs,
       type='l',
       col='black',
       ylab="Average path length", 
       xlab="Fraction of removed nodes", 
       ylim = c(min(sapply(results, min)), max(sapply(results, max)))) 
  lines(results$`fight-club`, col="red", type='l')
  lines(results$date, col="purple", type='l')
  lines(results$party, col='blue', type='l')
  lines(results$random, col='green', type='l')
  # TODO NON RIESCO A FAR APPARIRE LA LEGENDA IN MANIERA COMPOSTA
  #legend('topleft', names(results[1:5]), lty=1, lwd=2, 
  #       col=c("black", "red", "green", "blue", "purple"),
  #       xjust=1, yjust=1, cex=0.5)
  # plot average sp against fraction of removed nodes
  plot(results$switch,
       type='l',
       col='blue',
       ylab="Average path length", 
       xlab="Fraction of removed nodes", 
       ylim = c(min(sapply(results, min)), max(sapply(results, max))))
  axis(1, at=seq(0, (max.frac/tot.nodes), 3/max.frac))
  lines(results$`non-switch`, col='red', type='l')
  lines(results$random, col='green', type='l')
  #legend('topleft',legend=names(results[5:7]), lty=1, lwd=2,
  #       col=c("blue", "red", "green"), xjust=1, yjust=1, cex=0.5)
  
  return(results)
}

##################### 7 - linear regression ###################


evaluate.pears.thresh <- function(diss.matrix) {
  # for each value of rho:
  # - filter genes
  # - get the degree distribution
  # - linear model (log(P(k)) ~ -gamma*log(k))
  # - plot R^2
  # - check connected component with check.connectivity (maybe to be modified)
  
  
}
