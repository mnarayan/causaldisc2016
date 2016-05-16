library(R.Matlab)
aMFG <- readMat('pcnets/SP_aMFG_4.1.mat)
pMFG <- readMat('pcnets/SP_pMFG_4.2.mat)

X_aMFG <- aMFG$Data
X_pMFG <- pMFG$Data

X <- X_pMFG[,,2]
library(CompareCausalNetworks)

methods <- c("pc")
runStability <- TRUE
sq <- ceiling(sqrt(length(methods)+1))
par(mfrow=c(ceiling((length(methods)+1)/sq),sq))
V <- colnames(read.csv('pcnets/roinames.csv'))
#Note that calling (conservative = TRUE), or maj.rule = TRUE, together with solve.confl = TRUE produces a fully order-independent output, see Colombo and Maathuis (2013).
pc.fit <- pc(suffStat = list(C=cor(X),n=nrow(X)), indepTest=gaussCItest, alpha=.01,labels=V, verbose=TRUE);
plot(pc.fit, main = "Subj1, PC")
skel.fit <- skeleton(suffStat = list(C=cor(X),n=nrow(X)),indepTest=gaussCItest, alpha=.01, p=ncol(X), verbose=FALSE) ##u2pd="relaxed",skel.method="stable", verbose=FALSE)
helpplot(skel.fit)

environment<- rep(1,dim(X)[1])
## loop over all methods and compute and print/plot estimate
for (method in methods){
  cat("\n result for method", method,"  ------  \n" )
  
  if(!runStability){
    # Option 1): use this estimator as a point estimate
    Ahat <- getParents(X, environment, method=method, alpha=0.1, pointConf = TRUE)
  }else{
    # Option 2): use a stability selection based estimator
    # with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment, EV=2, method=method, alpha=0.1)
  }
  
  # print and plot estimate (point estimate thresholded if numerical estimates
  # are returned)
  print(Ahat)
  if(!runStability)
    plot(as_graphnel(graph_from_adjacency_matrix(as.matrix(Ahat))))
    #plotGraphEdgeAttr(Ahat, plotStabSelec = FALSE, labels = labels,
    #                 thres.point = 0.05,
    #                  main=paste("POINT ESTIMATE FOR METHOD\n", toupper(method)))
  else
    plot(as_graphnel(graph_from_adjacency_matrix(as.matrix(Ahat))))
    #plotGraphEdgeAttr(Ahat, plotStabSelec = TRUE, labels = labels,
    #                  thres.point = 0, main = paste("STABILITY SELECTION
    # ESTIMATE\n FOR METHOD", toupper(method)))
}
