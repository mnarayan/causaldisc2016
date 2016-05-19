
# Test function for GIES ------------------------------
require(R.matlab)
Alldata <- readMat('~/datasets/RestingAndStimData_5.16.mat')

xdim1 <- dim(Alldata$Data.0)
xdim2 <- dim(Alldata$Data.1)
xdim3 <- dim(Alldata$Data.2)
xdim3 <- dim(Alldata$Data.3)
xdim4 <- dim(Alldata$Data.4)


summary(target.index)
print(target.index[0])
print(target.index[xdim1[1]+5])
print(target.index[xdim1[1]+xdim2[1]+5])

stabMatright <- matrix(0,0,xdim1[2])
for(subjno in 1:xdim1[3]){
    print(subjno)
   # t.start <- proc.time()                         # start timer
   # t.date <- format(Sys.time(), "%Y_%m_%d_%H.%M") # time stamp for file saving
    nresamples <-50
    X <- rbind(Alldata$Data.0[,,subjno],Alldata$Data.1[,,subjno],Alldata$Data.2[,,subjno])
    print(dim(X))
    targetNodes <- list(0,27,28)
    print(targetNodes)
    targetInt <- c(rep(1,xdim1[1]),rep(2,xdim2[1]),rep(3,xdim3[1]))
    summary(targetInt)
    outputFile <- paste('~/causaldisc2016/RestingRightStim_',subjno,format(Sys.time(),'_Stability_GIES_%Y_%m_%d_%H%M'))

    stabMat<- SubjStabilityGIES(X, nresamples=nresamples,targets=targetNodes, target.index=targetInt, filename=outputFile)
    stabMatright<- rbind(stabMatright,stabMat)
}
writeMat('~/causaldisc2016/RightStim_GIES.mat',stabMat=stabMatright)

stabMatleft <- matrix(0,0,xdim1[2])
for(subjno in 1:xdim1[3]){
    print(subjno)
   # t.start <- proc.time()                         # start timer
   # t.date <- format(Sys.time(), "%Y_%m_%d_%H.%M") # time stamp for file saving
    nresamples <-50
    X <- rbind(Alldata$Data.0[,,subjno],Alldata$Data.3[,,subjno],Alldata$Data.4[,,subjno])
    targetNodes <- list(0,9,10)
    print(targetNodes)
    targetInt <- c(rep(1,xdim1[1]),rep(2,xdim2[1]),rep(3,xdim3[1]))
    summary(targetInt)
    outputFile <- paste('~/causaldisc2016/RestingLeftStim_',subjno,format(Sys.time(),'_Stability_GIES_%Y_%m_%d_%H%M'))
    
    stabMat<- SubjStabilityGIES(X, nresamples=nresamples,targets=targetNodes, target.index=targetInt, filename=outputFile)
    stabMatleft<- rbind(stabMatleft,stabMat)
}    
writeMat('~/causaldisc2016/LeftStim_GIES.mat',stabMat=stabMatleft)



stabMatRest <- matrix(0,0,xdim1[2])
for(subjno in 1:xdim1[3]){
    print(subjno)
   # t.start <- proc.time()                         # start timer
   # t.date <- format(Sys.time(), "%Y_%m_%d_%H.%M") # time stamp for file saving
    nresamples <-50
    X <- Alldata$Data.0[,,subjno]
    outputFile <- paste('~/causaldisc2016/Resting_',subjno,format(Sys.time(),'_Stability_GES_%Y_%m_%d_%H%M'))    
    stabMat<- SubjStabilityGES(X, nresamples=nresamples, filename=outputFile)
    stabMatRest<- rbind(stabMatleft,stabMat)
}    
writeMat('~/causaldisc2016/Resting_Stability_GES.mat',stabMat=stabMatRest)


###################### OLDER CODE, Try BackShift Later ##################
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
