matrixcov <- function(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman")){
xdim <- dim(x)
tmpCov <- matrix(0,xdim[2],xdim[2])
    for(ii in 1:xdim[3]){
        tmpVar <- var(x[,,ii])
        diag(tmpVar)
        tmpCov <- tmpCov + cov(x[,,ii] %*% solve(sqrt(diag(diag(tmpVar)))),method=method)
    }
    tmpCov <- tmpCov/xdim[3]
    return(tmpCov)
}

matrixcorr <- function(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman")){
xdim <- dim(x)
Sighat <- matrixcov(x,method=method)
Sigvar <- diag(Sighat)
Sigvar
Varinv <- solve(diag(sqrt(Sigvar)))
Sighat2 <- Varinv %*%Sighat%*% t(Varinv)
    return(Sighat2)
}


# Call PCAlg across subject resamples, aggregate stability matrix and then export stability adjacency matrix
StabilityPC <- function(X, nresamples=50,
                  indepTest = gaussCItest, ## indep.test: partial correlations
                  maj.rule = TRUE,
                  alpha=0.01, labels, verbose = FALSE,
                  solve.confl=TRUE){

    # Package dependencies ------------------------------
    require(pcalg)
    require(Rgraphviz)
    require(igraph)
    # Initial Parameters ------------------------------
    xdim <- dim(X)
    xdim
    subSize <- floor(xdim[3]/2)
    stabMat <- matrix(0,xdim[2],xdim[2])
    roinames <- read.csv('~/datasets/roinames.csv')
    #c(colnames(roinames))
    nodelabels <- c(colnames(roinames))
    colnames(stabMat) <- nodelabels[1:xdim[2]]
    rownames(stabMat) <- nodelabels[1:xdim[2]]
    # Resampling/Stability Matrices ------------------------------
    for(b in 1:nresamples){
        print(b)
        subIndices <- sample(xdim[3], subSize, replace = FALSE)
        subX <- X[,,subIndices] # subsample the data set
        subXdim <- dim(subX)

        ## PC Algorithm ------------------------------
        ## estimate CPDAG on subsample
        Sighat <- matrixcov(subX,method="pearson")
        roinames <- read.csv('~/datasets/roinames.csv')
        nodelabels <- c(colnames(roinames))
        suffStat <- list(C = Sighat, n = subXdim[1]*subXdim[3])

        ## estimate CPDAG
        pc.fit <- pc(suffStat,
                  indepTest=gaussCItest, ## indep.test: partial correlations
                  maj.rule = TRUE,
                  alpha=0.01, labels = nodelabels[1:subXdim[2]], verbose = FALSE,
                  solve.confl=TRUE)

        stabMat <- stabMat + as.matrix(as_adjacency_matrix(graph_from_graphnel(attributes(pc.fit)$graph)))
        #print(stabMat[1:5,1:5])

    }
    write.csv(stabMat,format(Sys.time(),'causaldisc2016/Resting_Stability_PC_%Y_%m_%d_%H%M'))
    return(stabMat)
}

# Call GES across subject resamples, aggregate stability matrix and then export stability adjacency matrix
StabilityGES <- function(X, nresamples=50,
                  indepTest = gaussCItest, ## indep.test: partial correlations
                  maj.rule = TRUE,
                  alpha=0.01, labels, verbose = FALSE,
                  solve.confl=TRUE){

    # Package dependencies ------------------------------
    require(pcalg)
    require(Rgraphviz)
    require(igraph)
    # Initial Parameters ------------------------------
    xdim <- dim(X)
    xdim
    subSize <- floor(xdim[3]/2)
    stabMat <- matrix(0,xdim[2],xdim[2])
    roinames <- read.csv('~/datasets/roinames.csv')
    #c(colnames(roinames))
    nodelabels <- c(colnames(roinames))
    colnames(stabMat) <- nodelabels[1:xdim[2]]
    rownames(stabMat) <- nodelabels[1:xdim[2]]
    # Resampling/Stability Matrices ------------------------------
    for(b in 1:nresamples){
        print(b)
        subIndices <- sample(xdim[3], subSize, replace = FALSE)
        subX <- X[,,subIndices] # subsample the data set
        subXdim <- dim(subX)

        ## PC Algorithm ------------------------------
        ## estimate CPDAG on subsample
        Sighat <- matrixcov(subX,method="pearson")
        roinames <- read.csv('~/datasets/roinames.csv')
        nodelabels <- c(colnames(roinames))
        suffStat <- list(C = Sighat, n = subXdim[1]*subXdim[3])

        ## estimate CPDAG
        ges.fit <- pc(suffStat,
                  indepTest=gaussCItest, ## indep.test: partial correlations
                  maj.rule = TRUE,
                  alpha=0.01, labels = nodelabels[1:subXdim[2]], verbose = FALSE,
                  solve.confl=TRUE)

        stabMat <- stabMat + as.matrix(as_adjacency_matrix(graph_from_graphnel(attributes(ges.fit)$graph)))
        #print(stabMat[1:5,1:5])

    }
    write.csv(stabMat,format(Sys.time(),'causaldisc2016/Resting_Stability_GES_%Y_%m_%d_%H%M'))
    return(stabMat)
}

# Call GIES across subject resamples, aggregate stability matrix and then export stability adjacency matrix
StabilityGIES <- function(X, nresamples=50,
                  targets, target.index, filename='ObsAndInt_Stability_GIES'){

    # Package dependencies ------------------------------
    require(pcalg)
    require(Rgraphviz)
    require(igraph)
    # Initial Parameters ------------------------------
    xdim <- dim(X)
    xdim
    subSize <- floor(xdim[3]/2)
    stabMat <- matrix(0,xdim[2],xdim[2])
    roinames <- read.csv('~/datasets/roinames.csv')
    #c(colnames(roinames))
    nodelabels <- c(colnames(roinames))
    colnames(stabMat) <- nodelabels[1:xdim[2]]
    rownames(stabMat) <- nodelabels[1:xdim[2]]
    # Resampling/Stability Matrices ------------------------------
    for(b in 1:nresamples){
        print(b)
        subIndices <- sample(xdim[3], subSize, replace = FALSE)
        subX <- X[,,subIndices] # subsample the data set
        subXdim <- dim(subX)

        ## GIES Algorithm ------------------------------
        ## estimate CPDAG on subsample
        subScore <- new("GaussL0penIntScore", subX, targets=targets, target.index=target.index)
        
        ## estimate CPDAG
        gies.fit <- gies(score)

        stabMat <- stabMat + as.matrix(as_adjacency_matrix(graph_from_graphnel(attributes(gies.fit)$essgraph)))
        #print(stabMat[1:5,1:5])

    }
    write.csv(stabMat,format(Sys.time(),'causaldisc2016/RestingStim_Stability_GIES_%Y_%m_%d_%H%M'))
    return(stabMat)
}



