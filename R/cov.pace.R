# Codes in the file mostly are adapted from fdapace package developed by
# COPYRIGHT HOLDER: Hans-Georg Mueller and Jane-Ling Wang
# ORGANIZATION: University of California, Davis
# YEAR: 2019
#
# The support to functional snippets are added by Zhenhua Lin (National University of Singapore)
#
# and specify as
#
# License: BSD_3_clause + file LICENSE
#
#
# Copyright (c) <YEAR>, <COPYRIGHT HOLDER>
#
#     Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.
#
# Neither the name of the <ORGANIZATION> nor the names of its
# contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#' PACE approach to covariance estimation
#' @param Lt a list (for irregular design) or a vector (for regular design)
#' @param Ly a list (for irregular design) for a matrix (for regular design). If \code{Ly} is a matrix, then \code{ncol(Ly)} must be equal to \code{length(Lt)}
#' @param bw bandwidth or bandwidth candidates
#' @param delta the snippet parameter, only used for irregular design
#' @keywords internal
cov.pace <- function(Lt,Ly,bw=NULL,newt=NULL,mu=NULL,tuning='GCV',weig='SUBJ',kernel='epanechnikov',delta=NULL,
                     binData='AUTO',numBins=NULL,nRegGrid=51)
{
    if(!is.list(Lt) || !is.list(Ly)) stop('regular design is not supported yet')

    if(is.null(mu)) mu <- meanfunc(Lt,Ly,method='PACE')
    if(is.null(delta)) delta <- estimate.delta(Lt)

    # Bin the data
    if ( binData != 'OFF'){
        BinnedDataset <- GetBinnedDataset(Ly,Lt,numBins)
        Ly = BinnedDataset$newy;
        Lt = BinnedDataset$newt;
        nRegGrid <- min(nRegGrid,BinnedDataset[['numBins']])
    }

    covobj <- create.cov.pace.obj(Lt,Ly,bw,mu,weig,method='PACE',delta=delta,nRegGrid=nRegGrid,kernel=kernel)
    obsGrid <- sort(unique( c(unlist(Lt))))
    if(is.null(newt))  newt <- seq(min(obsGrid), max(obsGrid),length.out=nRegGrid)
    obj <- GetSmoothedCovarSurface(covobj, obsGrid, newt, kernel, useBinnedCov=FALSE)
    covobj$fitted <- obj$smoothCov
    covobj$bw <- obj$bw
    covobj$newt <- obj$outGrid
    class(covobj) <- 'covfunc'
    return(covobj)
}

create.cov.pace.obj <- function(Lt,Ly,bw,mu,weig,method,delta=1,nRegGrid=51,kernel='epanechnikov')
{
    eval.mean <- function(mu,x)
    {
        if(is.function(mu)) return(mu(x))
        else return(predict(mu,x))
    }

    if(is.null(mu)) mu <- meanfunc(Lt,Ly,method='PACE')
    if(is.null(weig)) weig <- 'SUBJ'
    weig <- get.weig.cov(Lt,Ly,scheme=weig)


    x1 <- list()
    x2 <- list()
    y <- list()
    w <- list()
    for(i in 1:length(Lt))
    {
        Ti <- pracma::meshgrid(Lt[[i]])
        mui <- eval.mean(mu,Lt[[i]])
        yc <- Ly[[i]] - mui
        y[[i]] <- as.vector(yc %*% t(yc))
        x1[[i]] <- as.vector(Ti$X)
        x2[[i]] <- as.vector(Ti$Y)
        w[[i]] <- rep(weig[i],length(y[[i]]))
    }

    x1 <- unlist(x1)
    x2 <- unlist(x2)
    y <- unlist(y)
    w <- unlist(w)

    # remove diagonal elements
    idx <- (x1 != x2)
    x1 <- x1[idx]
    x2 <- x2[idx]
    y <- y[idx]
    w <- w[idx]

    R <- list(kernel=kernel,x=cbind(x1,x2),y=y,weig=w,nRegGrid=nRegGrid,
              bw=bw,Lt=Lt,Ly=Ly,mu=mu,method=method,delta=delta)
    class(R) <- 'covfunc'
    return(R)
}

GetSmoothedCovarSurface <- function(covobj, obsGrid, regGrid, kernel='epanechnikov', useBinnedCov=FALSE) {

    methodBwCov <- 'GCV'

    # get the truncation of the output grids.
    buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
    rangeGrid <- range(regGrid)
    minGrid <- rangeGrid[1]
    maxGrid <- rangeGrid[2]
    cutRegGrid <- regGrid[regGrid > minGrid -
                              buff &
                              regGrid < minGrid + diff(rangeGrid) +
                              buff]

    rcov <- list(tPairs=covobj$x,cxxn=covobj$y,win=covobj$weig,dataType='Sparse')

    bw <- covobj$bw
    if(is.null(bw) || length(bw) > 1)
    {
        gcvObj <- GCVLwls2DV2(obsGrid, regGrid, kernel, rcov, covobj$Lt, bw, covobj$delta)
        bw <- gcvObj$h
    }



    #    if (!useBinnedCov) {

    smoothCov <- Lwls2D(bw, kernel, xin=rcov$tPairs, yin=rcov$cxxn,
                        xout1=cutRegGrid, xout2=cutRegGrid, win=covobj$weig,
                        delta=covobj$delta)


    res <- list(rawCov = rcov,
                smoothCov = (smoothCov + t(smoothCov)) / 2,
                bw = bw,
                outGrid = cutRegGrid)
    return(res)
}


# The following function is adapted from fdapace package
GCVLwls2DV2 <- function(obsGrid, regGrid, kernel, rcov, Lt, bw, delta) {

    verbose <- T

    ## Get minimal bandwidth and range
    r <- diff(range(obsGrid)) * sqrt(2) # sqrt(2) because the window is circular.
    minBW <- fdapace::BwNN(Lt, k= 3, onlyCov = TRUE)['cov']

    h0 <- minBW

    if (is.null(h0)){
        stop('the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n')
    }

    if (kernel == 'gauss') {
        h0 = h0 * 0.2;
    }

    ## Get Candidate Bandwidths
    h0 <- min(h0, r/4)
    if (h0 < r/4) {
        q <- (r / (4 * h0)) ^ (1/9)
    } else if (h0 < r/2) {
        q <- (r / (2 * h0)) ^ (1/9)
    } else if (h0 < r) {
        q <- (r / h0) ^ (1/9)
    } else {
        stop('Data is too sparse. The minimal bandwidth is the range of data')
    }

    if(is.null(bw)) bw <- (q ^ seq(0,9,length.out = 10)) * h0 # from h0 to r / 4


    opth <- bw[1]

    leave <- FALSE
    iter <- 0
    maxIter <- 1

    minBWInvalid <- FALSE
    while (!leave && iter < maxIter) {
        if (minBWInvalid){
            minBW <- bw[1]
        }

        #Scores <- rep(Inf, length(bw))
        Scores <- matrix(Inf, nrow = length(bw), ncol = 2); colnames(Scores) <- c('SUM','SD');
        # try the bandwidths large to small in order to save time due to sparseness in the windows.
        for (i in rev(seq_along(bw))) {
            h <- bw[i]

            #            if (class(rcov) == 'BinnedRawCov') {
            #                Scores[i,'SUM'] <- getGCVscoresV2(h, kernel, rcov$tPairs, rcov$meanVals, win=rcov$count, regGrid=regGrid, RSS=rcov$RSS, verbose=verbose)
            #            }
            #            else {
            Scores[i,'SUM'] <- getGCVscoresV2(h, kernel, rcov$tPairs, rcov$cxxn, regGrid, delta)
            #            }

            if (is.infinite(Scores[i,'SUM'])) {
                minBWInvalid <- TRUE
                if (i < length(bw)) {
                    if (minBWInvalid) {
                        minBW <- bw[i + 1]
                        minBWInvalid <- FALSE
                    }
                }
                break; # This will help break out of the loop if the BW is too small to make sense
            }
        }


        if(is.infinite(min(Scores[,'SUM']))){
            opth <- max(bw)
            optgcv <- Inf
        } else {
            {
                ind <- which.min(Scores[,'SUM'])
                opth <- bw[ind]
                optgcv <- Scores[ind,'SUM']
            }
        }

        ## Check that what we found is coherent.
        if (opth >= r - 1e-12) {
            minBW <- r
            leave <- TRUE
            stop('Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.')
        }
        if ( (abs(opth - max(bw)) > 1e-12) && !is.infinite(optgcv))
            leave <- TRUE
        else if (is.infinite(optgcv)) {
            if (verbose)
                warning('Data is too sparse, retry with larger bandwidths!')
            h0 <- max(bw) * 1.01
        } else if ( (abs(opth - max(bw)) > 1e-12) ) {
            warning('Optimal bandwidth not found in the candidate bandwidths. Retry with larger bandwidths')
            h0 <- max(bw)
        }
        if (!leave) {
            newr <- seq(0.5, 1, by=0.05) * r # ??? this can be quite slow
            ind <- which(newr > h0)[1]
            q <- (newr[ind] / h0) ^ (1/9)
            bw <- q ^ (0:9) * h0
            if (verbose) {
                message('New bwuserCov candidates:\n')
                print(bw)
            }
            iter <- iter + 1
        }
    } # The "while (!leave && iter < maxIter) ..." end here

    ret <- list(h=opth, gcv=optgcv, minBW=minBW)
    message(sprintf('optimal h=%f',opth))
    return(ret)

}


# The following function is adapted from fdapace package
getGCVscoresV2 <- function(h, kernel, xin, yin, regGrid, delta, win=NULL) {

    verbose <- T
    if(is.null(win))
        win <- rep(1, length(yin))

    fit <- tryCatch(Lwls2D(h, kernel, xin=xin, yin=yin,
                               win=win, xout1=regGrid,
                               xout2=regGrid,delta=delta), error=function(err) {
            if (verbose) {
                message('Invalid bandwidth. Try enlarging the window size.\n')
            }
            return(Inf)
        })


    # workaround for degenerate case.
    if (any(is.nan(fit)) || !is.matrix(fit))
        return(Inf)



    obsFit <- pracma::interp2(regGrid, regGrid, fit, xin[, 1], xin[, 2])

    idx <- !is.infinite(obsFit)
    idx <- !is.nan(obsFit) & idx

    # residual sum of squares
    res <- sum((yin[idx] - obsFit[idx]) ^ 2 * win[idx])
    N <- sum(idx)

    # kernel at x=0
    if(kernel=='epanechnikov')
        k0 <- 0.75
    else if(kernel=='gauss')
        k0 <- 0.398942280401433
    else stop('unsupported kernel')

    r <- diff(range(xin[, 1]))
    bottom <- max(1 - 3 * (1 / N) * (r * k0 / h)^2, 0)
    GCV <- res / bottom^2

    return(GCV)
}


GetBinNum = function(n, m ){

    # Get the number of bins
    # n : number of curves
    # m : median or max value of number of time-points

    numBin = NULL;
    if (m <= 20){
        return(NULL)
    }


    if (m >400){
        numBin = 400;
    }

    if (n > 5000){
        mstar = max(20,(((5000-n)*19)/2250)+400);
        if (mstar < m){
            numBin = ceiling(mstar);
        }
        else
        {
            return(NULL)
        }
    }

    if( is.null(numBin) ) {
        message('No binning is needed!\n');
    }

    return(numBin)
}

GetBinnedCurve <- function(x, y, M = 10, limits = c(min(x),max(x)) ){

    # x : 1 * n vector of values to be binned
    # y : 1 * n vector of values corresponding to x
    # M : positive interger denotes the number of bins to be use  or
    #             positive value of the width of each equidistant bin
    # isMnumBin : use number of bin (TRUE) or bandwith h (FALSE)
    # nonEmptyOnly : output only non-empty bins (TRUE)
    # limits : lower and upper limit of domain x (a0, b0)


    # Function 'GetBinnedCurve' starts here
    if( M<0 ){
        stop("GetBinnedCurve is aborted because the argument: M is invalid!\n");
    }


    h <- diff(limits) / M
    xx <- c(limits[1],seq(limits[1] + h / 2, limits[2] - h / 2, length.out=M - 1),limits[2])
    N <- length(xx);
    midpoint <- seq(limits[1], limits[2], length.out=M)

    zList = GetBins(x,y,xx)
    newy = zList$newy
    count = zList$count

    midpoint = midpoint[count > 0]
    newy = newy[count > 0]
    count = count[count > 0]
    M = length(midpoint)


    res = list(midpoint = midpoint, newy = newy, count = count, numBin = M, binWidth = h)
    return(res);
}

# Auxilary function 'GetBins'
GetBins <-  function(x,y, xx){
    N <- length(xx)
    count = rep(0,N-1);
    newy = count;

    for (i in  2:(N-1)){
        ids = ((x >= xx[i-1]) &  (x < xx[i]));
        if (all(ids == 0)){
            count[i-1] = 0;
            newy[i-1] = NaN;
        } else {
            count[i-1] =   sum(ids);
            newy[i-1] =  mean(y[ids]);
        }
    }

    # print('GetBins used')
    # for the last bin, include the left and right end point
    ids =  ((x >= xx[i]) &(x <= xx[i+1]));
    if (all(ids == 0)){
        count[i] = 0;
        newy[i] = NaN;
    }else {
        count[i] = sum(ids);
        newy[i] = mean(y[ids]);
    }

    zList = list(newy = newy, count = count)
    return( zList );
}


GetBinnedDataset <- function (Ly, Lt, numBins=NULL){

    # Bin the data 'y'
    # y : n-by-1 list of vectors
    # t : n-by-1 list of vectors

    BinDataOutput <- list( newy=NULL, newt=NULL);

    tt = unlist(Lt);
    a0 = min(tt);
    b0 = max(tt);

    n = length(Lt);
    mi = sapply(Lt,length);

    m = median(mi)


    # Determine the number of bins automatically if numBins is null
    if (is.null(numBins)){
        numBins = GetBinNum(n,m)
        if (is.null(numBins)){
            BinDataOutput$newt = Lt
            BinDataOutput$newy = Ly
            return( BinDataOutput )
        }
    }

    numBins = ceiling(numBins)

    resList <- lapply(1:n, function(i)
        GetBinnedCurve(Lt[[i]], Ly[[i]], numBins, c(a0, b0)))
    BinDataOutput[['newt']] <- lapply(resList, `[[`, 'midpoint')
    BinDataOutput[['newy']] <- lapply(resList, `[[`, 'newy')


    result <- list( 'newt' = BinDataOutput$newt, 'newy' = BinDataOutput$newy,
                    numBins = numBins)
    return(result)
}

# the following function is adapted from fdapace
Lwls2D <- function(bw, kernel='epanechnikov', xin, yin, win=NULL, xout1=NULL, xout2=NULL,
                   subset=NULL,method = ifelse(kernel == 'gauss', 'plain', 'sort2'),
                   delta=1) {


    if (length(bw) == 1){
        bw <- c(bw, bw)
    }
    if (!is.matrix(xin) ||  (dim(xin)[2] != 2) ){
        stop('xin needs to be a n by 2 matrix')
    }
    # xin <- matrix(xin, ncol=2) # This causes unexcepted/wrong results.

    if (is.null(win)){
        win <- rep(1, nrow(xin))
    }
    if (!(nrow(xin) == length(yin) && nrow(xin) == length(win))) {
        stop('The length of xin, yin, and win (if specified) must be the same')
    }
    if (!is.null(subset)) {
        xin <- xin[subset, ]
        yin <- yin[subset]
        win <- win[subset]
    }

    if (is.null(xout1))
        xout1 <- sort(unique(xin[, 1]))

    if (is.null(xout2))
        xout2 <- sort(unique(xin[, 2]))

    # For passing numerics into the cpp smoother.
    storage.mode(bw) <- 'numeric'
    storage.mode(xin) <- 'numeric'
    storage.mode(yin) <- 'numeric'
    storage.mode(win) <- 'numeric'
    storage.mode(xout1) <- 'numeric'
    storage.mode(xout2) <- 'numeric'


    if(method == 'plain'){
        ret <- csmoothcov(bw, kernel, xin, yin, win,xout1, xout2, delta)
    } else if (method == 'sort2'){
        ord <- order(xin[, 1])
        xin <- xin[ord, ]
        yin <- yin[ord]
        win <- win[ord]
        # browser()
        ret <- csmoothcov(bw, kernel, xin, yin, win,xout1, xout2, delta)
    }


    return(ret)
}
