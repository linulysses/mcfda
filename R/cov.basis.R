#' Covariance estimation by basis expansion
#' @param Lt a list (for irregular design) or a vector (for regular design)
#' @param Ly a list (for irregular design) for a matrix (for regular design). If \code{Ly} is a matrix, then \code{ncol(Ly)} must be equal to \code{length(Lt)}
#' @keywords internal
cov.basis <- function(Lt,Ly,p=NULL,lam=NULL,ext=NULL,newt=NULL,mu=NULL,tuning='lle',weig=NULL,maxIt=3,domain=c(0,1))
{
    if(is.null(p)) p <- seq(5,7)
    if(is.null(lam)) lam <- 10^seq(-8,1,length.out=10)
    if(is.null(ext)) ext <- c(0,0.05,0.1)


    if(is.list(Ly))
    {
        n <- length(Ly)
        datatype <- 'irregular'
        if(is.null(weig)) weig <- get.weig.cov(Lt,Ly,'OBS')
    }
    else
    {
        n <- nrow(Ly)
        datatype <- 'regular'
        weig <- rep(1/n,n)
    }


    if(is.null(mu)) mu <- meanfunc(Lt,Ly,method='FOURIER')

    if(is.null(domain))
    {
        tmp <- unique(unlist(Lt))
        domain <- c(min(tmp),max(tmp))
    }

    if(domain[1] != 0 || domain[2] != 1) # standardize tobs
    {
        if(datatype=='irregular')
            Lt <- lapply(Lt,function(x) (x-domain[1])/(domain[2]-domain[1]))
        else
            Lt <- (Lt-domain[1])/(domain[2]-domain[1])
    }

    R <- list(C=NULL,p=p,lam=lam,ext=ext,mu=mu,datatype=datatype,domain=domain,
              maxIt=maxIt,Lt=Lt,Ly=Ly,weig=weig,method='FOURIER')
    class(R) <- 'covfunc'

    if(is.list(Ly) && is.list(Lt)) # irregular design
    {
        Lr <- lapply(1:n,function(i){
            mui <- predict(mu,Lt[[i]])
            (Ly[[i]]-mui) %*% t(Ly[[i]]-mui)
        })
        delta <- estimate.delta(Lt)

        R$Lr <- Lr
        R$delta <- delta

    }
    else if(is.vector(Lt) && is.matrix(Ly)) # regular design
    {
        mui <- predict(mu,Lt)
        Lr <- lapply(1:n,function(i){
            (Ly[i,]-mui) %*% t(Ly[i,]-mui)
        })

        R$Lr <- Lr
        R$delta <- 1
    }
    else stop('unsupported data type')

    if(length(p) > 1 || length(lam) > 1 || length(ext) > 1)
    {
        if(tuning=='CV') stop('unsupported tuning method')
        else if(tuning=='GCV') stop('unsupported tuning method')
        else
        {
            R <- cov.basis.lle(R,p,lam,ext)$R # select tuning parameters
        }
    }

    lam <- R$lam
    p <- R$p
    ext <- R$ext

    if(datatype=='irregular')
    {
        if(ext > 0) Ls <- shrink(Lt,ext)
        else Ls <- Lt
        auxmat <- comp.aux.mat(p,Ls)
        C <- estimate.cov.coef(Ls,Lr,auxmat,lam,weig,maxIt)
    }
    else
    {
        Lt <- list(Lt)
        if(ext > 0) Ls <- shrink(Lt,ext)
        else Ls <- Lt
        auxmat <- comp.aux.mat(p,Ls)

        Lr <- list(Reduce('+',Lr)/length(Lr))
        C <- estimate.cov.coef(Ls,Lr,auxmat,lam,weig,maxIt)
    }


    R$C <- C

    if(!is.null(newt)) R$fitted <- predict.covfunc(R,newt)

    return(R)

}

cov.basis.lle <- function(R,p,lam,ext)
{

    if(R$datatype == 'irregular')
    {
        xcov <- cov.pace(R$Lt,R$Ly,bw=NULL,newt=NULL,mu=R$mu,
                         tuning='GCV',weig='SUBJ',kernel='epanechnikov',
                         delta=R$delta)
        newt <- xcov$newt
        tmp <- pracma::meshgrid(newt)
        x1 <- as.vector(tmp$X)
        x2 <- as.vector(tmp$Y)
        idx <- abs(x1-x2) < R$delta
        idx <- idx & (x1!=x2)
        y <- as.vector(xcov$fitted)[idx]

        Lr <- R$Lr
        Lt <- R$Lt
        Ly <- R$Ly
    }
    else
    {
        newt <- R$Lt
        n <- length(R$Lr)
        tr <- sample(n,ceiling(0.75*n))
        te <- setdiff(1:n,tr)

        tmp <- R$Lr[te]
        y <- Reduce('+',tmp) / length(te)
        y <- as.vector(y)
        idx <- rep(TRUE,length(y))

        Lr <- list(Reduce('+',R$Lr[tr])/length(tr))
        Lt <- list(R$Lt)

    }



    # find p and lambda first, then ext
    R$ext <- 0
    if(length(p) > 1 || length(lam) > 1)
    {
        err <- matrix(Inf,length(p),length(lam))
        for(u in 1:length(p))
        {
            auxmat <- comp.aux.mat(p[u],Lt)
            R$p <- p[u]
            for(v in 1:length(lam))
            {
                C <- estimate.cov.coef(Lt,Lr,auxmat,lam[v],R$weig,maxIt=R$maxIt)
                R$C <- C
                R$lam <- lam[v]
                covhat <- predict.covfunc(R,newt)
                yhat <- as.vector(covhat)
                yhat <- yhat[idx]
                err[u,v] <- mean((y-yhat)^2)
            }
        }

        min.idx <- which(err==min(err),arr.ind=T)
        p <- p[min.idx[1]]
        lam <- lam[min.idx[2]]
    }

    R$p <- p
    R$lam <- lam

    if(length(ext) > 1)
    {
        err <- rep(0,length(ext))
        for(u in 1:length(ext))
        {
            Ls <- shrink(Lt,ext[u])
            auxmax <- comp.aux.mat(p,Ls)
            C <- estimate.cov.coef(Lt,Lr,auxmat,lam,R$weig,maxIt=R$maxIt)
            R$ext <- ext[u]
            covhat <- predict.covfunc(R,newt)
            yhat <- as.vector(covhat)
            yhat <- yhat[idx]
            err[u] <- mean((y-yhat)^2)
        }
        ext <- ext[which.min(err)]
    }
    else
    {
        R$ext <- ext
    }

    return(list(R=R,ext=ext,p=p,lam=lam))
}


comp.aux.mat <-function(p,Lt)
{
    if(is.list(Lt))
    {
        n <- length(Lt)
        B <- lapply(Lt, function(v) evaluate.basis(p,grid=v,type = 'FOURIER'))
    }
#    else if(is.vector(Lt))
#    {
#        B <- evaluate.basis(p,grid=Lt,type = 'FOURIER')
#    }
    else stop('unsupported type of Lt')

    m <- 1000
    pts <- seq(0.5/m,1-0.5/m,length.out=m)
    W <- deriv.fourier(p,grid=pts,r=2)
    W <- t(W) %*% W / m
    V <- deriv.fourier(p,grid=pts,r=1)
    V <- t(V) %*% V / m
    U <- deriv.fourier(p,grid=pts,r=0)
    U <- t(U) %*% U / m

    return(list(B=B,W=W,V=V,U=U))
}

shrink <- function(Lt,ext)
{
    if(is.list(Lt))
    {
        lapply(Lt,function(v) ext + (1-ext*2)*v)
    }
    else if(is.vector(Lt)) ext + (1-ext*2)*Lt
    else stop('unsupported type of Lt')
}

estimate.cov.coef <- function(Lt,Lr,auxmat,lam,weight,maxIt=3,method='manifold')
{
    Chat <- list()
    if(method == 'naive')
    {
        if(length(lam)==1) Chat <- minimizeL(Lt,Lr,auxmat,lam,weight)$C0
        else Chat <- lapply(lam,function(lam0) minimizeL(Lt,Lr,auxmat,lam0,weight))$C0
    }
    else if(method == 'manifold')
    {
        if(length(lam)==1)
        {
            Chat <- minimizeL(Lt,Lr,auxmat,lam,weight)$C0
            eigval <- eigen(Chat)$values
            rho <- min(eigval)
            if(rho < 0)
            {
                tmp <- Chat + abs(rho*1.1)*diag(rep(1,dim(Chat)[1]))
                tmp <- t(chol(tmp))
                Lmin <- minimizeQ(Lt,Lr,auxmat,lam,weight,
                                  L0=tmp,maxIt=maxIt)$L
                Chat <- Lmin %*% t(Lmin)
            }
        }
        else stop('multiple lambda is not supported yet.')
    }
    else stop('unsupported method')

    return(Chat)
}

minimizeL <- function(Lt,Lr,auxmat,lam,weight)
{
    n <- length(Lt)
    p <- dim(auxmat$U)[1]

    ## compute D and E
    E <- matrix(0,p*p,1)
    D <- matrix(0,p*p,p*p)
    for(i in 1:n)
    {
        Ni <- length(Lt[[i]])
        if(Ni > 1)
        {
            Ri <- Lr[[i]]
            Bi <- auxmat$B[[i]]

            # Term I
            E <- E - weight[i] * kronecker(t(Bi),t(Bi)) %*% matrix(Ri,length(Ri),1)
            # Term II
            tmp <- t(Bi) %*% Bi
            D <- D + weight[i] * kronecker(tmp,tmp)
            # Term III
            for(j in 1:Ni)
            {
                bij <- Bi[j,,drop=F]
                tmp <- t(bij) %*% bij
                D <- D - weight[i] * kronecker(tmp,tmp)
            }
            # Term IV
            for(j in 1:Ni)
            {
                bij <- Bi[j,,drop=F]
                E <- E + weight[i] * Ri[j,j] * kronecker(t(bij),t(bij))
            }
        }
    }

    # penalty term
    D <- D + lam * (kronecker(auxmat$W,auxmat$U) + kronecker(auxmat$V,auxmat$V))
    D <- 2*D
    E <- 2*E

    ## end of compution of DE

    D0 <- matrix(0,p*(p+1)/2,p*p)

    k <- 1
    for(ii in 1:p)
        for(jj in 1:ii)
        {
            k1 <- p * (jj-1) + ii;
            k2 <- p * (ii-1) + jj;

            if(ii != jj)
                D0[k,] <- D[k1,] + D[k2,]
            else
                D0[k,] <- D[k1,]

            k <- k + 1
        }

    D1 <- matrix(0,p*(p+1)/2,p*(p+1)/2)
    E1 <- matrix(p*(p+1)/2,1)
    k <- 1

    for (ii in 1:p)
        for (jj in 1:ii)
        {
            k1 = p * (jj-1) + ii
            k2 = p * (ii-1) + jj

            if (ii != jj)
            {
                D1[,k] <- D0[,k1] + D0[,k2]
                E1[k] <- E[k1] + E[k2]
            }
            else
            {
                D1[,k] <- D0[,k1]
                E1[k]  <- E[k1]
            }

            k <- k + 1
        }

    C0 <- solve(D1,-E1)
    C0 <- vec.to.lt(C0,p)
    C0 <- C0 + t(C0)
    diag(C0) <- diag(C0)/2

    return(list(C0=C0,D=D1,E=E1))


}

minimizeQ <- function(Lt,Lr,auxmat,lam,weight,L0=NULL,maxIt=10,verbose=F)
{
    p <- dim(auxmat$U)[1]
    if(length(lam) == 1)
    {
        if(is.null(L0))
        {
            C0 <- minimizeL(Lt,Lr,auxmat,lam,weight)$C0
            eigval <- eigen(C0)$values
            rho <- min(eigval)
            tmp <- Chat + abs(rho*1.1)*diag(rep(1,dim(Chat)[1]))
            L0 <- t(chol(tmp))
        }
    }

    L <- L0
    for(k in 1:maxIt)
    {
        if(verbose) message(sprintf('%d/%d ..done!',k,maxIt))
        H <- RhessianQ(Lt,Lr,matrix(0,p,p),auxmat$B,auxmat$U,auxmat$V,auxmat$W,lam,weight,L)$H
        #H <- hess.Q(Lt,Lr,matrix(0,p,p),auxmat,lam,weight,L)
        g <- RgradQ(Lt,Lr,L,auxmat$B,auxmat$U,auxmat$V,auxmat$W,lam,weight)
        #g <- grad.Q(Lt,Lr,L,auxmat,lam,weight)
        g <- lt.to.vec(g)
        eta <- solve(H,-g)
        eta <- vec.to.lt(eta,p)
        L <- expM(L,eta)
    }

    return(list(L=L,eta=eta,H=H))
}

# minimizeQ <- function(Lt,Lr,auxmat,lam,weight,L0=NULL,maxIt=10,verbose=F)
# {
#     p <- dim(auxmat$U)[1]
#     if(length(lam) == 1)
#     {
#         if(is.null(L0))
#         {
#             C0 <- minimizeL(Lt,Lr,auxmat,lam,weight)$C0
#             eigval <- eigen(C0)$values
#             rho <- min(eigval)
#             tmp <- Chat + abs(rho*1.1)*diag(rep(1,dim(Chat)[1]))
#             L0 <- t(chol(tmp))
#         }
#     }
#
#     L <- L0
#     for(k in 1:maxIt)
#     {
#         if(verbose) message(sprintf('%d/%d ..done!',k,maxIt))
#         H <- hess.Q(Lt,Lr,matrix(0,p,p),auxmat,lam,weight,L)
#         g <- grad.Q(Lt,Lr,L,auxmat,lam,weight)
#         g <- lt.to.vec(g)
#         eta <- solve(H,-g)
#         eta <- vec.to.lt(eta,p)
#         L <- expM(L,eta)
#     }
#
#     return(list(L=L,eta=eta,H=H))
# }


# use gradQ in Rhessian.cpp for speedup
# grad.Q <- function(Lt,Lr,L,auxmat,lam,weight)
# {
#     require(expm)
#     require(matrixcalc)
#
#     if(is.null(weight)) weight <- get.weight(Lt)
#
#     n <- length(Lt)
#     p <- dim(L)[1]
#
#     C <- L %*% t(L)
#
#     B <- auxmat$B
#     U <- auxmat$U
#     W <- auxmat$W
#     V <- auxmat$V
#
#     gQ <- Reduce('+',lapply(1:n,function(i){
#
#         Ni <- length(Lt[[i]])
#         if(Ni > 1)
#         {
#             Ri <- Lr[[i]]
#             Bi <- B[[i]]
#
#             tmp <- -4*(t(Bi) %*% Ri %*% Bi %*% L)
#
#             tmp <- tmp + 4*(t(Bi) %*% Bi %*% C %*% t(Bi) %*% Bi %*% L)
#
#             tmp <- tmp - 4*Reduce('+', lapply(1:Ni,function(j){
#                 bij <- Bi[j,,drop=F]
#                 t(bij) %*% bij %*% C %*% t(bij) %*% bij %*% L
#             }))
#
#             tmp <- tmp + 4* Reduce('+', lapply(1:Ni,function(j){
#                 bij <- Bi[j,,drop=F]
#                 Ri[j,j] * t(bij) %*% bij %*% L
#             }))
#             return(weight[i] * tmp)
#         }
#         else return(0)
#     }))
#
#     gQ <- gQ + 2 * lam * (U %*% C %*% W %*% L + W %*% C %*% U %*% L) +
#         4 * lam * V %*% C %*% V %*% L
#     return(gQ)
# }

expM <- function(X,S)
{
    (X-diag(diag(X))) + (S-diag(diag(S))) + diag(diag(X)*exp(diag(S)/diag(X)))
}

# use Rhessian.cpp instead for speedup
# hess.Q <- function(Lt,Lr,S,auxmat,lam,weight,X)
# {
#     n <- length(Lt)
#     L <- expM(X,S)
#
#     p <- dim(L)[1]
#     C <- L %*% t(L)
#     d <- p*(p+1)/2
#     H <- matrix(0,d,d)
#
#     U <- auxmat$U
#     W <- auxmat$W
#     V <- auxmat$V
#
#     G <- grad.Q(Lt,Lr,L,auxmat,lam,weight)
#
#     for(i in 1:n)
#     {
#         Ni <- length(Lt[[i]])
#         if(Ni > 1)
#         {
#             Ri <- Lr[[i]]
#             Bi <- auxmat$B[[i]]
#
#             E1 <- t(Bi) %*% Ri %*% Bi
#             E2 <- t(Bi) %*% Bi
#
#             for(j1 in 1:p)
#             {
#                 for(j2 in 1:j1)
#                 {
#                     j <- j1*(j1-1)/2 + j2
#                     for(k1 in 1:p)
#                     {
#                         for(k2 in 1:k1)
#                         {
#                             k <- k1*(k1-1)/2 + k2
#
#                             if(j1==j2) coef <- exp(S[j1,j2]/X[j1,j2])
#                             else coef <- 1
#                             if(k1==k2) coef <- coef * exp(S[k1,k2]/X[k1,k2])
#                             coef <- coef * weight[i]
#
#                             if(j2==k2)
#                             {
#                                 # Term I
#                                 H[j,k] <- H[j,k] - coef*E1[j1,k1]
#
#                                 # Term II
#                                 tmp1 <- t(L[,k2]) %*% E2 %*% L[,j2]
#                                 tmp2 <- (E2[j1,] %*% L[,k2]) * (E2[k1,] %*% L[,j2])
#                                 tmp3 <- E2[j1,] %*% C %*% E2[,k1]
#                                 H[j,k] <- H[j,k] + coef *(E2[j1,k1]*tmp1 + tmp2 + tmp3)
#                             }
#                             else
#                             {
#                                 # Term I: zero, no action
#                                 # Term II
#                                 tmp1 <- t(L[,k2]) %*% E2 %*% L[,j2]
#                                 tmp2 <- (E2[j1,] %*% L[,k2]) * (E2[k1,] %*% L[,j2])
#                                 H[j,k] <- H[j,k] + coef*(E2[j1,k1]*tmp1 + tmp2)
#                             }
#
#                             # Term III and IV
#                             for(m in 1:Ni)
#                             {
#                                 Eim <- t(Bi[m,,drop=F]) %*% Bi[m,]
#                                 if(j2 == k2)
#                                 {
#                                     # Term III
#                                     tmp1 <- t(L[,k2]) %*% Eim %*% L[,j2]
#                                     tmp2 <- (Eim[j1,] %*% L[,k2])*(Eim[k1,] %*% L[,j2])
#                                     tmp3 <- Eim[j1,] %*% C %*% Eim[,k1]
#                                     H[j,k] <- H[j,k] - coef*(Eim[j1,k1]*tmp1 + tmp2 + tmp3)
#
#                                     # Term IV
#                                     H[j,k] <- H[j,k] + coef*Ri[m,m]*Eim[j1,k1]
#                                 }
#                                 else
#                                 {
#                                     # Term III
#                                     tmp1 <- t(L[,k2]) %*% Eim %*% L[,j2]
#                                     tmp2 <- (Eim[j1,] %*% L[,k2])*(Eim[k1,] %*% L[,j2])
#                                     H[j,k] <- H[j,k] - coef*(Eim[j1,k1]*tmp1 + tmp2)
#
#                                     # Term IV: zero, no action
#                                 }
#                             }
#                         }
#                     }
#                 }
#             }
#         }
#     }
#
#     H <- 4*H
#
#     for(j in 1:p)
#     {
#         H[j*(j+1)/2,j*(j+1)/2] <- H[j*(j+1)/2,j*(j+1)/2] + G[j,j]*exp(S[j,j]/X[j,j])/X[j,j]
#     }
#
#     # penalty term
#     LH <- list()
#     for(q in 1:length(lam))
#     {
#         tmp <- H
#         for(j1 in 1:p)
#         {
#             for(j2 in 1:j1)
#             {
#                 j <- j1*(j1-1)/2 + j2
#                 for(k1 in 1:p)
#                 {
#                     for(k2 in 1:k1)
#                     {
#                         k <- k1*(k1-1)/2 + k2
#
#                         if(j1==j2) coef <- exp(S[j1,j2]/X[j1,j2])
#                         else coef <- 1
#
#                         if(j2==k2)
#                         {
#                             # Term I
#                             tmp1 <- U[j1,k1] * (t(L[,k2]) %*% W %*% L[,j2])
#                             tmp2 <- (U[j1,] %*% L[,k2]) * (W[k1,] %*% L[,j2])
#                             tmp3 <- U[j1,] %*% C %*% W[,k1]
#                             tmp[j,k] <- tmp[j,k] + 2*lam*coef*(tmp1 + tmp2 + tmp3)
#
#                             # Term II
#                             tmp1 <- W[j1,k1] * (t(L[,k2]) %*% U %*% L[,j2])
#                             tmp2 <- (W[j1,] %*% L[,k2]) * (U[k1,] %*% L[,j2])
#                             tmp3 <- W[j1,] %*% C %*% U[,k1]
#                             tmp[j,k] <- tmp[j,k] + 2*lam*coef*(tmp1 + tmp2 + tmp3)
#
#                             # Term III
#                             tmp1 <- V[j1,k1] * (t(L[,k2]) %*% V %*% L[,j2])
#                             tmp2 <- (V[j1,] %*% L[,k2]) * (V[k1,] %*% L[,j2])
#                             tmp3 <- V[j1,] %*% C %*% V[,k1]
#                             tmp[j,k] <- tmp[j,k] + 4*lam*coef*(tmp1 + tmp2 + tmp3)
#                         }
#                         else
#                         {
#                             # Term I
#                             tmp1 <- U[j1,k1] * (t(L[,k2]) %*% W %*% L[,j2])
#                             tmp2 <- (U[j1,] %*% L[,k2]) * (W[k1,] %*% L[,j2])
#                             tmp[j,k] <- tmp[j,k] + 2*lam*coef*(tmp1 + tmp2)
#
#                             # Term II
#                             tmp1 <- W[j1,k1] * (t(L[,k2]) %*% U %*% L[,j2])
#                             tmp2 <- (W[j1,] %*% L[,k2]) * (U[k1,] %*% L[,j2])
#                             tmp[j,k] <- tmp[j,k] + 2*lam*coef*(tmp1 + tmp2)
#
#                             # Term III
#                             tmp1 <- V[j1,k1] * (t(L[,k2]) %*% V %*% L[,j2])
#                             tmp2 <- (V[j1,] %*% L[,k2]) * (V[k1,] %*% L[,j2])
#                             tmp[j,k] <- tmp[j,k] + 4*lam*coef*(tmp1 + tmp2)
#                         }
#                     }
#                 }
#             }
#         }
#         LH[[q]] <- tmp
#     }
#
#     if(length(lam) == 1) LH <- LH[[1]]
#
#     return(LH)
# }

vec.to.lt <- function(v,p)
{
    k <- 1
    L <- matrix(0,p,p)
    for(i in 1:p)
        for(j in 1:i)
        {
            L[i,j] <- v[k]
            k <- k + 1
        }

    return(L)
}

lt.to.vec <- function(L)
{
    p <- dim(L)[1]
    k <- 1
    v <- matrix(0,p*(p+1)/2)
    for(i in 1:p)
    {
        for(j in 1:i)
        {
            v[k] <- L[i,j]
            k <- k + 1
        }
    }
    return(v)
}
