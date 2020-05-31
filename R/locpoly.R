#' one-dimensional local polynomial smoother
#' @param x a vector of observed values of the predictor
#' @param y a vector of observed values of the response
#' @param h the bandwidth, either a numeric or NULL; if NULL, the automatically selected by five-fold CV
#' @param newx a vector of new values of the predictor; set to \code{x} by default
#' @param degree the degree of local polynomials, nonnegative integer
#' @param weight a vector of real numbers obtaining the weight for each observations; if NULL, then equal weight
#' @param kernel a kernel from the package locpol, supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic", "sigmoid" and "silverman"
#' @param sorted indicator of whether the array x and newx are sorted in increasing order 
#' @return a list of the following elements
#' \item{h}{the selected bandwidth or input bandwidth}
#' \item{x}{the original input}
#' \item{y}{the original input}
#' \item{weight}{the original input}
#' \item{kernel}{the original input}
#' \item{degree}{the original input}
#' \item{yhat}{the estimated response at \code{newx}}
#' @example 
#' x <- regular.grid()
#' y <- sin(2*pi*x) + 0.1*rnorm(length(x))
#' yhat <- lp1D(x,y)
#' @export
lp1D <- function(x,y,h=NULL,newx=x,degree=1,weight=NULL,kernel='epanechnikov',sorted=F)
{
    if(is.null(weight)) weight <- rep(1,length(y))
    
    if(is.null(h)) h <- bw.lp1D(x,y,weight,kernel,degree,method='cv')
    
    if(sorted==F)
    {
        ord1 <- sort(x,index.return=T)$ix
        x <- x[ord1]
        y <- y[ord1]
        
        
        ord2 <- sort(newx,index.return=T)$ix
        
        yhat <- rep(0,length(newx))
        
        yhat[ord2] <- csmoothmean(x=x,
                                  z=y,
                                  w=weight,
                                  h=h,
                                  kernel=kernel,
                                  d=degree,
                                  newx=newx[ord2])
    }
    else
    {
        yhat <- csmoothmean(x=x,
                            z=y,
                            w=weight,
                            h=h,
                            kernel=kernel,
                            d=degree,
                            newx=newx)
    }
    
    list(yhat=yhat,x=x,y=y,kernel=kernel,h=h,degree=degree,weight=weight)
}



#' select bandwidth for one-dimensional local linear smoother
#' @param x a vector of observed values of the predictor
#' @param y a vector of observed values of the response
#' @param weight a vector of real numbers obtaining the weight for each observations
#' @param kernel a kernel from the package locpol
#' @param degree the degree of local polynomials
#' @param method selection method, supported is 'cv'
#' @param K the number of folds, a parameter for CV method
#' @param H candidate values of the bandwidth; if NULL, then automathcally generated
#' @return the selected bandwidth
#' @example 
#' x <- regular.grid()
#' y <- sin(2*pi*x) + 0.1*rnorm(length(x))
#' h <- bw.lp1D(x,y,rep(1,length(x)),kernel='epanechnikov',degree=1)
#' lp1D(x,y,h=h)
#' @export
bw.lp1D <- function(x,y,weight,kernel,degree,method='cv',K=5,H=NULL)
{
    bw <- NULL
    if(method=='cv')
    {
        a <- min(x)
        b <- max(x)
        
        if(is.null(H))  H <- 10^seq(-2,0,length.out=20) * (b-a)/3
        #max.it <- 5
        #it <- 0
        n <- length(y)
        
        #while(it < max.it)
        #{
        E <- matrix(0,length(H),K)
        tmp <- cv.partition(n,K)
        for(j in 1:K)
        {
            tridx <- tmp$training[[j]]
            teidx <- tmp$test[[j]]
            for(i in 1:length(H))
            {
                yhat <- lp1D (x[tridx],y[tridx],H[i],newx=x[teidx],degree=degree,weight=weight,kernel=kernel)$yhat
                E[i,j] <- E[i,j] + sum((yhat-y[teidx])^2)
            }
        }
        
        E <- apply(E,1,sum)
        bw <- H[which.min(E)]
        
        #}
    }
    else stop('bw must be cv or plug.in or thumb')
    
    return(bw)
}
