

#' one-dimensional local linear smoother
#' @param x a vector of observed values of the predictor
#' @param y a vector of observed values of the response
#' @param newx a vector of new values of the predictor
#' @param bw the bandwidth, either a numeric or an element from \code{c('cv','plug.in','thumb')}
#' @param weig a vector of real numbers obtaining the weight for each observations
#' @param kernel a kernel from the package locpol
#' @param deg the degree of local polynomials
#' @return a list of the following elements
#' \item{bw}{the selected bandwidth or input bandwidth}
#' \item{x}{the original input}
#' \item{y}{the original input}
#' \item{weig}{the original input}
#' \item{kernel}{the original input}
#' \item{deg}{the original input}
loclin1D <- function(x,y,newx,bw='cv',
                            weig=rep(1,length(y)),
                            kernel=EpaK,deg=1)
{
    if(!is.numeric(bw))
        bw <- compute.bw.1D(x,y,bw,weig,kernel,deg)

    fitted <- locpol::locPolSmootherC(x, y, newx, bw, deg,kernel, weig=weig)$beta0
    while(any(is.na(fitted) | is.infinite(fitted)))
    {
        warning('bw is two small. increase it by multiplier 1.2')
        bw <- bw * 1.2
        fitted <- locpol::locPolSmootherC(x, y, newx, bw, deg,kernel, weig=weig)$beta0
    }
    return(list(bw=bw,fitted=fitted,x=x,y=y,weig=weig,kernel=kernel,deg=deg))
}

#' select bandwidth for one-dimensional local linear smoother
#' @param x a vector of observed values of the predictor
#' @param y a vector of observed values of the response
#' @param method selection method
#' @param weig a vector of real numbers obtaining the weight for each observations
#' @param kernel a kernel from the package locpol
#' @param deg the degree of local polynomials
#' @return the selected bandwidth
#' @keywords internal
compute.bw.1D <- function(x,y,method=c('cv','plug.in','thumb'),weig,kernel,deg)
{
    bw <- NULL
    if(method=='cv')
    {
        piBwSel <- locpol::pluginBw(x, y, deg, kernel,weig=weig)
        thBwSel <- locpol::thumbBw(x, y, deg, kernel,weig=weig)
        #interval <- c(min(piBwSel/2,thBwSel/2),max(thBwSel*2,piBwSel*2))
        #interval <- c(piBwSel/2,piBwSel*2)
        interval <- c(thBwSel/10,thBwSel*20)
        max.it <- 10
        it <- 0
        best.bw <- piBwSel
        while(it < max.it)
        {
            bw <- locpol::regCVBwSelC(x, y, deg, kernel, weig=weig, interval=interval)
            if(bw > interval[1] && bw < interval[2]) break
            else if (bw == best.bw) break
            else if (bw < interval[1]){
                interval[2] <- interval[1]
                interval[1] <- interval[1] / 2
            }
            else
            {
                interval[1] <- interval[2]
                interval[2] <- interval[2] * 2
            }
            it <- it + 1
        }
    }
    else if(method=='plug.in')
    {
        bw <- locpol::pluginBw(x, y, deg, kernel,weig=weig)
    }
    else if(method=='thumb')
    {
        bw <- locpol::thumbBw(x, y, 1, kernel,weig=weig)
    }
    else stop('bw must be cv or plug.in or thumb')

    return(bw)
}
