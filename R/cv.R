#' the function mimic the cvpartition in MATLAB
#' @param n sample size
#' @param K the number of cv folders
#' @return a list of the follow elements
#' \item{training}{a list of \code{K} elements. Each element is a vector of indices of samples for training}
#' \item{test}{a list of \code{K} elements. Each element is a vector of indices of samples for test}
#' \item{num.test.sets}{an integer signifying the number of test sets}
#' @export cv.partition
cv.partition <- function(n,K)
{
    if(K <= 1) stop('K must be larger than 1')

    I <- sample.int(n)
    a <- n %% K
    b <- (n-a)/K

    trunk.size <- c(rep(b+1,a),rep(b,K-a))
    start.idx <- c(1,1+cumsum(trunk.size[1:(K-1)]))
    end.idx <- start.idx + trunk.size - 1

    test <- lapply(1:K,function(k){
        I[start.idx[k]:end.idx[k]]
    })

    training <- lapply(1:K,function(k){
        I[-(start.idx[k]:end.idx[k])]
    })

    return(list(training=training,test=test,num.test.sets=length(test)))
}
