get.weig.mean <- function(t,y,scheme=c('OBS','SUBJ'))
{
    scheme <- toupper(scheme)
    scheme <- match.arg(scheme)

    weig <- NULL

    if(scheme == 'OBS')
    {
        if(is.list(t))
        {
            mi <- sapply(t,length)
            weig <- rep(1/sum(mi),length(t))
        }
        else
        {
            n <- nrow(y)
            m <- ncol(y)
            weig <- rep(1/(n*m),n)
        }
    }
    else if(scheme == 'SUBJ')
    {
        if(is.list(t))
        {
            n <- length(t)
            weig <- sapply(t,function(s) 1/(n*length(s)))
        }
        else
        {
            n <- nrow(y)
            m <- ncol(y)
            weig <- rep(1/(n*m),n)
        }
    }
    else stop(paste0('scheme ',scheme, ' is not recognized.'))

    return(weig)
}

get.weig.cov <- function(t,y,scheme=c('OBS','SUBJ'))
{
    scheme <- toupper(scheme)
    scheme <- match.arg(scheme)

    weig <- NULL

    if(scheme == 'OBS')
    {
        if(is.list(t))
        {
            mi <- sapply(t,length)
            m <- sum(mi*(mi-1))
            weig <- rep(1/m,length(t))
        }
        else
        {
            n <- nrow(y)
            m <- ncol(y)
            weig <- rep(1/(n*(m*(m-1))),n)
        }
    }
    else if(scheme == 'SUBJ')
    {
        if(is.list(t))
        {
            n <- length(t)
            weig <- sapply(t,function(s) 1/(n*(length(s)*(length(s)-1))))
        }
        else
        {
            n <- nrow(y)
            m <- ncol(y)
            weig <- rep(1/(n*(m*(m-1))),n)
        }
    }
    else stop(paste0('scheme ',scheme, ' is not recognized.'))

    return(weig)
}

