
## project points in x onto the intersection of c^tx=1 and d^tx=1
## expects x in row major format
project.line <- function(x,c,d) {
    m <- length(c) # ==length(d)==ncol(x)
    lambda <- sapply(list(list(c=c,d=d),list(c=d,d=c)),simplify=FALSE,FUN=
                     function(coefs) with(coefs,
                                         2 * ( (1-x%*%d)%*%(t(c)%*%d) - (1-x%*%c)%*%(t(d)%*%d) ) / as.numeric((t(c)%*%d)^2-(t(c)%*%c)*(t(d)%*%d))
                                         ))
    lambda <- do.call(cbind,lambda)
    y <- x + lambda[,1]%*%t(c) / 2 + lambda[,2]%*%t(d) / 2
    return(y)
}


## grid of points on the probability simplex in R^d with k points on
## each edge
## todo: d==1 case
simplex.grid <- function(k,d) {
    if(d==1)return(1)
    ## choose d-1 bars out of k+d-1 stars and bars
    bars <- combn(k+d-1,d-1)
    tuples <- apply(bars,2,diff)
    tuples <- rbind(bars[1,],tuples) - 1
    tuples <- rbind(tuples,k-colSums(tuples))
    wts <- tuples / k
    return(wts)
}


#' Sample points on the intersection of the probability simplex with
#' $\\{c^tx=1\\}$.
#'
#' This sampling routine is relatively fast but non-uniform. See the
#' referenced manuscript for details.

#' @param n the number of points to sample
#' @param c the coefficient vector
#' @return a matrix of points in the intersection
#' @export

p.sampler.1 <- function(n,c,shape=1) {
    if(anyNA(c))return(NA)
    d <- c - 1
    d.p <- d[d>=0]; d.m <- -d[d<0]
    if(length(d.p)*length(d.m)==0)return(NA)
    us <- replicate(n, {
        u.p <- rbeta(length(d.p),shape,shape)
        u.m <- rbeta(length(d.m),shape,shape)
        u.m <- u.m * c(u.p%*%d.p / (u.m%*%d.m))
        u <- numeric(length(d))
        u[d>=0] <- u.p
        u[d<0] <- u.m
        u
    })
    us <- t(us) / rowSums(t(us))
}

#' Sample points on the intersection of the probability simplex with
#' $\\{c^tx=1\\}$.
#'
#' This sampling routine is relatively slow but allows for
#' deterministic or stochastic sampling.  See the referenced
#' manuscript for details.

#' @param n the number of points to sample
#' @param c the coefficient vector
#' @param concentration dirichlet parameter; if NULL then
#'     a deterministic sample is used
#' @return a matrix of points in the intersection
#' @export

p.sampler.2 <- function(n,c,concentration=1,tol=1e-6) {
    if(anyNA(c))return(NA)
    m <- length(c)
    A <- rbind(c,-diag(m),rep(1,m))
    b <- c(1,rep(0,m),1)
    vv <- vertexenum::enumerate.vertices(A=A,b=b,warn=TRUE)
    idx <- which(abs(vv%*%c-1 + rowSums(vv)-1) < tol)
    if(length(idx)<1)return(NA)
    vv <- vv[idx,,drop=FALSE]
    if(is.null(concentration)) {
        n.vertices <- function(n,k)gamma(k+n)/gamma(n)/gamma(k+1)
        f <- function(k)n.vertices(m,k)-n
        k <- uniroot(f,interval=c(0,1),extendInt='yes')$root
        k <- round(k)
        stopifnot(k>0)
        w <- t(simplex.grid(k=k,d=nrow(vv)))
    } else {        
        w <- matrix(rgamma(n*nrow(vv),concentration,1),ncol=nrow(vv))
        w <- w/rowSums(w)
    }
    u <- w%*%vv
}



p.sampler.singular <- function(n,c) {
    m <- length(c)
    zeros.idx <- which(abs(c)<.Machine$double.eps)
    if(length(zeros.idx)==0)return(NA)
    n.vertices <- function(n,k)gamma(k+n)/gamma(n)/gamma(k+1)
    f <- function(k)n.vertices(m,k)-n
    k <- uniroot(f,interval=c(0,1),extendInt='yes')$root
    k <- round(k)
    stopifnot(k>0)
    simplex <- t(simplex.grid(k=k,d=length(zeros.idx)))
    out <- matrix(0,nrow(simplex),m)
    out[,zeros.idx] <- simplex
    return(out)
}



## projection
p.sampler.3 <- function(n,c,tol=12,scale.factor=NULL) {
    m <- length(c)
    if(!is.null(scale.factor)) {
        c0 <- (0:(m-1))/(m-1)
        theta <- c0[-1]/c[-1]
        dist <- min(abs(min(c0)-theta),abs(max(c0)-theta))/(max(c0)-min(c0))
        multiplier <- 1/dist
        n <- round(n*multiplier*scale.factor)
        }
    
    n.vertices <- function(n,k)gamma(k+n)/gamma(n)/gamma(k+1)
    f <- function(k)n.vertices(m,k)-n
    k <- uniroot(f,interval=c(0,1),extendInt='yes')$root
    k <- round(k)
    stopifnot(k>0)
    p.simplex <- t(simplex.grid(k=k,d=m))

    y <- project.line(p.simplex,c,rep(1,m))
    keep.idx <- apply(y,1,function(y.i)sum(y.i<0 | y.i>1)==0)
    y.keep <- y[keep.idx,,drop=FALSE]
    ## print(dim(unique(round(y.keep,tol))))
    unique(round(y.keep,tol))
}


## projection onto a subsimplex
p.sampler.5 <- function(n,c,tol=12,scale.factor=NULL) {
    m <- length(c)
    if(!is.null(scale.factor)) {
        c0 <- (0:(m-1))/(m-1)
        theta <- c0[-1]/c[-1]
        dist <- min(abs(min(c0)-theta),abs(max(c0)-theta))/(max(c0)-min(c0))
        multiplier <- 1/dist
        n <- round(n*multiplier*scale.factor)
        }
    
    n.vertices <- function(n,k)gamma(k+n)/gamma(n)/gamma(k+1)
    f <- function(k)n.vertices(m-1,k)-n
    k <- uniroot(f,interval=c(0,1),extendInt='yes')$root
    k <- round(k)
    stopifnot(k>0)
    p.simplex <- t(simplex.grid(k=k,d=m-1))
    p.simplex <- cbind(0,p.simplex)

    y <- project.line(p.simplex,c,rep(1,m))
    keep.idx <- apply(y,1,function(y.i)sum(y.i<0 | y.i>1)==0)
    y.keep <- y[keep.idx,,drop=FALSE]
    unique(round(y.keep,tol))
}



#' Center and studentize observed multinomial proportions
#' 
#' @param p.obs the observed proportions
#' @param p the null proportions
#' @param c the coefficient vector applied to the proportions
#' @return a vector of test statistics
#' @export

test.stat.1 <- function(p.obs,p,c) {
    ## browser()
    ## cat('.')
    var.c <- t(c)%*%(diag(as.numeric(p.obs)) - p.obs%*%t(p.obs))%*%c
    abs(c%*%(p.obs-p)) / sqrt(var.c)
}

## centered but not studentized
test.stat.3 <- function(p.obs,p,c) {
    abs(c%*%(p.obs-p))
}



#' Center and studentize observed multinomial proportions with an
#' added regularization term
#' 
#' @param p.obs the observed proportions
#' @param p the null proportions
#' @param c the coefficient vector applied to the proportions
#' @return a vector of test statistics
#' @export

test.stat.5 <- function(p.obs,p,c,n,lambda=.1) {
    x <- p.obs*n
    m <- length(p)
    p.bar <- (x+1/m) / (n+1)
    Sigma.hat <- diag(as.numeric(p.bar)) - p.bar%*%t(p.bar)
    regularization <- t(p.obs-p)[-1]%*%solve(Sigma.hat[-1,-1])%*%(p.obs-p)[-1]
    var.c <- t(c)%*%Sigma.hat%*%c
    abs(c%*%(p.obs-p)) / sqrt(var.c) + lambda*regularization
}

p.val.mc <- function(p,p.obs,n,n.ref.samples,c,T) {
    if(!isTRUE(abs(sum(p)-1)<1e-8))return(NA) 
    p.star <- rmultinom(n.ref.samples,n,prob=p) / n
    T.star <- apply(p.star,2,T,p=p,c=c)
    T.hat <- as.numeric(T(p.obs,p,c))
    mean(T.star >= T.hat)
}


#' Confidence interval for a linear function of a multinomial
#' parameter
#' 
#' This routine computes p-values for a hypothesis test at a grid of
#' points in the intersection of the probability simplex in $R^m$ and
#' $\\{c^tx=\\theta\\}$. The object it returns, of class multinom.linear,
#' may be used to obtain confidence intervals for $\\theta=c^t p$,
#' where $p$ is the vector of multinomial parameters. See the
#' referenced manuscript for notation and further details on the
#' algorithm.
#'
#' test.stat is a function for computing a test statitsic. It must
#' accept parameters p.obs with the same semantics as p.obs in this
#' routine, p corresponding to the a null multinomial parameter, i.e.,
#' a probability, c and n as in this routine. Supplied are:
#'
#' test.stat.1 : centered and studentized
#' test.stat.5 : centered and studentized with a regularization term
#'
#' p.sampler is a routine for computing a set of points in the
#' intersection of the probability simpler with $\\{c^tx=1\\}$. The only
#' obligatory parameters are n, corresponding to the number of points,
#' and c. Supplied are:
#'
#' p.sampler.1: faster stochastic sampler, non-uniform
#' p.sampler.2: slower deterministic or stochastic sampler, using vertex enumeration



#' @param p.obs the observed proportions
#' @param c the coefficients of the linear combination, assumed nonnegative
#' @param n the multinomial count
#' @param theta the grid of theta values; if null then theta.resolution evenly spaced will be used
#' @param n.ref.samples the number of monte carlo samples for computing the empirical p-values
#' @param test.stat routine to compute a test statistic
#' @param p.sampler routine returning a set of points in $\\{c^tx=1\\}$
#' @param theta.resolution the numer of points theta, default is 50
#' @param p.resolution the number of points p at each theta
#' @param ... additional parameters passed to p.sampler
#' @return an object of class multinom.linear
#' @export
#' @examples

#' ## first put the airpollution data in a convenient format: n rows and
#' ## m columns of 0s and 1s
#' airpollution.split <- split(bild::airpollution,bild::airpollution$id)
#' wheeze <- lapply(airpollution.split, function(df)
#'     do.call(rbind,rep(list(df$wheeze),unique(df$counts))))
#' wheeze <- do.call(rbind,wheeze)

#' ## call multinom.linear to carry out the main computations
#' alpha <- .05
#' m <- ncol(wheeze)
#' c <- (0:m)/m
#' counts <- rowSums(wheeze)
#' x <- sapply(0:m, function(i)sum(counts==i))
#' n <- sum(x)
#' test.stat <- function(p.obs,p,c)test.stat.5(p.obs,p,c,n,lambda=.05)
#' ml <- multinom.linear(p.obs=x/n,c=c,n=n,theta=NULL,n.ref.samples=1e2,test.stat=test.stat,p.sampler=p.sampler.1,theta.resolution=200,p.resolution=20)

#' ## print the confidence interval
#' confint(ml)

#' ## plot the p-values used to generate the confidence interval
#' plot(ml)

multinom.linear <- function(p.obs,c,n,theta,n.ref.samples=1e2,test.stat,p.sampler,theta.resolution=50,p.resolution=NULL,...) {
    ## browser()
    if(is.null(theta)) theta <- seq(min(c),max(c),len=theta.resolution) else theta <- sort(unique(theta))
    p.by.theta <- lapply(theta, function(theta.i) {
        p.sampler(n=p.resolution,c=c/theta.i,...)
    })
    p.by.theta[[which(abs(theta)<.Machine$double.eps)]]  <-  p.sampler.singular(n=p.resolution,c=c,...)
    
    p.val.bins <- lapply(p.by.theta, function(p){
        if(anyNA(p))return(NaN)
        apply(p,1,p.val.mc,p.obs=p.obs,n=n,n.ref.samples=n.ref.samples,c=c,T=test.stat)
    })
    max.pvals <- sapply(p.val.bins,max)
    max.pvals[is.na(max.pvals)] <- -Inf

    out <- structure(list(theta=theta,p=p.by.theta,p.val.p=p.val.bins,p.val.theta=max.pvals,data.name=paste(deparse(substitute(p.obs))),p.obs=p.obs),class='multinom.linear')
    return(out)
}



#' Compute a confidence interval by test inversion given p-values on a
#' grid of points in the parameter space

#' @param ml an object of class multinom.linear
#' @param alpha level of the confidence interval
#' @return a matrix with 2 columns, corresponding to the starting and
#'     ending theta values that make up the confidence interval
#' @export
#' @examples
#' ## see the documentation for multinom.linear for an example

confint.multinom.linear <- function(ml,alpha=.05) {
    rle0 <- rle(ml$p.val.theta>alpha)
    accept.runs <- which(rle0$values)
    ends <- cumsum(rle0$lengths)[accept.runs]
    starts <- ends - rle0$lengths[accept.runs] + 1
    cbind(start=ml$theta[starts],end=ml$theta[ends])
}

#' Plot the p-values for the hypothesis test that is inverted to form
#' a confidence interval

#' @param ml an object of class multinom.linear
#' @param alpha level of the test; indicated in the plot by a horizontal line
#' @param ... additional parameters passed to plot
#' @return plot output
#' @export
#' @examples
#' ## see the documentation for multinom.linear for an example
#' 
plot.multinom.linear <- function(ml,alpha=.05,...){
    theta <- ml$theta
    p <- ml$p
    p.val.theta <- ml$p.val.theta
    plot(theta,p.val.theta,ylim=c(0,1),xlab=expression(theta),ylab='p-value',...)
    abline(h=alpha)
    }



#' Carry out a hypothesis test of the null that a particular value
#' theta generated the sample

#' @param ml an object of class multinom.linear
#' @param theta.null the null value of the parameter to be tested
#' @return an object of class htest excpt for "statistic" entry, since
#'     multinom.linear doesnt currently expose the test statistics,
#'     only the p-values
#' @export
#' @examples
#' ## see the documentation for multinom.linear for an example
#' 

test.multinom.linear <- function(ml,theta.null) {
    idx <- which.min(abs(ml$theta-c(theta.null)))
    p.val.null <- ml$p.val.theta[idx]
    
    out <- with(ml,
                list(null.value=c(theta=theta.null),alternative='two-sided',method='Monte carlo exact test',estimate=c(p=as.numeric(p.obs)),data.name=data.name,statistic=NA,parameters=c(),p.value=p.val.null)
                )
    class(out) <- 'htest'
    return(out)
}



