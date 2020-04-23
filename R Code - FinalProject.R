
#--------------------------------------------------------------------------#
# R code for ST 502 Project
# Authors Anna Tomkins, Ashwin Bhargava, Meredith Saunders                               
#--------------------------------------------------------------------------#
#Plot for Poisson
x<-0:10
plot(x, dpois(x, .1),type="o",main="Distribution of the Poisson for different values of lambda"
     ,ylab="dpois")
lines(x, dpois(x, .5),col=2,type="o")
lines(x, dpois(x, 2),col=3,type="o")
lines(x, dpois(x, 5),col=4,type="o")
lines(x, dpois(x, 10),col=5,type="o")
legend(6,.8,c("lambda=0.1","lambda=0.5","lambda=2","lambda=5","lambda=10"),lty=1,col=1:5)


pl.pois <- function (theta, x, PRS=c("Default0", "Default1", "Score-Balanced"), b=0) {
  #NEsted intersection
  sb.order <- function(X, lambda, d) {
    
    x <- X
    if (lambda == 0) return(structure(x + 1, names = x))
    if(missing(d)) d <- dpois(x, lambda)
    r <- integer(length(x))
    dd <- x - lambda
    ES <- SD <- SV <- 0
    for (i in 1:length(r)) {
      
      if (all(dd > 0) || all(dd < 0)) {
        
        r[i:length(r)] <- X[order(abs(dd))] + 1
        break
        
      }
      ddp <- dd
      ddp[dd < 0] <- max(dd)
      ddn <- dd
      ddn[dd > 0] <- min(dd)
      K <- c(which.min(ddp), which.max(ddn))
      V <- SV + ((X[K] - lambda)**2 - lambda) * d[K]
      if(any(V > 0)) {
        
        if(all(V > 0)) stop("Whoops: something's wrong!") else K[V > 0] <- K[V < 0]
        
      }
      z <- (ES + d[K] * dd[K]) / (SD + d[K])
      k <- K[which.min(abs(z))]
      ES <- ES + d[k] * dd[k]
      SD <- SD + d[k]
      SV <- SV + ((X[k] - lambda)**2 - lambda) * d[k]
      r[i] <- X[k] + 1
      X <- X[-k]
      dd <- dd[-k]
      d <- d[-k]
      
    }
    if (length(unique(r)) != length(x)) stop("Whoops: something's wrong!")
    names(r) <- x[r]
    return(r)
    
  }
  
  pl.default <- function (theta, x) {
    
    pl <- numeric(length(theta))
    for(i in seq_along(theta)) {
      
      pl.1 <- max(1 - 2 * pgamma(theta[i], x), 0)
      pl.2 <- max(2 * pgamma(theta[i], x + 1) - 1, 0)
      pl[i] <- 1 - pl.1 - pl.2
      
    }
    pl <- data.frame(theta=theta, pl=pl)
    attr(pl, "x") <- x
    return(pl)
    
  }
  
  default.order <- function(X, lambda, d) {
    
    p <- ppois(X, lambda = lambda)
    k <- which(p >= 0.5)[1]
    if(k >= length(p)) stop("Whoops: k >= length(X)")
    s <- (k+1):length(p)
    p[s] <- 1 - p[s-1]
    r <- order(p, decreasing=TRUE)
    names(r) <- X[r]
    return(r)
    
  }
  
  im.order <- if(PRS=="Score-Balanced") sb.order else
    if(PRS=="Default0") return(pl.default(theta, x)) else
      if(PRS=="Default1") default.order else
        stop("unknown PRS", PRS)
  
  pl <- numeric(length(theta))
  for(i in seq_along(theta)) {
    
    lambda <- theta[i]
    d <- dpois(X <- 0:(2*max(as.integer(qpois(0.999999, lambda)), x+1)), lambda)
    r <- im.order(X, lambda, d)
    names(d) <- X
    a <- d[r]
    a <- c("-1"=0,a)
    k <- which(names(a)==x)
    pl[i] <- 1 - (bel <- sum(a[1:(k-1)])) # cumsum(a)[k-1]
    # EB-SB effect: put all conflict mass on lambda=b
    # if lambda is on the boundary and there are any conflict cases
    if(lambda == b && ppois(X[r[1]] - 1, lambda) >= ppois(x, b)) pl[i] = 1.0
    
  }
 pl <- data.frame(theta=theta, pl=pl)
 attr(pl, "x") <- x
  if(length(theta) == 1) attr(pl, "r") <- as.integer(names(r))
  return(pl)
  
}



#Figure 3
#x=0
t<-seq(0.0001, 6, by = .5)
plot(pl.pois(t,0),xlim=c(0, 6))
lines(pl.pois(t,0,PRS="Default0"),col="gray",lty=2)
lines(pl.pois(t,0,PRS="Score-Balanced"),col="gray")
abline(h =0.1, col="lightgray",lty=2)
abline(v = 0, col="lightgray",lty=2)
lines(pl.pois(t,0),type="S")

#x=3
t<-seq(0, 11, by = .5)
plot(pl.pois(t,3))
lines(pl.pois(t,3,PRS="Default0"),col="gray",lty=2)
lines(pl.pois(t,3,PRS="Default1"),col="gray")
lines(pl.pois(t,3,PRS="Score-Balanced"),col="gray")
abline(h =0.1, col="lightgray",lty=2)
abline(v = 3, col="lightgray",lty=2)
lines(pl.pois(t,3),type="S")


#x=7
t<-seq(.5, 16, by = .5)
plot(pl.pois(t,7))
lines(pl.pois(t,7,PRS="Default0"),col="gray",lty=2)
lines(pl.pois(t,7,PRS="Default1"),col="gray")
lines(pl.pois(t,7,PRS="Score-Balanced"),col="gray")
abline(h =0.1, col="lightgray",lty=2)
abline(v = 7, col="lightgray",lty=2)
lines(pl.pois(t,7),type="S")


#x=10
t<-seq(.5, 20, by = .5)
plot(pl.pois(t,10))
lines(pl.pois(t,10,PRS="Default0"),col="gray",lty=2)
lines(pl.pois(t,10,PRS="Default1"),col="gray")
lines(pl.pois(t,10,PRS="Score-Balanced"),col="gray")
abline(h =0.1, col="lightgray",lty=2)
abline(v = 10, col="lightgray",lty=2)
lines(pl.pois(t,10),type="S")


############################
#Figure 2
############################
pl.pois <- function (theta, x, PRS=c("Default0", "Default1", "Score-Balanced"), b=0) {
  
  sb.order <- function(X, lambda, d) {
    
    x <- X
    if (lambda == 0) return(structure(x + 1, names = x))
    if(missing(d)) d <- dpois(x, lambda)
    r <- integer(length(x))
    dd <- x - lambda
    ES <- SD <- SV <- 0
    for (i in 1:length(r)) {
      
      if (all(dd > 0) || all(dd < 0)) {
        
        r[i:length(r)] <- X[order(abs(dd))] + 1
        break
        
      }
      ddp <- dd
      ddp[dd < 0] <- max(dd)
      ddn <- dd
      ddn[dd > 0] <- min(dd)
      K <- c(which.min(ddp), which.max(ddn))
      V <- SV + ((X[K] - lambda)**2 - lambda) * d[K]
      if(any(V > 0)) {
        
        if(all(V > 0)) stop("Whoops: something's wrong!") else K[V > 0] <- K[V < 0]
        
      }
      z <- (ES + d[K] * dd[K]) / (SD + d[K])
      k <- K[which.min(abs(z))]
      ES <- ES + d[k] * dd[k]
      SD <- SD + d[k]
      SV <- SV + ((X[k] - lambda)**2 - lambda) * d[k]
      r[i] <- X[k] + 1
      X <- X[-k]
      dd <- dd[-k]
      d <- d[-k]
      
    }
    if (length(unique(r)) != length(x)) stop("Whoops: something's wrong!")
    names(r) <- x[r]
    return(r)
    
  }
  
  pl.default <- function (theta, x) {
    
    pl <- numeric(length(theta))
    for(i in seq_along(theta)) {
      
      pl.1 <- max(1 - 2 * pgamma(theta[i], x), 0)
      pl.2 <- max(2 * pgamma(theta[i], x + 1) - 1, 0)
      pl[i] <- 1 - pl.1 - pl.2
      
    }
    return(pl)
    
  }
  
  default.order <- function(X, lambda, d) {
    
    p <- ppois(X, lambda = lambda)
    k <- which(p >= 0.5)[1]
    if(k >= length(p)) stop("Whoops: k >= length(X)")
    s <- (k+1):length(p)
    p[s] <- 1 - p[s-1]
    r <- order(p, decreasing=TRUE)
    names(r) <- X[r]
    return(r)
    
  }
  
  im.order <- if(PRS=="Score-Balanced"){
    sb.order
  } else{
    if(PRS=="Default0"){
      return(pl.default(theta, x))
    } else{
      if(PRS=="Default1"){ 
        default.order
      } else stop("unknown PRS", PRS)
    }
  }
  
  
  lambda <- theta
  d <- dpois(X <- 0:(2*max(as.integer(qpois(0.999999, lambda)), x+1)), lambda)
  r <- im.order(X, lambda, d)
  names(d) <- X
  a <- d[r]
  a <- c("-1"=0,a)
  k <- which(names(a)==x)
  pl <- 1 - (bel <- sum(a[1:(k-1)])) # cumsum(a)[k-1]
  # EB-SB effect: put all conflict mass on lambda=b
  # if lambda is on the boundary and there are any conflict cases
  if(lambda == b && ppois(X[r[1]] - 1, lambda) >= ppois(x, b)) pl = 1.0
  
  
  
  return(pl)  
}



pl.pois = Vectorize( pl.pois, vectorize.args='x')



#Figure 2 
x <- rpois(n=100000, lambda=7)

#plot1 theta = 4
theta <- 4
p1 <- pl.pois(theta,x,PRS="Score-Balanced")
p2 <- pl.pois(theta,x,PRS="Default1")
p3 <- pl.pois(theta,x,PRS="Default0")

p1 <- sort(p1)
p2 <- sort(p2)
p3 <- sort(p3)



plot(ecdf(p1),xlim=c(0,1),verticals = TRUE,cex.points = par("cex"),do.points=FALSE,col="gray",col.01line = NULL,main="theta=4",xlab="Plausibility",ylab="CDF")
lines(ecdf(p3),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL,col="gray")
lines(ecdf(p2),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL)
abline(a=0,b=1)
abline(h =0.1, col="gray",lty=2,lwd=1)
abline(v = 0.1, col="gray",lty=2,lwd=1)

#plot1 theta = 6
theta <- 6
p1 <- pl.pois(theta,x,PRS="Score-Balanced")
p2 <- pl.pois(theta,x,PRS="Default1")
p3 <- pl.pois(theta,x,PRS="Default0")

p1 <- sort(p1)
p2 <- sort(p2)
p3 <- sort(p3)



plot(ecdf(p1),xlim=c(0,1),verticals = TRUE,cex.points = par("cex"),do.points=FALSE,col="gray",col.01line = NULL,main="theta=6",xlab="Plausibility",ylab="CDF")
lines(ecdf(p3),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL,col="gray")
lines(ecdf(p2),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL)
abline(a=0,b=1)
abline(h =0.1, col="gray",lty=2,lwd=1)
abline(v = 0.1, col="gray",lty=2,lwd=1)


#plot1 theta = 7
theta <- 7
p1 <- pl.pois(theta,x,PRS="Score-Balanced")
p2 <- pl.pois(theta,x,PRS="Default1")
p3 <- pl.pois(theta,x,PRS="Default0")

p1 <- sort(p1)
p2 <- sort(p2)
p3 <- sort(p3)



plot(ecdf(p1),xlim=c(0,1),verticals = TRUE,cex.points = par("cex"),do.points=FALSE,col="gray",col.01line = NULL,main="theta=7",xlab="Plausibility",ylab="CDF")
lines(ecdf(p3),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL,col="gray")
lines(ecdf(p2),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL)
abline(a=0,b=1)
abline(h =0.1, col="gray",lty=2,lwd=1)
abline(v = 0.1, col="gray",lty=2,lwd=1)

#plot1 theta = 8
theta <- 8
p1 <- pl.pois(theta,x,PRS="Score-Balanced")
p2 <- pl.pois(theta,x,PRS="Default1")
p3 <- pl.pois(theta,x,PRS="Default0")

p1 <- sort(p1)
p2 <- sort(p2)
p3 <- sort(p3)



plot(ecdf(p1),xlim=c(0,1),verticals = TRUE,cex.points = par("cex"),do.points=FALSE,col="gray",col.01line = NULL,main="theta=8",xlab="Plausibility",ylab="CDF")
lines(ecdf(p3),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL,col="gray")
lines(ecdf(p2),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL)
abline(a=0,b=1)
abline(h =0.1, col="gray",lty=2,lwd=1)
abline(v = 0.1, col="gray",lty=2,lwd=1)

#plot1 theta = 10
theta <- 10
p1 <- pl.pois(theta,x,PRS="Score-Balanced")
p2 <- pl.pois(theta,x,PRS="Default1")
p3 <- pl.pois(theta,x,PRS="Default0")

p1 <- sort(p1)
p2 <- sort(p2)
p3 <- sort(p3)



plot(ecdf(p1),xlim=c(0,1),verticals = TRUE,cex.points = par("cex"),do.points=FALSE,col="gray",col.01line = NULL,main="theta=10",xlab="Plausibility",ylab="CDF")
lines(ecdf(p3),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL,col="gray")
lines(ecdf(p2),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL)
abline(a=0,b=1)
abline(h =0.1, col="gray",lty=2,lwd=1)
abline(v = 0.1, col="gray",lty=2,lwd=1)

#plot1 theta = 12
theta <- 12
p1 <- pl.pois(theta,x,PRS="Score-Balanced")
p2 <- pl.pois(theta,x,PRS="Default1")
p3 <- pl.pois(theta,x,PRS="Default0")

p1 <- sort(p1)
p2 <- sort(p2)
p3 <- sort(p3)



plot(ecdf(p1),xlim=c(0,1),verticals = TRUE,cex.points = par("cex"),do.points=FALSE,col="gray",col.01line = NULL,main="theta=12",xlab="Plausibility",ylab="CDF")
lines(ecdf(p3),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL,col="gray")
lines(ecdf(p2),xlim=c(0,1),verticals = TRUE,do.points=FALSE,col.01line = NULL)
abline(a=0,b=1)
abline(h =0.1, col="gray",lty=2,lwd=1)
abline(v = 0.1, col="gray",lty=2,lwd=1)





        