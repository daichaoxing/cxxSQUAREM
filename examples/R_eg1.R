#Before Running
#1. Install Rcpp package and SQUAREM package on R
#2. Put SQUAREM.cpp file in the current working directory
#3. Run

clc <- function() cat(rep("\n",50)) #clear console
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}
toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}


clc()
set.seed(2015)
library(Rcpp)
library(SQUAREM)

#Example from original package SQUAREM
poissmix.em <- function(p,y) {
# The fixed point mapping giving a single E and M step of the EM algorithm

pnew <- rep(NA,3)
i <- 0:(length(y)-1)
zi <- p[1]*exp(-p[2])*p[2]^i / (p[1]*exp(-p[2])*p[2]^i + (1 - p[1])*exp(-p[3])*p[3]^i)
pnew[1] <- sum(y*zi)/sum(y)
pnew[2] <- sum(y*i*zi)/sum(y*zi)
pnew[3] <- sum(y*i*(1-zi))/sum(y*(1-zi))
p <- pnew
return(pnew)
}



poissmix.loglik <- function(p,y) {
# Objective function whose local minimum is a fixed point 
# negative log-likelihood of binary poisson mixture
i <- 0:(length(y)-1)
loglik <- y*log(p[1]*exp(-p[2])*p[2]^i/exp(lgamma(i+1)) + 
		(1 - p[1])*exp(-p[3])*p[3]^i/exp(lgamma(i+1)))
return ( -sum(loglik) )
}

poissmix.dat <- data.frame(death=0:9, freq=c(162,267,271,185,111,61,27,8,3,1))
y <- poissmix.dat$freq
tol <- 1.e-08

p0 <- c(runif(1),runif(2,0,4))  # random starting value
tic()
pf1 <- fpiter(p=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, control=list(tol=tol))
toc()


tic()
pf2 <- squarem(par=p0, y=y, fixptfn=poissmix.em, objfn=poissmix.loglik, 
               control=list(tol=tol,trace=TRUE))
toc()

sourceCpp("R_SQUAREM.cpp")
tic()
pf3 <- cxxSQUAREM(p0,y,fixptfn=poissmix.em,objfn=poissmix.loglik)           
toc()
all.equal(pf2,pf3)

