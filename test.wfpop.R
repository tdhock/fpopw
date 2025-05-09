library(Rcpp)
library(devtools)

## compile Rcpp, build and installexample
## registration if not already done
#tools::package_native_routine_registration_skeleton('fpopw',con="fpopw/src/init.c",,FALSE)
compileAttributes("fpopw")
document("fpopw")
check("fpopw")
build("fpopw")
install("fpopw")

library(fpopw)
## test selection and smt using Fpop and Fpsn
## Fpop
n <- 1000
x0 <- rnorm(n)
res <- Fpop(x0, 2*log(n))
smt <- getSMT(res)

## Fpop_w
n <- 1000
x0 <- rnorm(n)
w0 <- 0.3 + runif(n)
res <- Fpop_w(x0, w0, 2*log(n))
smt <- getSMT(res)
weighted.mean(x0, w0) == smt[1] ## if no change should be equal

## Fpsn
x0 <- rnorm(1000) + rep(c(0, 2, -1), c(100, 400, 500))
res <- Fpsn(x0, 100)

## model selection
K.var <- select_Fpsn(res, method="givenVariance") 
K.biggest <- select_Fpsn(res, method="biggest.S3IB") 
K.notbiggest <- select_Fpsn(res, method="notbiggest.S3IB") 

## get smt line
smt <- getSMT(res, K.var)
plot(x0, pch=20, col="blue"); lines(smt, col="red", lwd=2)
##


##############################################################
## check fpsn and fpsn no mem give the same costs
##############################################################
n <- 3000; pri <- 0.1

identique <- TRUE
i.rep <- 1

## checking that we get the same changepoint vectors as gfpop
while(identique & i.rep < 100){
 if(i.rep %% 10) cat(".") else cat("|")
 Y <- rnorm(n)
 Weights <- runif(n)/(1-pri) + pri
 Kmax.sn <- 100
 system.time(res.sn <- Fpsn_w(Y, Weights, Kmax=Kmax.sn))
 system.time(res.sn_nomem <- Fpsn_w_nomemory(Y, Weights, Kmax=Kmax.sn))
 identique <- all.equal(res.sn_nomem$J.est, res.sn$J.est)
 i.rep <- i.rep + 1
}


##############################################################
## check same as gfpop on a few simulated example
##############################################################
library(gfpop)
getWaitGraph <- function(min.len, penalty, K=Inf){
  emptyGraph <- graph()
  if(min.len == 1){
    return(graph(penalty = penalty, type = "std", K=K)) ## classic graph
  }
  
  if(min.len >= 2){
    emptyGraph <- graph()
    emptyGraph <- graph(emptyGraph, Edge(1, 1, "null", K=K))
    emptyGraph <- graph(emptyGraph, Edge(1, 2, penalty=penalty, "std", K=K))
    if(min.len >=3){
    for(i in 3:(min.len)){
      emptyGraph <- graph(emptyGraph, Edge((i-1), i, "null", K=K))
    }}
    emptyGraph <- graph(emptyGraph, Edge(min.len, 1, "null", K=K))
    return(emptyGraph)
  }
  
}



n <- 10^5
identique <- TRUE
i.rep <- 1

## checking that we get the same changepoint vectors as gfpop
while(identique & i.rep < 100){
 if(i.rep %% 10) cat(".") else cat("|")
 x <- rnorm(n) + rep(c(0, 3, 2, -1), each=n/4)
 w <- 10^-3 + runif(n)/(1-10^-3)

 ans <- gfpop(data = x,  weights=w,
                   mygraph = getWaitGraph(min.len=1, penalty = 1*log(n)), 
                   type = "mean")


 res <- fpopw::Fpop_w(x, w, 1*log(n), mini = min(x), maxi = max(x))

 identique <- all.equal(ans$changepoints, res$t.est)
 i.rep <- i.rep +1 
}

##



