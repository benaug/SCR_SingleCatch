e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.multiCatch <- function(N=NA,p0=NA,sigma=NA,X=NA,buff=NA){
  xlim <- range(X[,1]) + c(-buff,buff)
  ylim <- range(X[,2]) + c(-buff,buff)
  s <- cbind(runif(N,xlim[1],xlim[2]),runif(N,ylim[1],ylim[2]))
  D <- e2dist(s,X)
  pd <- p0*exp(-D*D/(2*sigma*sigma))
  lambda <- -log(1-pd)
  J <- nrow(X)
  K2D <- matrix(1,J,K) #trap operation. assuming perfect here
  #simulate capture times
  times <- array(Inf,dim=c(N,J,K))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:K){
        if(lambda[i,j]>0&K2D[j,k]>0){
          times[i,j,k] <- rexp(1,lambda[i,j])
        }
      }
    }
  }
  y <- array(0,dim=c(N,J,K))
  for(i in 1:N){
    for(k in 1:K){
      if(any(times[i,,k]<1)){
        j.cap <- which(times[i,,k]==min(times[i,,k]))
        y[i,j.cap,k] <- 1       
      }
    }
  }
  y.full <- y
  captured <- which(rowSums(y)>0) #discard uncaptured individuals
  y <- y[captured,,] #observed data
  n.cap <- nrow(y) #number captured
  
  #plot data
  par(mfrow=c(1,1),ask=FALSE)
  plot(X,pch=4,xlim=xlim,ylim=ylim,xlab="X",ylab="Y")
  y2D <- apply(y.full,c(1,2),sum)
  for(i in 1:N){
    trapcaps <- which(y2D[i,]>0)
    points(s[i,1],s[i,2],pch=16)
    for(j in trapcaps){
      lines(x=c(s[i,1],X[j,1]),y=c(s[i,2],X[j,2]),lwd=2,col="black")
    }
  }
  return(list(y=y,n.cap=n.cap,X=X,xlim=xlim,ylim=ylim,K2D=K2D))
}