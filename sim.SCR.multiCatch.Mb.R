e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.multiCatch.Mb <- function(N=NA,p0.p=NA,p0.c=NA,sigma=NA,X=NA,buff=NA){
  xlim <- range(X[,1]) + c(-buff,buff)
  ylim <- range(X[,2]) + c(-buff,buff)
  s <- cbind(runif(N,xlim[1],xlim[2]),runif(N,ylim[1],ylim[2]))
  D <- e2dist(s,X)
  pd.p <- p0.p*exp(-D*D/(2*sigma*sigma))
  pd.c <- p0.c*exp(-D*D/(2*sigma*sigma))
  lambda.p <- -log(1-pd.p)
  lambda.c <- -log(1-pd.c)
  J <- nrow(X)
  K2D <- matrix(1,J,K) #trap operation. assuming perfect here
  #simulate capture times
  times.p <- times.c <- array(Inf,dim=c(N,J,K))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:K){
        if(lambda.p[i,j]>0&K2D[j,k]>0){
          times.p[i,j,k] <- rexp(1,lambda.p[i,j])
        }
        if(lambda.c[i,j]>0&K2D[j,k]>0){
          times.c[i,j,k] <- rexp(1,lambda.c[i,j])
        }
      }
    }
  }
  times.use <- times.p
  y <- y.state <- array(0,dim=c(N,J,K))
  for(i in 1:N){
    for(k in 1:K){
      if(any(times.use[i,,k]<1)){
        j.cap <- which(times.use[i,,k]==min(times.use[i,,k]))
        y[i,j.cap,k] <- 1
        if(k<K){
          y.state[i,j.cap,(k+1):K] <- 1
          times.use[i,j.cap,(k+1):K] <- times.c[i,j.cap,(k+1):K]
        }
      }
    }
  }
  y.full <- y
  y.state.full <- y.state
  captured <- which(rowSums(y)>0) #discard uncaptured individuals
  y <- y[captured,,] #observed data
  y.state <- y.state[captured,,]
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
  return(list(y=y,y.state=y.state,n.cap=n.cap,X=X,xlim=xlim,ylim=ylim,K2D=K2D))
}