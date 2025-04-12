e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.multiCatch.Dcov.Mb <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,xlim=NA,ylim=NA,res=NA,
                                       p0.p=NA,p0.c=NA,sigma=NA,X=NA){
  #get expected N
  cellArea <- res^2
  lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
  lambda.N <- sum(lambda.cell)
  #simulate realized N
  N <- rpois(1,lambda.N)
  
  #recreate some Dcov things so we can pass fewer arguments into this function
  x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
  y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
  dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
  cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
  n.cells <- nrow(dSS)
  n.cells.x <- length(x.vals)
  n.cells.y <- length(y.vals)
  
  # simulate a population of activity centers
  pi.cell <- lambda.cell/sum(lambda.cell)
  #zero out non-habitat
  pi.cell[InSS==0] <- 0
  s.cell <- sample(1:n.cells,N,prob=pi.cell,replace=TRUE)
  #distribute activity centers uniformly inside cells
  s <- matrix(NA,nrow=N,ncol=2)
  for(i in 1:N){
    tmp <- which(cells==s.cell[i],arr.ind=TRUE) #x and y number
    s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]+res/2])
    s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]+res/2])
  }
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
  #plot data
  par(mfrow=c(1,1),ask=FALSE)
  image(x.vals,y.vals,matrix(lambda.cell*InSS,n.cells.x,n.cells.y),xlab="X",ylab="Y",col=cols1)
  y2D <- apply(y.full,c(1,2),sum)
  for(i in 1:N){
    trapcaps <- which(y2D[i,]>0)
    points(s[i,1],s[i,2],pch=16)
    for(j in trapcaps){
      lines(x=c(s[i,1],X[j,1]),y=c(s[i,2],X[j,2]),lwd=2,col="black")
    }
  }
  points(X,pch=4)
  return(list(y=y,y.state=y.state,n.cap=n.cap,X=X,K2D=K2D,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N))
}