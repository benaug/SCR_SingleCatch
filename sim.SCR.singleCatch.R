e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.singleCatch <- function(N=NA,p0=NA,sigma=NA,X=NA,buff=NA){
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
  times.true <- times #store these
  #times less than 1 are latent captures (if traps didn't fill up)
  times[times>1] <- Inf
  
  #convert capture times to the observed capture history
  y <- array(0,dim=c(N,J,K))
  store.n.animals <- store.n.traps <- rep(0,K)
  y.order <- array(Inf,dim=c(N,J,K))
  for(k in 1:K){
    n.total.cap <- sum(times[,,k]<Inf)
    n.animals <- sum(1*((rowSums(times[,,k]<Inf))>0))
    n.traps <- sum(1*((colSums(times[,,k]<Inf))>0))
    store.n.animals[k] <- n.animals
    store.n.traps[k] <- n.traps
    order.idx <- 1
    while(n.animals>0&n.traps>0){
      this.cap <- which(times[,,k]==min(times[,,k]),arr.ind=TRUE)
      y[this.cap[1],this.cap[2],k] <- 1
      y.order[this.cap[1],this.cap[2],k] <- order.idx
      focal.lambda <- lambda[this.cap[1],this.cap[2]]
      other.lambdas <- lambda[which(times[,,k]<Inf)]
      other.lambdas <- other.lambdas[-which(other.lambdas==focal.lambda)]
      times[this.cap[1],,k] <- Inf
      times[,this.cap[2],k] <- Inf
      n.animals <- sum(1*((rowSums(times[,,k]<Inf))>0))
      n.traps <- sum(1*((colSums(times[,,k]<Inf))>0))
      order.idx <- order.idx + 1
    }
  }
  
  y.true <- 1*(times.true<=1) #all captures that would have happened if traps weren't single catch
  captured <- which(rowSums(y)>0) #discard uncaptured individuals
  y.obs <- y[captured,,] #observed data
  n.cap <- nrow(y.obs) #number captured
  y.order <- y.order[captured,,] #latent capture order info here
  
  tmp <- which(y.obs==1,arr.ind=TRUE)
  obs.j <- tmp[,2]
  obs.k <- tmp[,3]
  obs.i <- order <- rep(NA,length(obs.j))
  for(l in 1:length(obs.j)){
    obs.i[l] <- which(y.obs[,obs.j[l],obs.k[l]]==1)
    order[l] <- y.order[obs.i[l],obs.j[l],obs.k[l]]
  }
  
  n.obs.cells <- as.numeric(colSums(table(obs.j,obs.k)))
  obs.i2D <- obs.j2D <- order2D <- matrix(NA,max(n.obs.cells),5)
  for(k in 1:K){
    obs.i2D[1:n.obs.cells[k],k] <- as.double(obs.i[obs.k==k]) #must be double for custom update
    obs.j2D[1:n.obs.cells[k],k] <- as.double(obs.j[obs.k==k])
    order2D[1:n.obs.cells[k],k] <- order[obs.k==k]
  }
  n.obs.cells.max <- max(n.obs.cells)
  
  #plot data
  par(mfrow=c(1,1),ask=FALSE)
  plot(X,pch=4,xlim=xlim,ylim=ylim,xlab="X",ylab="Y",main="Observed (black) and latent (gold) captures")
  y.true2D <- apply(y.true,c(1,2),sum)
  for(i in 1:N){
    trapcaps <- which(y.true2D[i,]>0)
    points(s[i,1],s[i,2],pch=16)
    for(j in trapcaps){
      lines(x=c(s[i,1],X[j,1]),y=c(s[i,2],X[j,2]),lwd=2,col="goldenrod") #latent events
    }
  }
  y.obs2D <- apply(y.obs,c(1,2),sum)
  for(i in 1:n.cap){
    trapcaps <- which(y.obs2D[i,]>0)
    for(j in trapcaps){
      lines(x=c(s[captured[i],1],X[j,1]),y=c(s[captured[i],2],X[j,2]),lwd=2,col="black")
    }
  }
  
  return(list(y.obs=y.obs,obs.i=obs.i,obs.j=obs.j,obs.k=obs.k,obs.i2D=obs.i2D,
              obs.j2D=obs.j2D,n.obs.cells=n.obs.cells,
              y.true=y.true,order=order,order2D=order2D,n.cap=n.cap,X=X,xlim=xlim,ylim=ylim,K2D=K2D))
}