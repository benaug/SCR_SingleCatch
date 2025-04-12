e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SCR.singleCatch.Dcov.Mb <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,xlim=NA,ylim=NA,res=NA,
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
  lambda.p <- -log(1-pd.p)
  pd.c <- p0.c*exp(-D*D/(2*sigma*sigma))
  lambda.c <- -log(1-pd.c)
  J <- nrow(X)
  K2D <- matrix(1,J,K) #trap operation. assuming perfect here
  #simulate capture times for first and subsequent captures
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
  times.true.p <- times.p #store these
  times.true.c <- times.c #store these
  #times less than 1 are latent captures (if traps didn't fill up)
  times.p[times.p>1] <- Inf
  times.c[times.c>1] <- Inf
  times.use <- times.p
  #convert capture times to the observed capture history
  y <- y.state <- array(0,dim=c(N,J,K))
  store.n.animals <- store.n.traps <- rep(0,K)
  y.order <- array(Inf,dim=c(N,J,K))
  for(k in 1:K){
    n.total.cap <- sum(times.use[,,k]<Inf)
    n.animals <- sum(1*((rowSums(times.use[,,k]<Inf))>0))
    n.traps <- sum(1*((colSums(times.use[,,k]<Inf))>0))
    store.n.animals[k] <- n.animals
    store.n.traps[k] <- n.traps
    order.idx <- 1
    while(n.animals>0&n.traps>0){
      this.cap <- which(times.use[,,k]==min(times.use[,,k]),arr.ind=TRUE)
      y[this.cap[1],this.cap[2],k] <- 1
      y.order[this.cap[1],this.cap[2],k] <- order.idx
      times.use[this.cap[1],,k] <- Inf
      times.use[,this.cap[2],k] <- Inf
      n.animals <- sum(1*((rowSums(times.use[,,k]<Inf))>0))
      n.traps <- sum(1*((colSums(times.use[,,k]<Inf))>0))
      order.idx <- order.idx + 1
    }
    if(k<K){
      #update y.state
      for(i in 1:N){
        for(j in 1:J){
          if((y[i,j,k]>0)){ #if captured, update y.state and capture times to subsequent
            y.state[i,j,(k+1):K] <- 1
            times.use[i,j,(k+1):K] <- times.c[i,j,(k+1):K]
          }
        }
      }
    }
  }
  
  times.true <- times.true.p
  times.true[y.state==1] <- times.true.c[y.state==1]
  
  y.true <- 1*(times.true<=1) #all captures that would have happened if traps weren't single catch
  captured <- which(rowSums(y)>0) #discard uncaptured individuals
  y.obs <- y[captured,,] #observed data
  y.state <- y.state[captured,,]
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
  image(x.vals,y.vals,matrix(lambda.cell*InSS,n.cells.x,n.cells.y),xlab="X",ylab="Y",
        main="Observed (black) and latent (gold) captures",col=cols1)
  points(X,pch=4,lwd=2)
  y.true2D <- apply(y.true,c(1,2),sum)
  for(i in 1:N){
    trapcaps <- which(y.true2D[i,]>0)
    for(j in trapcaps){
      lines(x=c(s[i,1],X[j,1]),y=c(s[i,2],X[j,2]),lwd=2,col="goldenrod") #latent events
    }
    points(s[i,1],s[i,2],pch=16)
  }
  y.obs2D <- apply(y.obs,c(1,2),sum)
  for(i in 1:n.cap){
    trapcaps <- which(y.obs2D[i,]>0)
    for(j in trapcaps){
      lines(x=c(s[captured[i],1],X[j,1]),y=c(s[captured[i],2],X[j,2]),lwd=2,col="black")
    }
  }
  
  return(list(y.obs=y.obs,y.state=y.state,obs.i=obs.i,obs.j=obs.j,obs.k=obs.k,obs.i2D=obs.i2D,
              obs.j2D=obs.j2D,n.obs.cells=n.obs.cells,
              y.true=y.true,order=order,order2D=order2D,n.cap=n.cap,X=X,K2D=K2D,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N))
}