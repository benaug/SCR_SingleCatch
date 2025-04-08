NimModel <- nimbleCode({
  p0.p ~ dunif(0,1) #first capture p0
  p0.c ~ dunif(0,1) #subsequent capture p0
  sigma ~ dunif(0,100)
  lambda.N ~ dunif(0,100) #expected N
  N ~ dpois(lambda.N) #realized N
  for(i in 1:M){ #N/z and y.true update under the hood
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    kern[i,1:J] <- GetKern(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma, z=z[i])
    pd.p[i,1:J] <- GetPd(kern=kern[i,1:J],p0=p0.p,J=J,z=z[i])
    pd.c[i,1:J] <- GetPd(kern=kern[i,1:J],p0=p0.c,J=J,z=z[i])
    lambda.p[i,1:J] <- -log(1-pd.p[i,1:J]) #capture time rate parameters, 1st cap
    lambda.c[i,1:J] <- -log(1-pd.c[i,1:J]) #capture time rate parameters, subsequent cap
    #partially latent true capture history.
    #all events that would have happened if traps did not fill up
    #all detection events that have realized capture time < 1
    y.true[i,1:J,1:K] ~ dBernoulliMatrixMb(pd.p=pd.p[i,1:J],pd.c=pd.c[i,1:J],y.state[i,1:J,1:K],K2D=K2D[1:J,1:K],z=z[i])
  }
  #model for captures we observe given partially latent captures and capture order
  #capture order is also latent and updated
  for(k in 1:K){
    y.obs[1:n.cap,1:J,k] ~ dThin(y.true=y.true[1:M,1:J,k],y.state=y.state[1:M,1:J,k],
                             lambda.p=lambda.p[1:M,1:J],lambda.c=lambda.c[1:M,1:J],
                             order=order2D[1:n.obs.cells[k],k],obs.i=obs.i2D[1:n.obs.cells[k],k],
                             obs.j=obs.j2D[1:n.obs.cells[k],k],n.cap=n.cap)
  }
  y.true.sum <- GetSum(y.true[1:M,1:J,1:K],z=z[1:M])
})