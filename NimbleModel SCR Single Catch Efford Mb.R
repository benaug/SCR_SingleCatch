NimModel <- nimbleCode({
  p0.p ~ dunif(0,1) #first capture p0
  p0.c ~ dunif(0,1) #subsequent capture p0
  sigma ~ dunif(0,100)
  lambda.N ~ dunif(0,1000) #expected N
  N ~ dpois(lambda.N) #realized N
  for(i in 1:M){ #N/z and y.true update under the hood
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    kern[i,1:J] <- GetKern(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma, z=z[i])
    lambda.p[i,1:J] <- -log1p(-p0.p*kern[i,1:J]) #capture time rate parameters, 1st cap
    lambda.c[i,1:J] <- -log1p(-p0.c*kern[i,1:J]) #capture time rate parameters, subsequent cap
    #partially latent true capture history.
    #all events that would have happened if traps did not fill up
    #all detection events that have realized capture time < 1
    y.true[i,1:J,1:K] ~ dBernoulliMatrixMb(kern=kern[i,1:J],p0.p=p0.p,p0.c=p0.c,
                                           y.state=y.state[i,1:J,1:K],K2D=K2D[1:J,1:K],
                                           K1D.p=K1D.p[i,1:J],K1D.c=K1D.c[i,1:J],z=z[i])
  }
  #model for captures we observe given partially latent captures and capture order
  #capture order is also latent and updated
  for(k in 1:K){
    y.obs[1:n.cap,1:J,k] ~ dThin(y.true=y.true[1:M,1:J,k],y.state=y.state[1:M,1:J,k],
                                 lambda.p=lambda.p[1:M,1:J],lambda.c=lambda.c[1:M,1:J],
                                 order=order2D[1:n.obs.cells.max,k],obs.i=obs.i2D[1:n.obs.cells.max,k],
                                 obs.j=obs.j2D[1:n.obs.cells.max,k],n.obs=n.obs.cells[k],n.cap=n.cap)
  }
  y.true.sum <- GetSum(y.true[1:M,1:J,1:K],z=z[1:M])
})