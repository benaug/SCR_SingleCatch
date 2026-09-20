NimModel <- nimbleCode({
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,100)
  lambda.N ~ dunif(0,1000)
  N ~ dpois(lambda.N) #realized N in state space
  for(i in 1:M){#N/z and y.true update under the hood
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    kern[i,1:J] <- GetKern(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma,z=z[i])
    lambda[i,1:J] <- -log1p(-p0*kern[i,1:J]) #capture time rate parameters
    #partially latent true capture history.
    #all events that would have happened if traps did not fill up
    #all detection events that have realized capture time < 1
    y.true[i,1:J,1:K] ~ dBernoulliMatrix(kern=kern[i,1:J],p0=p0,
                                         K2D=K2D[1:J,1:K],K1D=K1D[1:J],z=z[i])
  }
  #model for captures we observe given partially latent captures and capture order
  #capture order is also latent and updated
  for(k in 1:K){
    y.obs[1:n.cap,1:J,k] ~ dThin(y.true=y.true[1:M,1:J,k],order=order2D[1:n.obs.cells.max,k],
                                 lambda=lambda[1:M,1:J],obs.i=obs.i2D[1:n.obs.cells.max,k],
                                 obs.j=obs.j2D[1:n.obs.cells.max,k],n.obs=n.obs.cells[k],n.cap=n.cap)
  }
  y.true.sum <- GetSum(y.true[1:M,1:J,1:K],z=z[1:M])
})
