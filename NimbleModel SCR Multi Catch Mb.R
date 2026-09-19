NimModel <- nimbleCode({
  p0.p ~ dunif(0,1) #first capture p0
  p0.c ~ dunif(0,1) #subsequent capture p0
  sigma ~ dunif(0,100)
  lambda.N ~ dunif(0,1000)
  N ~ dpois(lambda.N) #realized N in state space
  for(i in 1:M){ #N/z update under the hood
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    kern[i,1:J] <- GetKern(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma, z=z[i])
    #detection data are trap of capture on each occasion, 0 if not captured
    y[i,1:K] ~ dObsMatrix(y.state=y.state[i,1:J,1:K],kern=kern[i,1:J],
                          p0.p=p0.p,p0.c=p0.c,K2D=K2D[1:J,1:K],
                          K1D.p=K1D.p[i,1:J],K1D.c=K1D.c[i,1:J],K=K,z=z[i])
  }
})
