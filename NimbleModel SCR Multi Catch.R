NimModel <- nimbleCode({
  p0 ~ dunif(0,1)
  sigma ~ dunif(0,100)
  lambda.N ~ dunif(0,1000)
  N ~ dpois(lambda.N) #realized N in state space
  for(i in 1:M){ #N/z update under the hood
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    pd[i,1:J] <- GetPd(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma,p0=p0,z=z[i])
    pd.multi[i,1:J,1:K] <- GetPdMulti(pd=pd[i,1:J], K2D=K2D[1:J,1:K],z=z[i])
    #detection data are trap of capture on each occasion, 0 if not captured
    y[i,1:K] ~ dObsMatrix(pd.multi=pd.multi[i,1:J,1:K],K2D=K2D[1:J,1:K],K=K,z=z[i])
  }
})
