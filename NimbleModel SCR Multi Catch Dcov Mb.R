NimModel <- nimbleCode({
  p0.p ~ dunif(0,1) #first capture p0
  p0.c ~ dunif(0,1) #subsequent capture p0
  sigma ~ dunif(0,100)
  #Density parameters
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #Density model
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space
  for(i in 1:M){ #N/z update under the hood
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~  dunif(xlim[1],xlim[2])
    s[i,2] ~  dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InSS=InSS[s.cell[i]])
    kern[i,1:J] <- GetKern(s=s[i,1:2],X=X[1:J,1:2],J=J,sigma=sigma, z=z[i])
    pd.p[i,1:J] <- GetPd(kern=kern[i,1:J],p0=p0.p,J=J,z=z[i])
    pd.c[i,1:J] <- GetPd(kern=kern[i,1:J],p0=p0.c,J=J,z=z[i])
    pd.multi[i,1:J,1:K] <- GetPdMulti(y.state[i,1:J,1:K],pd.p=pd.p[i,1:J],
                                      pd.c=pd.c[i,1:J], K2D=K2D[1:J,1:K],z=z[i])
    #detection data are trap of capture on each occasion, 0 if not captured
    y[i,1:K] ~ dObsMatrix(pd.multi=pd.multi[i,1:J,1:K],K2D=K2D[1:J,1:K],K=K,z=z[i])
  }
})
