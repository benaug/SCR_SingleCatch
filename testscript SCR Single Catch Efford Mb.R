#This version has a behavioral response to capture
#The assumption is that an individual has to be actually
#captured, not "latently captured" to undergo a trap response.
#Therefore, the capture states (first/subsequent) are observed from
#the observed capture history.
#Observation model parameters are estimated quite imprecisely.
#I haven't tested this version, yet.

library(nimble)
library(coda)
source("sim.SCR.singleCatch.Mb.R") # data simulator
source("NimbleModel SCR Single Catch Efford Mb.R") #nimble model file
source("NimbleFunctions SCR Single Catch Efford Mb.R") #nimble functions and custom updates
source("sSampler.R") #custom activity center update

#Simulate some data
N <- 50 #Abundance
p0.p <- 0.25 #baseline detection, first capture
p0.c <- 0.50 #baseline detection, subsequent capture
sigma <- 0.50 #spatial scale
K <- 5 #number of occasions
X <- as.matrix(expand.grid(2:8,2:8)) #traps 7 x 7 here
buff <- 2 #state space buffer around maximal trap extent

data <- sim.SCR.singleCatch.Mb(N=N,p0.p=p0.p,p0.c=p0.c,sigma=sigma,X=X,buff=buff)

#What is the observed data?
str(data$y.obs) #observed capture history
str(data$y.state) #capture states (first vs subsequent)

#Fit model
M <- 100 #set data augmentation limit
J <- nrow(data$X) #number of traps

#initialize y.true as y.obs
y.true.init <- array(0,dim=c(M,J,K))
y.true.init[1:data$n.cap,,] <- data$y.obs
#augment y.state. y.state is observed data
y.state <- array(0,dim=c(M,J,K))
y.state[1:data$n.cap,,] <- data$y.state
#initialize z and N. N.init must be consistent with z.init
z.init <- 1*(rowSums(y.true.init)>0)
N.init <- sum(z.init)
#initialize s
s.init <- cbind(runif(M,data$xlim[1],data$xlim[2]),runif(M,data$xlim[1],data$xlim[2]))
y.true.init2D <- apply(y.true.init,c(1,2),sum)
idx <- which(rowSums(y.true.init2D)>0) #switch for those caught or with latent captures on initialization
for(i in idx){
  trps <- matrix(X[y.true.init2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,] <- trps
  }
}
#initialize capture order
order2D.init <- matrix(NA,max(data$n.obs.cells),K)
for(k in 1:K){
  order2D.init[1:data$n.obs.cells[k],k] <- sample(1:data$n.obs.cells[k],data$n.obs.cells[k],replace=FALSE)
}

#initial values for nimble
Niminits <- list(lambda.N=N.init,p0.p=runif(1,0.1,0.5),p0.c=runif(1,0.3,0.9),sigma=runif(1,0.5,1),
                 s=s.init,y.true=y.true.init,z=z.init,N=N.init,
                 order2D=order2D.init)

#constants for nimble
constants <- list(M=M,J=J,K=K,K2D=data$K2D,xlim=data$xlim,ylim=data$ylim,n.cap=data$n.cap,
                  obs.i2D=data$obs.i2D,obs.j2D=data$obs.j2D,n.obs.cells=data$n.obs.cells)

#supply data to nimble
Nimdata <- list(X=data$X,y.obs=data$y.obs,y.state=y.state)

# set parameters to monitor
parameters <- c('p0.p','p0.c','sigma','lambda.N','N','y.true.sum')
nt <- 1 #thinning rate
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('p0.p','p0.c','sigma','lambda.N')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,nodes=config.nodes)

#add sampler for y.true
y.ups <- 2 #no idea what is optimal. 1 might be fine choice.
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'ySampler',control = list(M=M,J=J,K=K,K2D=data$K2D,y.obs=data$y.obs,n.cap=data$n.cap,
                                                 obs.i=data$obs.i,obs.j=data$obs.j,obs.k=data$obs.k,
                                                 obs.i2D=data$obs.i2D,obs.j2D=data$obs.j2D,n.obs.cells=data$n.obs.cells,
                                                 y.ups=y.ups),
                silent = TRUE)

#add sampler for N/z
z.ups <- round(M*0.25) # how many z proposals per iteration? 25% of M generally seems good, but no idea what is optimal
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M),
                silent = TRUE)

#probably worth adding this block sampler
conf$addSampler(target = c("p0.p","sigma"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

#add activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=data$xlim,ylim=data$ylim,
                                                 scale=1,adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#if p0 and (especially) sigma too small, starting logProb will return NaN. Don't run, nimble will crash.
#raise sigma and maybe p0 init.
Cmodel$calculate()

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

sum(data$y.true) #true number of latent captures

cor(mvSamples[250:nrow(mvSamples),]) #posterior correlation

#can plot captures from final iteration
par(mfrow=c(1,1),ask=FALSE)
plot(X,pch=4,xlim=data$xlim,ylim=data$ylim,xlab="X",ylab="Y",
     main="Observed (black) and latent (gold) captures")
y.true2D <- apply(Cmodel$y.true,c(1,2),sum)
for(i in 1:M){
  if(Cmodel$z[i]==1){
    trapcaps <- which(y.true2D[i,]>0)
    s.tmp <- Cmodel$s[i,]
    points(s.tmp[1],s.tmp[2],pch=16)
    for(j in trapcaps){
      lines(x=c(s.tmp[1],X[j,1]),y=c(s.tmp[2],X[j,2]),lwd=2,col="goldenrod") #latent events
    }
  }
}
y.obs2D <- apply(data$y.obs,c(1,2),sum)
for(i in 1:data$n.cap){
  trapcaps <- which(y.obs2D[i,]>0)
  s.tmp <- Cmodel$s[i,]
  for(j in trapcaps){
    lines(x=c(s.tmp[1],X[j,1]),y=c(s.tmp[2],X[j,2]),lwd=2,col="black") #observed events
  }
}
