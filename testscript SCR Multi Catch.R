#competing hazards model for multi catch traps of Efford and Borchers 2009
# https://link.springer.com/chapter/10.1007/978-0-387-78151-8_11
# Efford, Murray G., David L. Borchers, and Andrea E. Byrom. 
# "Density estimation by spatially explicit capture–recapture: likelihood-based methods."
# Modeling demographic processes in marked populations (2009): 255-269.

library(nimble)
library(coda)
source("sim.SCR.multicatch.R") # data simulator
source("NimbleModel SCR Multi Catch.R") #nimble model file
source("NimbleFunctions SCR Multi Catch.R") #nimble functions and custom updates
source("sSampler.R") #custom activity center update

#Simulate some data
N <- 50 #Abundance
p0 <- 0.25 #baseline detection
sigma <- 0.50 #spatial scale
K <- 5 #number of occasions
X <- as.matrix(expand.grid(2:8,2:8)) #traps 7 x 7 here
buff <- 2 #state space buffer around maximal trap extent

data <- sim.SCR.multiCatch(N=N,p0=p0,sigma=sigma,X=X,buff=buff) #simulates complete trap operation

#What is the observed data?
str(data$y) #n x J x K detection history
#after augmenting n up to M, we will convert to M x K detection history, recording trap of capture on occasion k, 0 if not captured

#Fit model
M <- 125 #set data augmentation limit
J <- nrow(data$X) #number of traps

#initialize y. M x J x K detection history here, used to initialize s. then convert to M x K.
y <- array(0,dim=c(M,J,K))
y[1:data$n.cap,,] <- data$y

#initialize z and N. N.init must be consistent with z.init
z.init <- 1*(rowSums(y)>0)
N.init <- sum(z.init)
#initialize s
s.init <- cbind(runif(M,data$xlim[1],data$xlim[2]),runif(M,data$xlim[1],data$xlim[2]))
y2D <- apply(y,c(1,2),sum)
idx <- which(rowSums(y2D)>0) #switch for those caught or with latent captures on initialization
for(i in idx){
  trps <- matrix(X[y2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,] <- trps
  }
}

#convert to trap of capture. 0 if not captured
y.use <- matrix(0,M,K)
for(i in 1:M){
  for(k in 1:K){
    trapcap <- which(y[i,,k]==1)
    if(length(trapcap)>0){
      y.use[i,k] <- trapcap
    }
  }
}

y <- y.use


#initial values for nimble
Niminits <- list(lambda.N=N.init,p0=runif(1,0.1,0.9),sigma=runif(1,0.5,1),
                 s=s.init,z=z.init,N=N.init)

#constants for nimble
constants <- list(M=M,J=J,K=K,K2D=data$K2D,xlim=data$xlim,ylim=data$ylim)

#supply data to nimble
Nimdata <- list(X=data$X,y=y)

# set parameters to monitor
parameters <- c('p0','sigma','lambda.N','N')
nt <- 1 #thinning rate
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('p0','sigma','lambda.N')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,nodes=config.nodes)

#add sampler for N/z
z.ups <- round(M*0.25) # how many z proposals per iteration? 25% of M generally seems good, but no idea what is optimal
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,K=K,inds.detected=which(rowSums(y)>0)),
                silent = TRUE)

#adding this block sampler might help, check posterior correlation
conf$addSampler(target = c("p0","sigma"),
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

Cmodel$calculate() #make sure starting logProb is finite

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

cor(mvSamples[250:nrow(mvSamples),]) #posterior correlation

