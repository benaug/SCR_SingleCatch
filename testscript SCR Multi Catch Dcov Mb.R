#competing hazards model for multi catch traps of Efford and Borchers 2009
# https://link.springer.com/chapter/10.1007/978-0-387-78151-8_11
# Efford, Murray G., David L. Borchers, and Andrea E. Byrom. 
# "Density estimation by spatially explicit capture–recapture: likelihood-based methods."
# Modeling demographic processes in marked populations (2009): 255-269.

library(nimble)
nimbleOptions(determinePredictiveNodesInModel = FALSE) #must run this line
library(coda)
source("sim.SCR.multicatch.Dcov.Mb.R") # data simulator
source("NimbleModel SCR Multi Catch Dcov Mb.R") #nimble model file
source("NimbleFunctions SCR Multi Catch Dcov Mb.R") #nimble functions and custom updates
source("sSampler Dcov.R") #custom activity center update

#Simulate some data
N <- 50 #Abundance
p0.p <- 0.25 #baseline detection, first detection
p0.c <- 0.50 #baseline detection, subsequent detection
sigma <- 0.50 #spatial scale
K <- 5 #number of occasions
X <- as.matrix(expand.grid(2:8,2:8)) #traps 7 x 7 here
buff <- 2 #state space buffer around maximal trap extent

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
# set.seed(13210) #pretty good one
set.seed(13216)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(100,100),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#Additionally, maybe we want to exclude "non-habitat" or limit the state space extent
#let's use a 3sigma buffer
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(0,length(D.cov))
dists <- e2dist(X,dSS.tmp)
min.dists <- apply(dists,2,min)
InSS[min.dists<(3*sigma)] <- 1
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="Habitat",col=cols1)
points(X,pch=4)

#Density covariates
D.beta0 <- -0.7
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell*InSS,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X,pch=4)

set.seed(399403) #setting new seed here since we set the same seed for D.cov above. Change to get new data set
data <- sim.SCR.multiCatch.Dcov.Mb(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                              xlim=xlim,ylim=ylim,res=res,p0.p=p0.p,p0.c=p0.c,sigma=sigma,X=X) #simulates complete trap operation

#What is the observed data?
str(data$y) #n x J x K detection history
str(data$y.state) #n x J x K capture state history
#after augmenting n up to M, we will convert to M x K detection history, recording trap of capture on occasion k, 0 if not captured

#Fit model
M <- 150 #set data augmentation limit
J <- nrow(data$X) #number of traps

#initialize y. M x J x K detection history here, used to initialize s. then convert to M x K.
y <- array(0,dim=c(M,J,K))
y[1:data$n.cap,,] <- data$y

#augment y.state
y.state <- array(0,dim=c(M,J,K))
y.state[1:data$n.cap,,] <- data$y.state

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

#If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
alldists <- e2dist(s.init,data$dSS)
alldists[,data$InSS==0] <- Inf
for(i in 1:M){
  this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
  if(data$InSS[this.cell]==0){
    cands <- alldists[i,]
    new.cell <- which(alldists[i,]==min(alldists[i,]))
    s.init[i,] <- data$dSS[new.cell,]
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
Niminits <- list(N=N.init,lambda.N=N.init,D0=sum(z.init)/(sum(data$InSS)*data$res^2),D.beta1=0,
                 p0.p=runif(1,0.1,0.9),p0.c=runif(1,0.1,0.9),
                 sigma=runif(1,0.5,1),s=s.init,z=z.init)

#constants for nimble
constants <- list(M=M,J=J,K=K,K2D=data$K2D,D.cov=data$D.cov,cellArea=data$cellArea,
                  n.cells=data$n.cells,xlim=data$xlim,ylim=data$ylim,res=data$res)

#supply data to nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(X=data$X,y=y,y.state=y.state,dummy.data=dummy.data,cells=cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('p0.p','p0.c','sigma','D0','D.beta1','lambda.N','N')
parameters2 <- c("lambda.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 5 #thinning rate for parameters2

start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('p0.p','p0.c','sigma')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2, thin2=nt2,,
                      useConjugacy = FALSE,nodes=config.nodes)

#add sampler for N/z
z.ups <- round(M*0.25) # how many z proposals per iteration? 25% of M generally seems good, but no idea what is optimal
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,K=K,inds.detected=which(rowSums(y)>0)),
                silent = TRUE)

#adding this block sampler might help, check posterior correlation
conf$addSampler(target = c("p0.p","sigma"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

#add activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerDcov',control=list(i=i,res=res,xlim=data$xlim,ylim=data$ylim,
                                                     n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                     scale=1,adaptive=TRUE),silent = TRUE)
}

#add blocked sampler for density parameters. AF_slice mixes better than RW_block, often more ESS/time
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)

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

data$N #true realized abundance
data$lambda.N #true expected abundance

cor(mvSamples[250:nrow(mvSamples),]) #posterior correlation

#Plot density surface
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10
#image will show posterior means
lambda.cell.post <- cellArea*mvSamples2[burnin2:nrow(mvSamples2),D0.idx]*mvSamples2[burnin2:nrow(mvSamples2),lambda.cell.idx]
lambda.cell.ests <- colMeans(lambda.cell.post)
#remove non-habitat
lambda.cell.ests[InSS==0] <- NA
lambda.cell[InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
