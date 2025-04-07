#this version has a spatial density covariate

library(nimble)
nimbleOptions(determinePredictiveNodesInModel = FALSE) #must run this line
library(coda)
source("sim.SCR.singleCatch.Dcov.R") # data simulator
source("NimbleModel SCR Single Catch Efford Dcov.R") #nimble model file
source("NimbleFunctions SCR Single Catch Efford Dcov.R") #nimble functions and custom updates
source("sSampler Dcov.R") #custom activity center update

#Simulate some data
p0 <- 0.5 #baseline detection
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

set.seed(399405) #setting new seed here since we set the same seed for D.cov above. Change to get new data set
data <- sim.SCR.singleCatch.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                                 xlim=xlim,ylim=ylim,res=res,p0=p0,sigma=sigma,X=X)

#What is the observed data?
str(data$y.obs) #observed capture history

data$N #realized N, make sure M is high enough

#Fit model
M <- 125 #set data augmentation limit
J <- nrow(data$X) #number of traps

#initialize y.true as y.obs
y.true.init <- array(0,dim=c(M,J,K))
y.true.init[1:data$n.cap,,] <- data$y.obs
#initialize z and N. N.init must be consistent with z.init
z.init <- 1*(rowSums(y.true.init)>0)
N.init <- sum(z.init)
#initialize s
xlim <- data$xlim
ylim <- data$ylim
s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
y.true.init2D <- apply(y.true.init,c(1,2),sum)
idx <- which(rowSums(y.true.init2D)>0) #switch for those actually caught
for(i in idx){
  trps <- matrix(X[y.true.init2D[i,]>0,1:2],ncol=2,byrow=FALSE)
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

#initialize capture order
order2D.init <- matrix(NA,max(data$n.obs.cells),K)
for(k in 1:K){
  order2D.init[1:data$n.obs.cells[k],k] <- sample(1:data$n.obs.cells[k],data$n.obs.cells[k],replace=FALSE)
}

#initial values for nimble
Niminits <- list(N=N.init,lambda.N=N.init,D0=sum(z.init)/(sum(data$InSS)*data$res^2),D.beta1=0,
                 s=s.init,y.true=y.true.init,z=z.init,order2D=order2D.init,
                 p0=runif(1,0.1,0.9),sigma=runif(1,0.5,1))

#constants for nimble
constants <- list(M=M,J=J,K=K,K2D=data$K2D,n.cap=data$n.cap,
                  obs.i2D=data$obs.i2D,obs.j2D=data$obs.j2D,n.obs.cells=data$n.obs.cells,
                  D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  xlim=data$xlim,ylim=data$ylim,res=data$res)

#supply data to nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y.obs=data$y.obs,X=data$X,dummy.data=dummy.data,cells=cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('D0','D.beta1','p0','sigma','lambda.N','N','y.true.sum')
parameters2 <- c("lambda.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 5 #thinning rate for paremeters2


start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('p0','sigma')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2, thin2=nt2,
                      useConjugacy = FALSE,nodes=config.nodes)

#add blocked sampler for denstity parameters. AF_slice mixes better than RW_block, often more ESS/time
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)

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
conf$addSampler(target = c("p0","sigma"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

#add activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerDcov',control=list(i=i,res=res,xlim=data$xlim,ylim=data$ylim,
                                                     n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                 scale=1,adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#if p0 and (expecially) sigma too small, starting logProb will return NaN. Don't run, nimble will crash.
#raise sigma and maybe p0 init.
Cmodel$calculate() 

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(1000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

data$N #true realized abundance
data$lambda.N #true expected abundance
sum(data$y.true) #true number of latent captures


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
