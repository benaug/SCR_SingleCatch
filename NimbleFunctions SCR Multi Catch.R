GetPd <- nimbleFunction(
  run = function(s = double(1), J = double(0), p0 = double(0), sigma = double(0), 
                 X = double(2), z = double(0)){ 
    returnType(double(1))
    if(z==0){
      pd <- rep(0,J)
    }else{
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      pd <- p0*exp(-d2/(2*sigma^2))
    }
    return(pd)
  }
)

#This version faster than one below
dObsMatrix <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K2D = double(2), K1D = double(1),
                 K = double(0), z = double(0), log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      J <- nimDim(K2D)[1]
      lambda <- rep(0,J)
      logProb <- 0
      #likelihood assuming no captures on any occasion
      for(j in 1:J){
        lambda[j] <- -log1p(-pd[j]) #numerically stable for small pd
        logProb <- logProb - lambda[j]*K1D[j]
      }
      #replace noncapture likelihood with capture likelihood where captured
      for(k in 1:K){
        if(x[k]>0){
          lambda.dot <- 0
          for(j in 1:J){
            if(K2D[j,k]==1){
              lambda.dot <- lambda.dot + lambda[j]
            }
          }
          if(K2D[x[k],k]==1){
            lambda.cap <- lambda[x[k]]
            logProb <- logProb + lambda.dot + log(lambda.cap) - log(lambda.dot) +
              log(1-exp(-lambda.dot))
          }else{
            return(-Inf)
          }
        }
      }
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rObsMatrix <- nimbleFunction(
  run = function(n = integer(0), pd = double(1), K2D = double(2),  K1D = double(1), K = double(0), z = double(0)) {
    returnType(double(1))
    K <- nimDim(K2D)[2]
    out <- rep(0,K)
    return(out)
  }
)

# dObsMatrix <- nimbleFunction(
#   run = function(x = double(1), pd = double(1), K2D = double(2), K = double(0), z = double(0),
#                  log = integer(0)) {
#     returnType(double(0))
#     if(z==0){
#       return(0)
#     }else{
#       J <- nimDim(K2D)[1]
#       lambda <- rep(0,J)
#       for(j in 1:J){
#         lambda[j] <- -log1p(-pd[j]) #numerically stable for small pd
#       }
#       logProb <- 0
#       for(k in 1:K){
#         lambda.dot <- 0
#         for(j in 1:J){
#           lambda.dot <- lambda.dot + lambda[j]*K2D[j,k]
#         }
#         if(x[k]>0){
#           lambda.cap <- lambda[x[k]]*K2D[x[k],k]
#           logProb <- logProb + log(lambda.cap) - log(lambda.dot) +
#             log(1-exp(-lambda.dot))
#         }else{
#           logProb <- logProb - lambda.dot
#         }
#       }
#       return(logProb)
#     }
#   }
# )
# 
# #make dummy random vector generator to make nimble happy
# rObsMatrix <- nimbleFunction(
#   run = function(n = integer(0), pd = double(1), K2D = double(2), K = double(0), z = double(0)) {
#     returnType(double(1))
#     K <- nimDim(K2D)[2]
#     out <- rep(0,K)
#     return(out)
#   }
# )

zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    z.ups <- control$z.ups
    M <- control$M
    K <- control$K
    inds.detected <- control$inds.detected
    #nodes used for update
    y.nodes <- model$expandNodeNames("y")
    N.node <- model$expandNodeNames("N")
    z.nodes <- model$expandNodeNames("z")
    pd.nodes <- model$expandNodeNames(paste("pd"))
    calcNodes <- c(N.node,z.nodes,pd.nodes,y.nodes)
  },
  run = function(){
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a captured individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn pd off
          model$calculate(pd.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(pd.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # copy(from = model, to = mvSaved, row = 1, nodes = z.nodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)