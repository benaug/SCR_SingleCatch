dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

GetKern <- nimbleFunction(
  run = function(s = double(1), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      kern <- exp(-d2/(2*sigma^2))
      return(kern)
    }
  }
)

dObsMatrix <- nimbleFunction(
  run = function(x = double(1), y.state = double(2), kern = double(1),
                 p0.p = double(0), p0.c = double(0),
                 K2D = double(2), K1D.p = double(1), K1D.c = double(1),
                 K = double(0), z = double(0), log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      J <- nimDim(K2D)[1]
      lambda.p <- rep(0,J)
      lambda.c <- rep(0,J)
      logProb <- 0
      
      #likelihood assuming no captures on any occasion
      for(j in 1:J){
        if(K1D.p[j]>0){
          lambda.p[j] <- -log1p(-p0.p*kern[j])
          logProb <- logProb - lambda.p[j]*K1D.p[j]
        }
        if(K1D.c[j]>0){
          lambda.c[j] <- -log1p(-p0.c*kern[j])
          logProb <- logProb - lambda.c[j]*K1D.c[j]
        }
      }
      #replace noncapture likelihood with capture likelihood where captured
      for(k in 1:K){
        if(x[k]>0){
          lambda.dot <- 0
          for(j in 1:J){
            if(K2D[j,k]==1){
              if(y.state[j,k]==0){
                lambda.dot <- lambda.dot + lambda.p[j]
              }else{
                lambda.dot <- lambda.dot + lambda.c[j]
              }
            }
          }
          
          if(K2D[x[k],k]==1){
            if(y.state[x[k],k]==0){
              lambda.cap <- lambda.p[x[k]]
            }else{
              lambda.cap <- lambda.c[x[k]]
            }
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
  run = function(n = integer(0), y.state = double(2), kern = double(1),
                 p0.p = double(0), p0.c = double(0),
                 K2D = double(2), K1D.p = double(1), K1D.c = double(1),
                 K = double(0), z = double(0)) {
    returnType(double(1))
    out <- rep(0,K)
    return(out)
  }
)

zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    z.ups <- control$z.ups
    M <- control$M
    inds.detected <- control$inds.detected
    #convert detected indices to an indicator once in setup
    ind.detected <- rep(0,M)
    ind.detected[inds.detected] <- 1
    #nodes used for update
    y.nodes <- model$expandNodeNames("y")
    N.node <- model$expandNodeNames("N")
    z.nodes <- model$expandNodeNames("z")
    kern.nodes <- model$expandNodeNames(paste("kern"))
    calcNodes <- c(N.node,z.nodes,kern.nodes,y.nodes)
  },
  run = function(){
    #build undetected on/off lists once, then update them after accepted proposals
    #detected individuals are never proposed off, avoiding automatic rejections
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(ind.detected[i]==0){
        if(model$z[i]==1){
          non.curr <- non.curr+1
          z.on[non.curr] <- i
        }else{
          noff.curr <- noff.curr+1
          z.off[noff.curr] <- i
        }
      }
    }
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        non.init <- non.curr
        if(non.init>0){ 
          pick.pos <- rcat(1,rep(1/non.init,non.init)) 
          pick <- z.on[pick.pos]
          N.init <- model$N[1] #needed for proposal/combinatorial correction
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<- model$N[1] - 1
          model$z[pick] <<- 0
          
          model$calculate(kern.nodes[pick]) #turn kern off
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- 0 
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) -
            (lp.initial.N + lp.initial.y) + log(non.init/N.init)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            #calculate y now to synchronize accepted logProb
            model$calculate(y.nodes[pick])
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["kern",1][pick,] <<- model[["kern"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            #move accepted individual from on list to off list
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["kern"]][pick,] <<- mvSaved["kern",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(N.node)
          }
        }
        
      }else{#add
        noff.init <- noff.curr
        if(noff.init>0){
          pick.pos <- rcat(1,rep(1/noff.init,noff.init))
          pick <- z.off[pick.pos]
          N.init <- model$N[1] #needed for proposal/combinatorial correction
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- 0
          
          #propose new N/z
          model$N[1] <<- model$N[1] + 1
          model$z[pick] <<- 1
          
          model$calculate(kern.nodes[pick]) #changed: turn kern on
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) -
            (lp.initial.N + lp.initial.y) + log((N.init+1)/(non.curr+1))
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["kern",1][pick,] <<- model[["kern"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            #move accepted individual from off list to on list
            z.off[pick.pos] <- z.off[noff.curr]
            z.off[noff.curr] <- 0
            noff.curr <- noff.curr-1
            non.curr <- non.curr+1
            z.on[non.curr] <- pick
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["kern"]][pick,] <<- mvSaved["kern",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    
    #copy back to mvSaved to update logProbs
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)