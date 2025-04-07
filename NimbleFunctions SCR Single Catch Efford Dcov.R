dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0), InSS = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(InSS==1){
      logProb <- log(pi.cell)
    }else{
      logProb <- -Inf
    }
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0), InSS = double(0)) {
    returnType(double(0))
    return(0)
  }
)
#sum up total number of latent capture events
GetSum <- nimbleFunction(
  run = function(y.true = double(3),z = double(1)){ 
    returnType(double(0))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    K <- nimDim(y.true)[3]
    y.true.sum <- 0
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            y.true.sum <- y.true.sum + y.true[i,j,k]
          }
        }
      }
    }
    return(y.true.sum)
  }
)

GetPd <- nimbleFunction(
  run = function(s = double(1), p0 = double(0), sigma = double(0), 
                 X = double(2), J = double(0),z = double(0)){ 
    returnType(double(1))
    if(z==0){
      return(rep(0,J))
    }else{
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      pd <- p0*exp(-d2/(2*sigma^2))
    }
    return(pd)
  }
)

dBernoulliMatrix <- nimbleFunction(
  run = function(x = double(2), pd = double(1), K2D = double(2), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){#skip calculation if z=0
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, size = K2D, p = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBernoulliMatrix <- nimbleFunction(
  run = function(n = integer(0), pd = double(1), K2D = double(2), z = double(0)) {
    returnType(double(2))
    J <- nimDim(K2D)[1]
    K <- nimDim(K2D)[2]
    out <- matrix(0,J,K)
    return(out)
  }
)

#used in pSmaller() below
integrand <- nimbleFunction(
  run = function(x = double(1),param = double(1)){ 
    returnType(double(1))
    n1 <- length(param)
    lambda1 <- param[1]
    lambda2 <- param[2:n1]
    n <- length(lambda2)
    n.x <- length(x)
    prod.term <- rep(1,n.x)
    for(j in 1:n.x){
      for(i in 1:n){
        prod.term[j] <- prod.term[j] * (exp(-lambda2[i] * x[j]) - exp(-lambda2[i]))
      }
    }
    for(j in 1:n.x){
      prod.term[j] <- prod.term[j] * exp(-lambda1 * x[j])
    }
    return(prod.term)
  })

#probability exponential RV right-truncated at 1 with parameter lambda1 is less than one or more other exponential
#RVs right-truncated at 1 with parameter(s) lambda2
pSmaller <- nimbleFunction(
  run = function(lambda1 = double(0), lambda2 = double(1), log = integer(0)) {
    returnType(double(0))
    param = c(lambda1,lambda2)
    #integral from 0 to 1
    integral <- nimIntegrate(integrand, lower = 0, upper = 1, param = param)[1]
    #denominator terms
    lambda.term <- 1 - exp(-lambda1)
    lambda2.prod <- prod(1 - exp(-lambda2))
    # prob <- (lambda1 / (lambda.term * lambda2.prod)) * integral
    logProb <- log(lambda1*integral) - log(lambda.term * lambda2.prod) #less likely to underflow
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  })

dThin <- nimbleFunction(
  run = function(x = double(2), y.true = double(2), lambda = double(2), obs.i = double(1),
                 obs.j = double(1), order = double(1), n.cap = double(0), log = integer(0)) {
    returnType(double(0))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    lambda.tmp <- lambda #so we don't write over lambda in c++
    logProb <- 0
    # for(o in 1:length(order)){ 
    for(o in 1:(length(order)-1)){ #we do not need the final logProb which is always 0
      idx <- which(order==o)[1]
      focal.lambda <- lambda.tmp[obs.i[idx],obs.j[idx]]
      n.other.lambdas <- sum(y.true==1&lambda.tmp<Inf)-1
      #excluding lambdas of 0. leads to nonfinite logProb, these inds will never be captured so they cannot get there first
      if(n.other.lambdas>0){
        other.lambdas <- rep(0,n.other.lambdas) #lambda < Inf is not using traps removed below on next loop iteration
        idx2 <- 1
        for(i in 1:M){
          for(j in 1:J){
            if(y.true[i,j]==1&lambda.tmp[i,j]<Inf){ #if a latent capture
              if(!(i==obs.i[idx]&j==obs.j[idx])){ #don't include focal
                other.lambdas[idx2] <- lambda.tmp[i,j]
                idx2 <- idx2 + 1
              }
            }
          }
        }
        logProb <- logProb + pSmaller(focal.lambda,other.lambdas,log=TRUE)
      } #else add logProb of 0. But we are just skipping the last index in the o loop
      #zero out this individual and trap
      lambda.tmp[obs.i[idx],] <- Inf
      lambda.tmp[,obs.j[idx]] <- Inf
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
    return(logProb)
  })

rThin <- nimbleFunction(
  run = function(n = integer(0), y.true = double(2), lambda = double(2), obs.i = double(1),
                 obs.j = double(1), order = double(1), n.cap = double(0)){
    returnType(double(2))
    J <- nimDim(y.true)[2]
    return(matrix(0,n.cap,J))
  }
)

ySampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    y.ups <- control$y.ups
    M <- control$M
    J <- control$J
    K <- control$K
    obs.i <- control$obs.i
    obs.j <- control$obs.j
    obs.k <- control$obs.k
    obs.i2D <- control$obs.i2D
    obs.j2D <- control$obs.j2D
    n.obs.cells <- control$n.obs.cells
    K2D <- control$K2D
    y.obs <- control$y.obs
    n.cap <- control$n.cap
    calcNodes <- model$getDependencies(c("y.true","y.obs","order2D"))
  },
  run = function(){
    y.true <- model$y.true
    z <- model$z
    pd <- model$pd
    lambda <- model$lambda
    order2D <- model$order2D
    ll.y <- array(0,dim=c(M,J,K)) #computing this instead of pulling out of model because it is 1D in model
    for(k in 1:K){
      for(j in 1:J){
        ll.y[,j,k] <- dbinom(y.true[,j,k],K2D[j,k],pd[,j],log=TRUE)
      }
    }
    ll.y.cand <- ll.y
    ll.y.obs <- model$logProb_y.obs[1,1,]
    ll.y.obs.cand <- ll.y.obs
    y.true.cand <- y.true
    n.obs.cells.all <- sum(n.obs.cells)
    for(up in 1:y.ups){ #update one or more times per iteration
      #update y.true for cells with y.obs=1
      for(c in 1:n.obs.cells.all){
        skip <- FALSE
        updown <- rbinom(1,1,0.5) #do we propose to turn on or off a y.true for this j-k? symmetric with p=0.5
        if(updown==1){ #propose to turn on a y.true. y.true must be 0 and z must be 1
          select.probs.for <- pd[,obs.j[c]]*(1-y.true[,obs.j[c],obs.k[c]])*z
          select.probs.for <- select.probs.for/sum(select.probs.for)
        }else{ #propose to turn off a y.true
          select.probs.for <- (1-pd[,obs.j[c]])*y.true[,obs.j[c],obs.k[c]]*z
          select.probs.for[obs.i[c]] <- 0 # cannot turn off observed guys
          sum.probs.for <- sum(select.probs.for)
          if(sum.probs.for==0){ #no one can be turned off
            skip <- TRUE
          }else{
            select.probs.for <- select.probs.for/sum.probs.for
          }
        }
        if(!skip){ #skip if no one to turn off
          select.cand <- rcat(1,prob=select.probs.for) #this is not a symmetric proposal
          #swap this y.true state. also symmetric
          if(updown==1){
            y.true.cand[select.cand,obs.j[c],obs.k[c]] <- 1
          }else{
            y.true.cand[select.cand,obs.j[c],obs.k[c]] <- 0
          }
          #update observation model likelihood
          ll.y.cand[select.cand,obs.j[c],obs.k[c]] <- 
            dbinom(y.true.cand[select.cand,obs.j[c],obs.k[c]],1,pd[select.cand,obs.j[c]],log=TRUE)
          #update thinning likelihood
          ll.y.obs.cand[obs.k[c]] <- dThin(x=y.obs[1:n.cap,1:J,obs.k[c]],y.true=y.true.cand[1:M,1:J,obs.k[c]],
                                           lambda=lambda[1:M,1:J],
                                           obs.i=obs.i2D[1:n.obs.cells[obs.k[c]],obs.k[c]],
                                           obs.j=obs.j2D[1:n.obs.cells[obs.k[c]],obs.k[c]],
                                           order=order2D[1:n.obs.cells[obs.k[c]],obs.k[c]],
                                           n.cap=n.cap,log=TRUE)
          #get backwards proposal probs
          if(updown==1){
            select.probs.back <- (1-pd[,obs.j[c]])*y.true.cand[,obs.j[c],obs.k[c]]*z
            select.probs.back[obs.i[c]] <- 0 # cannot turn off observed guys
            select.probs.back <- select.probs.back/sum(select.probs.back)
          }else{
            select.probs.back <- pd[,obs.j[c]]*(1-y.true.cand[,obs.j[c],obs.k[c]])*z
            select.probs.back <- select.probs.back/sum(select.probs.back)
          }
          
          logProb.curr <- ll.y.obs[obs.k[c]] +  ll.y[select.cand,obs.j[c],obs.k[c]]
          logProb.cand <- ll.y.obs.cand[obs.k[c]] +  ll.y.cand[select.cand,obs.j[c],obs.k[c]]
          log_MH_ratio <-  (logProb.cand + log(select.probs.back[select.cand])) - (logProb.curr + log(select.probs.for[select.cand]))
          
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[select.cand,obs.j[c],obs.k[c]] <- y.true.cand[select.cand,obs.j[c],obs.k[c]]
            ll.y[select.cand,obs.j[c],obs.k[c]]  <- ll.y.cand[select.cand,obs.j[c],obs.k[c]]
            ll.y.obs[obs.k[c]] <- ll.y.obs.cand[obs.k[c]]
          }else{
            y.true.cand[select.cand,obs.j[c],obs.k[c]] <- y.true[select.cand,obs.j[c],obs.k[c]]
            ll.y.cand[select.cand,obs.j[c],obs.k[c]]  <- ll.y[select.cand,obs.j[c],obs.k[c]]
            ll.y.obs.cand[obs.k[c]] <- ll.y.obs[obs.k[c]]
          }
        }
      }
      
      #update y.true cells with y.obs=0
      for(k in 1:K){ #loop over occasions
        for(i in 1:n.obs.cells[k]){ #loop over captured individuals
          this.i <- obs.i2D[i,k]
          this.j <- obs.j2D[i,k]
          skip <- FALSE
          updown <- rbinom(1,1,0.5) #do we propose to turn on or off a y.true for this j-k? symmetric with p=0.5
          if(updown==1){ #propose to turn on a y.true. y.true must be 0
            select.probs.for <- pd[this.i,]*(1-y.true[this.i,,k])
            select.probs.for <- select.probs.for/sum(select.probs.for)
          }else{ #propose to turn off a y.true
            select.probs.for <- (1-pd[this.i,])*y.true[this.i,,k]
            select.probs.for[this.j] <- 0 # cannot turn off trap where this guy was observed
            sum.probs.for <- sum(select.probs.for)
            if(sum.probs.for==0){ #no one can be turned off
              skip <- TRUE
            }else{
              select.probs.for <- select.probs.for/sum.probs.for
            }
          }
          if(!skip){ #skip if no one to turn off
            select.cand <- rcat(1,prob=select.probs.for) #this is not a symmetric proposal
            #swap this y.true state. also symmetric
            if(updown==1){
              y.true.cand[this.i,select.cand,k] <- 1
            }else{
              y.true.cand[this.i,select.cand,k] <- 0
            }
            #update observation model likelihood
            ll.y.cand[this.i,select.cand,k] <- dbinom(y.true.cand[this.i,select.cand,k],1,pd[this.i,select.cand],log=TRUE)
            #update thinning likelihood
            ll.y.obs.cand[k] <- dThin(x=y.obs[1:n.cap,1:J,k],y.true=y.true.cand[1:M,1:J,k],
                                      lambda=lambda[1:M,1:J],
                                      obs.i=obs.i2D[1:n.obs.cells[k],k],
                                      obs.j=obs.j2D[1:n.obs.cells[k],k],
                                      order=order2D[1:n.obs.cells[k],k],
                                      n.cap=n.cap,log=TRUE)
            #get backwards proposal probs
            if(updown==1){
              select.probs.back <- (1-pd[this.i,])*y.true.cand[this.i,,k]
              select.probs.back[this.j] <- 0 # cannot turn off trap where this guy was observed
              select.probs.back <- select.probs.back/sum(select.probs.back)
            }else{
              select.probs.back <- pd[this.i,]*(1-y.true.cand[this.i,,k])
              select.probs.back <- select.probs.back/sum(select.probs.back)
            }
            
            logProb.curr <- ll.y.obs[k] + ll.y[this.i,select.cand,k]
            logProb.cand <- ll.y.obs.cand[k] + ll.y.cand[this.i,select.cand,k]
            log_MH_ratio <-  (logProb.cand + log(select.probs.back[select.cand])) - (logProb.curr + log(select.probs.for[select.cand]))
            accept <- decide(log_MH_ratio)
            if(accept){
              y.true[this.i,select.cand,k] <- y.true.cand[this.i,select.cand,k]
              ll.y[this.i,select.cand,k]  <- ll.y.cand[this.i,select.cand,k]
              ll.y.obs[k] <- ll.y.obs.cand[k]
            }else{
              y.true.cand[this.i,select.cand,k] <- y.true[this.i,select.cand,k]
              ll.y.cand[this.i,select.cand,k]  <- ll.y[this.i,select.cand,k]
              ll.y.obs.cand[k] <- ll.y.obs[k]
            }
          }
        }
      }
      #now update order
      order2D.cand <- order2D
      for(k in 1:K){
        #symmetric proposal
        select.probs <- rep(1/n.obs.cells[k],n.obs.cells[k])
        swap1 <- rcat(1,prob=select.probs)
        swap2 <- rcat(1,prob=select.probs)
        if(swap1!=swap2){
          order2D.cand[swap1,k] <- order2D[swap2,k]
          order2D.cand[swap2,k] <- order2D[swap1,k]
          ll.y.obs.cand[k] <- dThin(x=y.obs[1:n.cap,1:J,k],y.true=y.true[1:M,1:J,k],
                                    lambda=lambda[1:M,1:J],
                                    obs.i=obs.i2D[1:n.obs.cells[k],k],
                                    obs.j=obs.j2D[1:n.obs.cells[k],k],
                                    order=order2D.cand[1:n.obs.cells[k],k],
                                    n.cap=n.cap,log=TRUE)
          log_MH_ratio <-  ll.y.obs.cand[k] - ll.y.obs[k]
          accept <- decide(log_MH_ratio)
          
          if(accept){
            order2D[swap1,k] <- order2D.cand[swap1,k]
            order2D[swap2,k] <- order2D.cand[swap2,k]
            ll.y.obs[k] <- ll.y.obs.cand[k]
          }else{
            order2D.cand[swap1,k] <- order2D[swap1,k]
            order2D.cand[swap2,k] <- order2D[swap2,k]
            ll.y.obs.cand[k] <- ll.y.obs[k]
          }
        }
      }
    }
    
    #put everything back into the model$stuff 
    model$y.true <<- y.true
    model$order2D <<- order2D
    model$calculate(calcNodes) #update dependencies, likelihoods
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    z.ups <- control$z.ups
    M <- control$M
    #nodes used for update
    y.nodes <- model$expandNodeNames("y.true")
    N.node <- model$expandNodeNames("N")
    z.nodes <- model$expandNodeNames("z")
    pd.nodes <- model$expandNodeNames(paste("pd"))
    lambda.nodes <- model$expandNodeNames(paste("lambda"))
    calcNodes <- c(N.node,z.nodes,pd.nodes,lambda.nodes,y.nodes)
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
        if(sum(model$y.true[pick,,])>0){ #is this individual captured?
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
    model$calculate(lambda.nodes) #not used in this update, but needs to be updated after we are done
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # copy(from = model, to = mvSaved, row = 1, nodes = z.nodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)