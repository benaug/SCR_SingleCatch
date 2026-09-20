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

#collapse repeated Bernoulli likelihood across occasions using K1D.p/K1D.c
dBernoulliMatrixMb <- nimbleFunction(
  run = function(x = double(2), kern = double(1), p0.p = double(0), p0.c = double(0),
                 y.state = double(2), K2D = double(2), K1D.p = double(1), K1D.c = double(1),
                 z = double(0), log = integer(0)) {
    returnType(double(0))
    if(z==0){#skip calculation if z=0
      return(0)
    }else{
      J <- nimDim(K2D)[1]
      K <- nimDim(K2D)[2]
      logProb <- 0
      for(j in 1:J){
        ncap.p <- 0
        ncap.c <- 0
        for(k in 1:K){
          if(x[j,k]==1){
            if(K2D[j,k]==0){
              return(-Inf)
            }
            if(y.state[j,k]==0){
              ncap.p <- ncap.p+1
            }else{
              ncap.c <- ncap.c+1
            }
          }
        }
        pd.p <- p0.p*kern[j]
        pd.c <- p0.c*kern[j]
        if(ncap.p>0){
          if(pd.p<=0){
            return(-Inf)
          }
          logProb <- logProb+ncap.p*log(pd.p)
        }
        nnon.p <- K1D.p[j]-ncap.p
        if(nnon.p>0){
          logProb <- logProb+nnon.p*log1p(-pd.p)
        }
        if(ncap.c>0){
          if(pd.c<=0){
            return(-Inf)
          }
          logProb <- logProb+ncap.c*log(pd.c)
        }
        nnon.c <- K1D.c[j]-ncap.c
        if(nnon.c>0){
          logProb <- logProb+nnon.c*log1p(-pd.c)
        }
      }
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBernoulliMatrixMb <- nimbleFunction(
  run = function(n = integer(0), kern = double(1), p0.p = double(0), p0.c = double(0),
                 y.state = double(2), K2D = double(2), K1D.p = double(1), K1D.c = double(1),
                 z = double(0)) {
    returnType(double(2))
    J <- nimDim(K2D)[1]
    K <- nimDim(K2D)[2]
    out <- matrix(0,J,K)
    return(out)
  }
)

#used in pSmaller() below
integrand <- nimbleFunction(
  run = function(x = double(1), param = double(1)){ 
    returnType(double(1))
    # param contains:
    # lambda1, lambda2[1:n], exp(-lambda2[1:n])
    n <- (length(param) - 1) / 2
    lambda1 <- param[1]
    #can split vector, but faster to not split it
    #lambda2 <- param[2:(n+1)]
    #exp.lambda2 <- param[(n+2):(2*n+1)]
    n.x <- length(x)
    prod.term <- exp(-lambda1 * x)
    for(j in 1:n.x){
      for(i in 1:n){
        #if vector split
        # prod.term[j] <- prod.term[j] * (exp(-lambda2[i] * x[j]) - exp.lambda2[i])
        #if not split
        prod.term[j] <- prod.term[j] * (exp(-param[i+1] * x[j]) - param[n+i+1])
      }
    }
    return(prod.term)
  })

#probability exponential RV right-truncated at 1 with parameter lambda1 is less than
#one or more other exponential RVs right-truncated at 1 with parameter(s) lambda2
pSmaller <- nimbleFunction(
  run = function(lambda1 = double(0), lambda2 = double(1), log = integer(0)) {
    returnType(double(0))
    if(lambda1 < 1e-8){
      lambda1 <- 1e-8
    }
    for(i in 1:length(lambda2)){
      if(lambda2[i] < 1e-8){
        lambda2[i] <- 1e-8
      }
    }
    exp.lambda2 <- exp(-lambda2)
    # pass precomputed exp(-lambda2) to integrand
    param <- c(lambda1,lambda2,exp.lambda2)
    #integral from 0 to 1
    integral <- nimIntegrate(integrand, lower = 0, upper = 1, param = param)[1]
    #denominator terms
    lambda.term <- 1 - exp(-lambda1)
    lambda2.prod <- prod(1 - exp.lambda2)
    logProb <- log(lambda1*integral) - log(lambda.term * lambda2.prod)
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  })


dThin <- nimbleFunction(
  run = function(x = double(2), y.true = double(2), y.state = double(2), lambda.p = double(2), lambda.c = double(2), 
                 order = double(1), obs.i = double(1), obs.j = double(1), n.obs = double(0), n.cap = double(0), log = integer(0)) { 
    returnType(double(0))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    #track availability
    i.available <- rep(1,M)
    j.available <- rep(1,J)
    #identify number of latent capture events
    n.latent <- 0
    for(i in 1:M){
      for(j in 1:J){
        if(y.true[i,j]==1){
          n.latent <- n.latent+1
        }
      }
    }
    #if no observed captures, any latent event would necessarily generate at least one real capture
    if(n.obs==0){ 
      if(n.latent>0){ 
        return(-Inf) 
      }else{ 
        return(0) 
      } 
    } 
    #with >=1 observed capture, there must be >=1 latent event
    if(n.latent==0){
      return(-Inf)
    }
    #identify indices of latent capture events once
    latent.i <- rep(0,n.latent)
    latent.j <- rep(0,n.latent)
    idx.latent <- 1
    
    for(i in 1:M){
      for(j in 1:J){
        if(y.true[i,j]==1){
          latent.i[idx.latent] <- i
          latent.j[idx.latent] <- j
          idx.latent <- idx.latent+1
        }
      }
    }
    #construct inverse capture order
    order.idx <- rep(0,n.obs) 
    for(ii in 1:n.obs){ 
      order.idx[order[ii]] <- ii
    }
    other.lambdas <- rep(0,n.latent)
    logProb <- 0
    
    for(o in 1:n.obs){ 
      idx <- order.idx[o]
      focal.i <- obs.i[idx]
      focal.j <- obs.j[idx]
      #observed capture must correspond to a latent capture
      if(y.true[focal.i,focal.j]==0){
        return(-Inf)
      }
      #observed individual and trap must still be available
      if(i.available[focal.i]==0){
        return(-Inf)
      }
      if(j.available[focal.j]==0){
        return(-Inf)
      }
      #choose first/subsequent-capture lambda from observed y.state
      if(y.state[focal.i,focal.j]==0){
        focal.lambda <- lambda.p[focal.i,focal.j]
      }else{
        focal.lambda <- lambda.c[focal.i,focal.j]
      }
      #collect all other currently competing latent events
      n.other.lambdas <- 0
      for(l in 1:n.latent){
        i <- latent.i[l]
        j <- latent.j[l]
        if(i.available[i]==1){
          if(j.available[j]==1){
            if(!(i==focal.i & j==focal.j)){
              n.other.lambdas <- n.other.lambdas+1
              if(y.state[i,j]==0){
                other.lambdas[n.other.lambdas] <- lambda.p[i,j]
              }else{
                other.lambdas[n.other.lambdas] <- lambda.c[i,j]
              }
            }
          }
        }
      }
      #probability focal event occurs before all competitors
      if(n.other.lambdas>0){
        logProb <- logProb + pSmaller(focal.lambda,other.lambdas[1:n.other.lambdas],log=TRUE)
      }
      #single-catch: individual and trap are unavailable afterward
      i.available[focal.i] <- 0
      j.available[focal.j] <- 0
    }
    
    #after the final observed capture, there cannot be
    #another latent event with both its individual and trap available
    for(l in 1:n.latent){
      i <- latent.i[l]
      j <- latent.j[l]
      if(i.available[i]==1){
        if(j.available[j]==1){
          return(-Inf)
        }
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rThin <- nimbleFunction(
  run = function(n = integer(0), y.true = double(2), y.state = double(2),
                 lambda.p = double(2),lambda.c = double(2), obs.i = double(1),
                 obs.j = double(1), order = double(1), n.obs = double(0), n.cap = double(0)){
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
    n.obs.cells.max <- control$n.obs.cells.max
    n.obs.cells.all <- sum(n.obs.cells)
    K2D <- control$K2D
    y.obs <- control$y.obs
    n.cap <- control$n.cap
    calcNodes <- model$getDependencies(c("y.true","order2D","y.obs"))
  },
  run = function(){
    y.true <- model$y.true
    y.state <- model$y.state
    z <- model$z
    kern <- model$kern 
    pd.p <- model$p0.p[1]*kern 
    pd.c <- model$p0.c[1]*kern 
    lambda.p <- model$lambda.p
    lambda.c <- model$lambda.c
    order2D <- model$order2D
    
    ll.y.obs <- model$logProb_y.obs[1,1,]
    
    for(up in 1:y.ups){ #update one or more times per iteration
      # update y.true for cells with y.obs=1
      for(c in 1:n.obs.cells.all){
        skip <- FALSE
        updown <- rbinom(1,1,0.5)
        this.j <- obs.j[c] 
        this.k <- obs.k[c] 
        pd.use <- pd.p[,this.j]*(1-y.state[,this.j,this.k]) +
          pd.c[,this.j]*y.state[,this.j,this.k] 
        
        if(updown==1){ #propose to turn on a y.true. y.true must be 0 and z must be 1
          select.probs.for <- pd.use*(1-y.true[,this.j,this.k])*z 
          sum.probs.for <- sum(select.probs.for) 
          if(sum.probs.for==0){ 
            skip <- TRUE 
          }else{
            select.probs.for <- select.probs.for/sum.probs.for 
          }
        }else{ #propose to turn off a y.true
          select.probs.for <- (1-pd.use)*y.true[,this.j,this.k]*z 
          select.probs.for[obs.i[c]] <- 0 # cannot turn off observed guys
          sum.probs.for <- sum(select.probs.for)
          if(sum.probs.for==0){
            skip <- TRUE
          }else{
            select.probs.for <- select.probs.for/sum.probs.for
          }
        }
        
        if(!skip){
          select.cand <- rcat(1,prob=select.probs.for)
          y.curr <- y.true[select.cand,this.j,this.k] 
          
          this.pd <- pd.use[select.cand] 
          if(y.curr==1){
            ll.y.curr <- log(this.pd) 
            y.true[select.cand,this.j,this.k] <- 0 
            ll.y.cand <- log1p(-this.pd) 
          }else{
            ll.y.curr <- log1p(-this.pd) 
            y.true[select.cand,this.j,this.k] <- 1 
            ll.y.cand <- log(this.pd) 
          }
          
          # update thinning likelihood
          ll.y.obs.cand <- dThin(x=y.obs[1:n.cap,1:J,this.k],y.true=y.true[1:M,1:J,this.k], 
                                 y.state=y.state[1:M,1:J,this.k],
                                 lambda.p=lambda.p[1:M,1:J],lambda.c=lambda.c[1:M,1:J],
                                 obs.i=obs.i2D[1:n.obs.cells.max,this.k],
                                 obs.j=obs.j2D[1:n.obs.cells.max,this.k],
                                 order=order2D[1:n.obs.cells.max,this.k],
                                 n.obs=n.obs.cells[this.k],n.cap=n.cap,log=TRUE) 
          
          #get backwards proposal probs
          if(updown==1){
            select.probs.back <- (1-pd.use)*y.true[,this.j,this.k]*z 
            select.probs.back[obs.i[c]] <- 0
            select.probs.back <- select.probs.back/sum(select.probs.back)
          }else{
            select.probs.back <- pd.use*(1-y.true[,this.j,this.k])*z 
            select.probs.back <- select.probs.back/sum(select.probs.back)
          }
          
          logProb.curr <- ll.y.obs[this.k]+ll.y.curr 
          logProb.cand <- ll.y.obs.cand+ll.y.cand 
          log_MH_ratio <- (logProb.cand+log(select.probs.back[select.cand])) -
            (logProb.curr+log(select.probs.for[select.cand]))
          
          accept <- decide(log_MH_ratio)
          if(accept){
            ll.y.obs[this.k] <- ll.y.obs.cand 
          }else{
            y.true[select.cand,this.j,this.k] <- y.curr #restore only changed cell
          }
        }
      }
      
      #update y.true cells with y.obs=0
      for(k in 1:K){
        if(n.obs.cells[k]>0){ 
          for(i in 1:n.obs.cells[k]){ 
            this.i <- obs.i2D[i,k]
            this.j <- obs.j2D[i,k]
            skip <- FALSE
            updown <- rbinom(1,1,0.5)
            
            #proposal probabilities use the observed behavioral state for each trap
            pd.use <- pd.p[this.i,]*(1-y.state[this.i,,k]) +
              pd.c[this.i,]*y.state[this.i,,k] 
            
            if(updown==1){ #propose to turn on a y.true. y.true must be 0
              select.probs.for <- pd.use*(1-y.true[this.i,,k])*K2D[,k] #state-specific and exclude closed traps
              sum.probs.for <- sum(select.probs.for) 
              if(sum.probs.for==0){ 
                skip <- TRUE 
              }else{
                select.probs.for <- select.probs.for/sum.probs.for 
              }
            }else{ #propose to turn off a y.true
              select.probs.for <- (1-pd.use)*y.true[this.i,,k] 
              select.probs.for[this.j] <- 0 # cannot turn off trap where this guy was observed
              sum.probs.for <- sum(select.probs.for)
              if(sum.probs.for==0){
                skip <- TRUE
              }else{
                select.probs.for <- select.probs.for/sum.probs.for
              }
            }
            
            if(!skip){
              select.cand <- rcat(1,prob=select.probs.for)
              y.curr <- y.true[this.i,select.cand,k] 
              
              this.pd <- pd.use[select.cand] 
              if(y.curr==1){
                ll.y.curr <- log(this.pd) 
                y.true[this.i,select.cand,k] <- 0 
                ll.y.cand <- log1p(-this.pd) 
              }else{
                ll.y.curr <- log1p(-this.pd) 
                y.true[this.i,select.cand,k] <- 1 
                ll.y.cand <- log(this.pd) 
              }
              
              #update thinning likelihood
              ll.y.obs.cand <- dThin(x=y.obs[1:n.cap,1:J,k],y.true=y.true[1:M,1:J,k], 
                                     y.state=y.state[1:M,1:J,k],
                                     lambda.p=lambda.p[1:M,1:J],lambda.c=lambda.c[1:M,1:J],
                                     obs.i=obs.i2D[1:n.obs.cells.max,k],
                                     obs.j=obs.j2D[1:n.obs.cells.max,k],
                                     order=order2D[1:n.obs.cells.max,k],
                                     n.obs=n.obs.cells[k],n.cap=n.cap,log=TRUE) 
              
              #get backwards proposal probs
              if(updown==1){
                select.probs.back <- (1-pd.use)*y.true[this.i,,k] 
                select.probs.back[this.j] <- 0
                select.probs.back <- select.probs.back/sum(select.probs.back)
              }else{
                select.probs.back <- pd.use*(1-y.true[this.i,,k])*K2D[,k] 
                select.probs.back <- select.probs.back/sum(select.probs.back)
              }
              
              logProb.curr <- ll.y.obs[k]+ll.y.curr 
              logProb.cand <- ll.y.obs.cand+ll.y.cand 
              log_MH_ratio <- (logProb.cand+log(select.probs.back[select.cand])) -
                (logProb.curr+log(select.probs.for[select.cand]))
              
              accept <- decide(log_MH_ratio)
              if(accept){
                ll.y.obs[k] <- ll.y.obs.cand 
              }else{
                y.true[this.i,select.cand,k] <- y.curr #restore only changed cell
              }
            }
          }
        } 
      }
      
      #now update order
      for(k in 1:K){
        if(n.obs.cells[k]>1){ 
          #symmetric proposal
          select.probs <- rep(1/n.obs.cells[k],n.obs.cells[k])
          swap1 <- rcat(1,prob=select.probs) 
          swap2 <- rcat(1,prob=select.probs) 
          if(swap1!=swap2){
            old1 <- order2D[swap1,k] 
            old2 <- order2D[swap2,k] 
            order2D[swap1,k] <- old2 
            order2D[swap2,k] <- old1 
            
            ll.y.obs.cand <- dThin(x=y.obs[1:n.cap,1:J,k],y.true=y.true[1:M,1:J,k], 
                                   y.state=y.state[1:M,1:J,k],
                                   lambda.p=lambda.p[1:M,1:J],lambda.c=lambda.c[1:M,1:J],
                                   obs.i=obs.i2D[1:n.obs.cells.max,k],
                                   obs.j=obs.j2D[1:n.obs.cells.max,k],
                                   order=order2D[1:n.obs.cells.max,k],
                                   n.obs=n.obs.cells[k],n.cap=n.cap,log=TRUE) 
            log_MH_ratio <- ll.y.obs.cand-ll.y.obs[k]
            accept <- decide(log_MH_ratio)
            
            if(accept){
              ll.y.obs[k] <- ll.y.obs.cand
            }else{
              order2D[swap1,k] <- old1 #restore only swapped elements
              order2D[swap2,k] <- old2 
            }
          }
        } 
      }
    }
    
    model$y.true <<- y.true
    model$order2D <<- order2D
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    z.ups <- control$z.ups
    M <- control$M
    #nodes used for update, calcNodes + z nodes
    y.nodes <- model$expandNodeNames("y.true")
    N.node <- model$expandNodeNames("N")
    z.nodes <- model$expandNodeNames("z")
    kern.nodes <- model$expandNodeNames(paste("kern"))
    lambda.p.nodes <- model$expandNodeNames(paste("lambda.p")) 
    lambda.c.nodes <- model$expandNodeNames(paste("lambda.c")) 
    calcNodes <- c(N.node,z.nodes,kern.nodes,lambda.p.nodes,lambda.c.nodes,y.nodes) 
  },
  run = function(){
    #build eligible on/off lists once and update them after accepted moves
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(model$z[i]==1){
        if(sum(model$y.true[i,,])==0){ #active individuals with latent captures cannot be turned off
          non.curr <- non.curr+1
          z.on[non.curr] <- i
        }
      }else{
        #if z=0, y.true must be all zero, so no need to scan y.true
        noff.curr <- noff.curr+1
        z.off[noff.curr] <- i
      }
    }
    
    for(up in 1:z.ups){
      updown <- rbinom(1,1,0.5)
      if(updown==0){#subtract
        non.init <- non.curr 
        if(non.init>0){ 
          pick.pos <- rcat(1,rep(1/non.init,non.init)) 
          pick <- z.on[pick.pos] 
          N.init <- model$N[1] 
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<- model$N[1]-1
          model$z[pick] <<- 0
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- 0 
          
          #MH step
          log_MH_ratio <- (lp.proposed.N+lp.proposed.y)-
            (lp.initial.N+lp.initial.y)+log(non.init/N.init) 
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            model$calculate(kern.nodes[pick]) #synchronize after acceptance
            model$calculate(lambda.p.nodes[pick]) 
            model$calculate(lambda.c.nodes[pick]) 
            model$calculate(y.nodes[pick]) 
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["kern",1][pick,] <<- model[["kern"]][pick,] 
            mvSaved["lambda.p",1][pick,] <<- model[["lambda.p"]][pick,] 
            mvSaved["lambda.c",1][pick,] <<- model[["lambda.c"]][pick,] 
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            
            #move accepted individual from on list to off list
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(N.node)
          }
        }
        
      }else{#add
        noff.init <- noff.curr 
        if(noff.init>0){ 
          pick.pos <- rcat(1,rep(1/noff.init,noff.init)) 
          pick <- z.off[pick.pos] 
          N.init <- model$N[1] 
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- 0 
          
          #propose new N/z
          model$N[1] <<- model$N[1]+1
          model$z[pick] <<- 1
          
          model$calculate(kern.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N+lp.proposed.y)-
            (lp.initial.N+lp.initial.y)+log((N.init+1)/(non.curr+1)) 
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            model$calculate(lambda.p.nodes[pick]) #synchronize after acceptance
            model$calculate(lambda.c.nodes[pick]) 
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["kern",1][pick,] <<- model[["kern"]][pick,] 
            mvSaved["lambda.p",1][pick,] <<- model[["lambda.p"]][pick,] 
            mvSaved["lambda.c",1][pick,] <<- model[["lambda.c"]][pick,] 
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
            model$calculate(y.nodes[pick]) #restore y logProb
            model$calculate(N.node)
          }
        }
      }
    }
    
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
