transprobs <- function(fp, fp_dist, mcs, msd_dist){
  ##### Get a matrix of the transition probabilities #####
  nplant <- nrow(fp)
  trans_probs <- array(data = 0, dim = c(nplant, nplant, nplant))
  for(i in 1:nplant){ # i indexes plant 0
    for(j in 1:nplant){ # j indexes plant 1
      if(i != j){
        rpk <- rep(0, nplant)
        for(k in 1:nplant){ # k indexes plant 2
          if(j != k){
            constancy <- as.numeric(fp$Species[j] == fp$Species[k]) - as.numeric(fp$Species[j] != fp$Species[k])
            distance <- (fp_dist[j,k]-msd_dist[1])/msd_dist[2] #mean and sd from "ALL paths 2015 2016 w constancy lambda.csv"
            spCast <- fp$Species[k] == "Castilleja"
            spPenst <- fp$Species[k] == "Penst"
            spLupine <- fp$Species[k] == "Lupinus"
            lambda <- fp$lambda[k]
            v1 <- c(fp$X[j] - fp$X[i], fp$Y[j] - fp$Y[i])
            v2 <- c(fp$X[k] - fp$X[j], fp$Y[k] - fp$Y[j])
            dot.prod <- v1 %*% v2
            cos.turn <- dot.prod/(norm(v1,type="2") * norm(v2,type="2"))
            
            rpk[k] <- exp(mcs['constancy']*constancy + mcs['distance']*distance + 
                            mcs['end.plant.spCast']*spCast + mcs['end.plant.spLupine']*spLupine +
                            mcs['end.plant.spPenst']*spPenst + mcs['lambda']*lambda +
                            mcs['cos.turn']*cos.turn + mcs['constancy:lambda']*constancy*lambda)
          }
          apk <- rpk/sum(rpk)
          trans_probs[i, j, ] <- apk
        }
      }
    }
  }
  return(trans_probs)
}


pathsim <- function(nsim, nstep, fp, fp_dist, mcs, pollenvector, start, msd_dist, trans_probs){
  # nsim: number of simulations to run
  # nstep: number of steps per simulation, either a single integer or an integer vector of length nsim
  # fp: a dataframe with columns named "Species", "X", "Y", "lambda", "uniqueID", and "prop".  Other columns will be ignored.
  # mcs: named fitted model coefficients
  # start: 'random plant' initializes the path on a random plant in fp.
  #         an integer initializes the path on the plant corresponding to that row of fp
  #         passing a specific list structure initializes the path at a random point on the patch edge.
  # pollenvector: a vector of pollen transfer probabilities. pollen_vector[i] gives the pollen
  #     transfer probability between plants with i-1 plants visited in between.
  # return probs: should the output include the matrix of plant-to-plant transition probabilities?
  # msd: mean and standard deviation of distances from the data object to which the model that yields mcs was fit
  
  nplant <- nrow(fp)
  if(length(nstep) == 1){nstep <- rep(nstep, nsim)}
  if(length(pollenvector) < max(nstep)){
    pollenvector <- c(pollenvector, rep(0, max(nstep) - length(pollenvector)))
  }
  
  ##### Simulate the actual paths #####
  paths <- matrix(NA, nrow = max(nsim), ncol = max(nstep) + 1)
  # Below, cp1 is the first plant in a path; cp2 is the second plant
  for(sim in 1:nsim){
    if(!is.list(start)){
      if(start == 'random plant'){
        cp1 <- sample(1:nplant, size=1)
      }else if(is.numeric(start)){
        cp1 <- start
      }
      XD <- fp$X[cp1]
      YD <- fp$Y[cp1]
      spD <- fp$Species[cp1]
      
      rel_probs_2 <- vector()
      for(j in 1:nplant){
        constancy <- as.numeric(spD == fp$Species[j]) - as.numeric(spD != fp$Species[j])
        d1 <- ((XD - fp$X[j])^2 + (YD - fp$Y[j])^2)^.5
        distance <- (d1-msd_dist[1])/msd_dist[2] 
        spCast <- fp$Species[j] == "Castilleja"
        spPenst <- fp$Species[j] == "Penst"
        spLupine <- fp$Species[j] == "Lupinus"
        lambda <- fp$lambda[j]
        cos.turn <- 0
        
        rel_probs_2[j] <- exp(mcs['constancy']*constancy + mcs['distance']*distance + 
                                mcs['end.plant.spCast']*spCast + mcs['end.plant.spLupine']*spLupine +
                                mcs['end.plant.spPenst']*spPenst + mcs['lambda']*lambda + 
                                mcs['cos.turn']*cos.turn + mcs['constancy:lambda']*constancy*lambda)
      }
      start_probs <- rel_probs_2/sum(rel_probs_2)
      cp2 <- sample.int(nplant, 1, prob = start_probs)
      
      
    }
    else{
      theta <- runif(1,0,2*pi)
      XD <- start[[1]]*cos(theta)+25
      YD <- start[[1]]*sin(theta)+25
      spD <- start[[2]]

      rel_probs_1 <- vector()
      for(j in 1:nplant){
        constancy <- as.numeric(spD == fp$Species[j]) - as.numeric(spD != fp$Species[j])
        if(spD == "none"){constancy <- 0}
        d1 <- ((XD - fp$X[j])^2 + (YD - fp$Y[j])^2)^.5
        distance <- (d1-msd_dist[1])/msd_dist[2] 
        spCast <- fp$Species[j] == "Castilleja"
        spPenst <- fp$Species[j] == "Penst"
        spLupine <- fp$Species[j] == "Lupinus"
        lambda <- fp$lambda[j]
        cos.turn <- 0
        
        rel_probs_1[j] <- exp(mcs['constancy']*constancy + mcs['distance']*distance + 
                                  mcs['end.plant.spCast']*spCast + mcs['end.plant.spLupine']*spLupine +
                                  mcs['end.plant.spPenst']*spPenst + mcs['lambda']*lambda + 
                                  mcs['cos.turn']*cos.turn + mcs['constancy:lambda']*constancy*lambda)
      }
      start_probs <- rel_probs_1/sum(rel_probs_1)
      cp1 <- sample.int(nplant, 1, prob = start_probs)
      
      
      XD1 <- fp$X[cp1]
      YD1 <- fp$Y[cp1]
      spD1 <- fp$Species[cp1]
      
      rel_probs_2 <- vector()
      for(j in 1:nplant){
        constancy <- as.numeric(spD1 == fp$Species[j]) - as.numeric(spD1 != fp$Species[j])
        d1 <- ((XD1 - fp$X[j])^2 + (YD1 - fp$Y[j])^2)^.5
        distance <- (d1-msd_dist[1])/msd_dist[2] 
        spCast <- fp$Species[j] == "Castilleja"
        spPenst <- fp$Species[j] == "Penst"
        spLupine <- fp$Species[j] == "Lupinus"
        lambda <- fp$lambda[j]
        v1 <- c(XD1 - XD, YD1 - YD)
        v2 <- c(fp$X[j] - XD1, fp$Y[j] - YD1)
        dot.prod <- v1 %*% v2
        cos.turn <- dot.prod/(norm(v1,type="2") * norm(v2,type="2"))
        
        rel_probs_2[j] <- exp(mcs['constancy']*constancy + mcs['distance']*distance + 
                                mcs['end.plant.spCast']*spCast + mcs['end.plant.spLupine']*spLupine +
                                mcs['end.plant.spPenst']*spPenst + mcs['lambda']*lambda + 
                                mcs['cos.turn']*cos.turn + mcs['constancy:lambda']*constancy*lambda)
      }
      rel_probs_2[cp1] <- 0
      start_probs <- rel_probs_2/sum(rel_probs_2)
      cp2 <- sample.int(nplant, 1, prob = start_probs)
    }
    
    cp <- c(cp1, cp2)

    ids <- rep(NA, max(nstep)+1)
    ids[1] <- cp1
    ids[2] <- cp2
    for(i in 3:(nstep[sim]+1)){
      c <- sample.int(nplant, 1, prob = trans_probs[cp[1], cp[2], ])
      cp[1] <- cp[2]
      ids[i] <- cp[2] <- c
    }
    paths[sim, ] <- ids
  }
  pollen <- fp[, c("Species", "uniqueID")]
  pollen$success <- 0
  for(i in 1:nsim){
    for(p1 in 1:nstep[i]){
      for(p2 in (p1+1):(nstep[i]+1)){
        if((paths[i, p1] != paths[i, p2])){
          if(fp$Species[paths[i, p1]] == fp$Species[paths[i, p2]]){
            pollen$success[paths[i, p2]] <- pollen$success[paths[i, p2]] + rbinom(1,1,pollenvector[p2 - p1])
          }
        }
      }
    }
  }
  return(list(paths=paths, pollen=pollen, trans_probs=trans_probs))
}