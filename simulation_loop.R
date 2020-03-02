# Generate dataframe of simulation parameters
desired_values_1 <- c(100,100,100,100,100,25,25,25,100,100,100,100,100,100,100,100,100,100,25,25,25,25,25,25)
desired_values_2 <- c(100,100,100,100,100,175,175,175,100,100,100,100,100,100,100,100,100,100,175,175,175,175,175,175)
nudiag_1 <- c(0.01,50,50,50,50,0.01,50,50,0.01,0.01,50,50,0.01,50,50,0.01,50,50,0.01,50,50,0.01,50,50)
nudiag_2 <- c(0.01,50,50,50,50,0.01,50,50,0.01,0.01,50,50,0.01,50,50,0.01,50,50,0.01,50,50,0.01,50,50)
rho_1 <- c(0,400,400,400,400,0,400,400,0,0,400,400,0,400,400,0,400,400,0,400,400,0,400,400)
rho_2 <- c(0,360,-360,360,-360,0,360,-360,0,0,360,-360,0,360,-360,0,360,-360,0,360,-360,0,360,-360)
rho_3 <- c(0,360,-360,360,-360,0,360,-360,0,0,360,-360,0,360,-360,0,360,-360,0,360,-360,0,360,-360)
rho_4 <- c(0,400,400,400,400,0,400,400,0,0,400,400,0,400,400,0,400,400,0,400,400,0,400,400)
plantnames_1 <- c("Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag",
                  "Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag","Astrag",
                  "Astrag","Astrag")
plantnames_2 <- c("Lupine","Lupine","Lupine","Penst","Penst","Lupine","Lupine","Lupine","Penst","Lupine","Lupine","Lupine",
                  "Lupine","Lupine","Lupine","Lupine","Lupine","Lupine","Lupine","Lupine","Lupine","Lupine","Lupine",
                  "Lupine")
new_dist <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0,0,3,3,3,-7,-7,-7,0,0,0,-7,-7,-7)

params <- data.frame(desired_values_1, desired_values_2, nudiag_1, nudiag_2, rho_1, rho_2, rho_3, rho_4, plantnames_1,
                     plantnames_2, new_dist, stringsAsFactors = F)
params <- params[-which(params$new_dist==3),]

patch.rad <- 20

load('/Users/JacobSocolar/Dropbox/Work/Bee_Path/bombus.new.352.R')

pt_amts <- .2*.8^(0:19)
mcs <- bombus.new.352$coefficients
msd_dist <- c(17.98465, 11.59513)

run_sims <- function(params, i){
  outputs <- list()
  for(j in 1:20){
    check1 <- check2 <- check3 <- T
    while(check1 | check2 | check3){
      patch <- patchsim(nudiag = c(params$nudiag_1[i], params$nudiag_2[i]),
                        rho = matrix(c(params$rho_1[i], params$rho_2[i], params$rho_3[i], params$rho_4[i]), nrow=2, byrow = T),
                        plantnames = c(params$plantnames_1[i], params$plantnames_2[i]),
                        desired_values = c(50^2/(3.14159*patch.rad^2)*params$desired_values_1[i], 50^2/(3.14159*patch.rad^2)*params$desired_values_2[i]),
                        radius = patch.rad)
      check1 <- abs(log(sum(patch$fp$Species == params$plantnames_1[i])/params$desired_values_1[i])) > .356
      check2 <- abs(log(sum(patch$fp$Species == params$plantnames_2[i])/params$desired_values_2[i])) > .356
      rabr <- sum(patch$fp$Species == params$plantnames_1[i])/sum(patch$fp$Species == params$plantnames_2[i])
      tabr <- params$desired_values_1[i]/params$desired_values_2[i]
      check3 <- abs(log(rabr/tabr)) > .105
    }
    trans_probs <- transprobs(fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs, msd_dist = msd_dist)
    mcs2 <- mcs
    if(!is.na(params$new_dist[i])){
      mcs2["distance"] <- params$new_dist[i]
    }
    paths <- pathsim(nsim = 1000, nstep = 10, fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs2, pollenvector = pt_amts, start = list(patch$pr,"none"), msd_dist = msd_dist, trans_probs = trans_probs, pt_type = "amt")
    outputs[[j]] <- paths
  }
  return(outputs)
}

library(doParallel)

simulations <- list()

s <- proc.time()                          # for keeping track of how long the job takes
for(rowi in 1:nrow(params)){
  print(rowi)
  cl <- parallel::makeCluster(6)            # make a cluster (called cl) with 6 cores. Note that this function comes from the package 'parallel' on which 'doParallel' depends
  doParallel::registerDoParallel(cl)        # some doParallel voodoo to get the cores ready for work
  clusterCall(cl, function() library(RandomFields))
  simulations[[rowi]] <- foreach(k = 1:6) %dopar% run_sims(params = params, i = rowi)
  parallel::stopCluster(cl)                # Closing down our cluster
}
elapsed <- proc.time() - s               # how long the job took

save(simulations, file = "/Users/JacobSocolar/Dropbox/Work/Bee_Path/full_simulations.Rdata")



  