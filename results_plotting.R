load("/Users/JacobSocolar/Dropbox/Work/Bee_Path/full_simulations.Rdata")
load('/Users/JacobSocolar/Dropbox/Work/Bee_Path/bombus.new.352.R')

##### simulation parameters (code block copied from simulation_loop.R) #####
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

pt_amts <- .2*.8^(0:19)
mcs <- bombus.new.352$coefficients
msd_dist <- c(17.98465, 11.59513)


##### Understanding the simulation output #####
class(simulations)
length(simulations)  # One element for each parameter combo

class(simulations[[1]])
length(simulations[[1]])  # Six times...

class(simulations[[1]][[1]])
length(simulations[[1]][[1]])  # ... times twenty replicates = 120 replicates per parameter combo

class(simulations[[1]][[1]][[1]])
length(simulations[[1]][[1]][[1]])  # pathsim output

simulations[[1]][[1]][[1]][[1]]  # matrix of what individuals are visited in each path
simulations[[1]][[1]][[1]][[2]]  # data.frame of amt of pollination success by individual
simulations[[1]][[1]][[1]][[3]]  # list of simulation parameters

##### function summarize delivery for a given parameter combination #####
total_pollen <- function(i, species){
  s <- simulations[[i]]
  v <- data.frame(mean = rep(NA,120), var = rep(NA, 120), skew = rep(NA, 120))
  counter <- 0
  for(j in 1:6){
    for(k in 1:20){
      counter <- counter + 1
      rep <- s[[j]][[k]][[2]]
      repSp <- rep[rep$Species == species, ]
      v$mean[counter] <- mean(repSp$success)
      v$var[counter] <- var(repSp$success)
      v$skew[counter] <- moments::skewness(repSp$success)
    }
  }
  return(v)
}

p2 <- params
p2$new_dist[is.na(p2$new_dist)] <- -3.65

sim_summary_1 <- sim_summary_2 <- list()
for(i in 1:nrow(p2)){
  sim_summary_1[[i]] <- cbind(data.frame(sp = rep(1,120), sp2 = rep(p2$plantnames_2[i], 120), abun = rep(p2$desired_values_2[i], 120),
                                   n1 = rep(p2$nudiag_1[i], 120), r2 = rep(p2$rho_2[i], 120)), dist = p2$new_dist[i], total_pollen(i, "Astrag"))
  sim_summary_2[[i]] <- cbind(data.frame(sp = rep(2,120), sp2 = rep(p2$plantnames_2[i], 120), abun = rep(p2$desired_values_2[i], 120),
                                         n1 = rep(p2$nudiag_1[i], 120), r2 = rep(p2$rho_2[i], 120)), dist = p2$new_dist[i], total_pollen(i, p2$plantnames_2[i]))
}


s <- rbind(do.call("rbind", sim_summary_1), do.call("rbind", sim_summary_2))
s$cv <- s$var/s$mean
s$nr <- as.factor(paste(s$n1, s$r2, sep = "_"))
s$cov <- factor(s$nr, levels = c("50_-360", "0.01_0", "50_360"), labels = c("-", "0", "+"))

s1 <- s[s$sp == 1 & s$sp2 == "Lupine" & s$dist == -3.65, ]
boxplot(mean ~ cov*abun, data = s1, notch = T, col = 'goldenrod', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov*abun, data = s1, notch = T)
#boxplot(cv ~ cov*abun, data = s1, notch = T)
#boxplot(skew ~ cov*abun, data = s1, notch = T)

s2 <- s[s$sp == 2 & s$sp2 == "Lupine" & s$dist == -3.65, ]
boxplot(mean ~ cov*abun, data = s2, notch = T, col = 'darkmagenta', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov*abun, data = s2, notch = T)
#boxplot(cv ~ cov*abun, data = s2, notch = T)
#boxplot(skew ~ cov*abun, data = s2, notch = T)

s3 <- s[s$sp == 1 & s$sp2 == "Penst" & s$dist == -3.65, ]
boxplot(mean ~ cov, data = s3, notch = T, col = 'goldenrod', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov, data = s3, notch = T)
#boxplot(cv ~ cov, data = s3, notch = T)
#boxplot(skew ~ cov, data = s3, notch = T)

s4 <- s[s$sp == 2 & s$sp2 == "Penst" & s$dist == -3.65, ]
boxplot(mean ~ cov, data = s4, notch = T, col = 'mediumpurple3', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov, data = s4, notch = T)
#boxplot(cv ~ cov, data = s4, notch = T)
#boxplot(skew ~ cov, data = s4, notch = T)








s1 <- s[s$sp == 1 & s$sp2 == "Lupine" & s$dist == 0, ]
boxplot(mean ~ cov*abun, data = s1, notch = T, col = 'goldenrod', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov*abun, data = s1, notch = T)
#boxplot(cv ~ cov*abun, data = s1, notch = T)
#boxplot(skew ~ cov*abun, data = s1, notch = T)

s2 <- s[s$sp == 2 & s$sp2 == "Lupine" & s$dist == 0, ]
boxplot(mean ~ cov*abun, data = s2, notch = T, col = 'darkmagenta', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov*abun, data = s2, notch = T)
#boxplot(cv ~ cov*abun, data = s2, notch = T)
#boxplot(skew ~ cov*abun, data = s2, notch = T)





s1 <- s[s$sp == 1 & s$sp2 == "Lupine" & s$dist == -7, ]
boxplot(mean ~ cov*abun, data = s1, notch = T, col = 'goldenrod', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov*abun, data = s1, notch = T)
#boxplot(cv ~ cov*abun, data = s1, notch = T)
#boxplot(skew ~ cov*abun, data = s1, notch = T)

s2 <- s[s$sp == 2 & s$sp2 == "Lupine" & s$dist == -7, ]
boxplot(mean ~ cov*abun, data = s2, notch = T, col = 'darkmagenta', ylab = 'mean pollination', xlab = "", xaxt = 'n')
#boxplot(var ~ cov*abun, data = s2, notch = T)
#boxplot(cv ~ cov*abun, data = s2, notch = T)
#boxplot(skew ~ cov*abun, data = s2, notch = T)

