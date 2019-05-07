####Remember, this is just the function file! Run this, but then run actual code in a different file where the parameters
#the go into this function are set for different simulation runs (i.e. don't be an idiot).
library(RandomFields)
inv.logit <- function(x){return(exp(x)/(1+exp(x)))}
patchsim <- function(nudiag, rho, plantnames = NULL, desired_values, radius){
  # nudiag: smoothness parameters for plant field
  # rho: covariance matrix for plants
  # names: character vector of plant species names
  # desired_values: numeric vector of desired *average* plant densities on a 50x50 meter scale
  # radius: radius of patch circle (in meters). Maximum allowable value for now is 25.
  if(!is.vector(nudiag)){stop('nudiag must be a vector')}
  if(!is.matrix(rho)){stop('rho must be a matrix')}
  if(length(nudiag) != nrow(rho)){stop('length nudiag and nrow(rho) differ')}
  if(nrow(rho) != ncol(rho)){stop('rho must be a valid covariance matrix')}
  if(any(eigen(rho)$values < 0)){stop('rho must be a valid covariance matrix')}
  
  nsp <- length(nudiag)
  if(is.null(plantnames)){plantnames <- paste0('plant', 1:nsp)}
  onesp <- F
  if(nsp == 1){
    onesp <- T
    nudiag <- c(nudiag, 1)
    rho <- matrix(c(rho[1,1],0,0,1), nrow = 2)
    plantnames[2] <- 'NULLPLANT'
    desired_values <- c(desired_values, 100)
    nsp <- 2
  }
  
  model <- RMparswmX(nudiag=nudiag, rho=rho) # Random field model
  x.seq <- y.seq <- seq(-10, 10, 0.04)
  z <- RFsimulate(model = model, x=x.seq, y=y.seq) # Simulated realization of field
  pixel_probs_unscaled <- data.frame(inv.logit(z@data[ ,1]))
  for(i in 2:nsp){
    pixel_probs_unscaled <- cbind(pixel_probs_unscaled, inv.logit(z@data[ ,i]))
  }
  pixel_probs <- pixel_probs_unscaled
  for(i in 1:nsp){
    pixel_probs[,i] <- pixel_probs_unscaled[,i] * desired_values[i]/sum(pixel_probs_unscaled[,i])
  }
  npixel <- nrow(pixel_probs)
  floral_locations <- data.frame(rbinom(npixel, 1, pixel_probs[,1]))
  for(i in 2:nsp){
    floral_locations <- cbind(floral_locations, rbinom(npixel, 1, pixel_probs[,i]))
  }
  names(floral_locations) <- plantnames
  nplants <- rowSums(floral_locations)
  problem_pixels <- which(nplants > 1)
  if(length(problem_pixels) > 0){
    for(i in 1:length(problem_pixels)){
      nsp <- nplants[problem_pixels[i]]
      keeper <- sample(c(1:4), 1, prob = floral_locations[problem_pixels[i], ])
      for(j in 1:4){
        floral_locations[problem_pixels[i], j] <- as.numeric(keeper == j)
      }
    }
  }
  plant_rows <- which(nplants > 0)
  floral_locations2 <- data.frame(Species = rep(NA, length(plant_rows)),
                                  X = 0, Y = 0)
  for(i in 1:length(plant_rows)){
    frow <- plant_rows[i]
    floral_locations2$Species[i] <- colnames(floral_locations)[which(floral_locations[frow, ] == 1)]
    floral_locations2$X[i] <- ((frow - 1) %/% 501) / 10
    floral_locations2$Y[i] <- ((frow - 1) %% 501) / 10
  }
  fp <- floral_locations2[((floral_locations2$X - 25)^2 + (floral_locations2$Y - 25)^2) < pr^2, ]
  fp_pattern <- spatstat::ppp(fp$X, fp$Y, c(0,50), c(0,50))
  fp_lambdas <- spatstat::density.ppp(fp_pattern, at="points", edge = F)
  fp$lambda <- fp_lambdas
  
  nplant <- nrow(fp)
  fp$uniqueID <- paste0("indiv", 1:nplant)
  
  fp_dist <- matrix(NA, nrow = nplant, ncol = nplant)
  dimnames(fp_dist) <- list(fp$uniqueID, fp$uniqueID)
  for(i in 1:nplant){
    for(j in 1:i){
      fp_dist[i,j] <- fp_dist[j,i] <- sqrt((fp$X[i]-fp$X[j])^2+
                                             (fp$Y[i]-fp$Y[j])^2)
    }
  }
  diag(fp_dist) <- NA
  
  props <- rep(0, nsp)
  for(i in 1:nsp){
    props[i] <- sum(fp$Species == plantnames[i])/nplant
  }
  
  fp$prop <- 0
  for(i in 1:nplant){
    fp$prop[i] <- props[which(plantnames == fp$Species[i])]
  }
  if(onesp){
    fp_dist <- fp_dist[fp$Species != 'NULLPLANT', fp$Species!='NULLPLANT']
    fp <- fp[fp$Species != 'NULLPLANT', ]
  }
  return(list(fp = fp, fp_dist=fp_dist, pr = pr))
}
