library(spatstat)
library(spatstat.utils)
library(sp)
library(gsl)

###### Functions to initialize bee path at random point on perimeter of convex hull surrounding patch #####

# From https://stat.ethz.ch/pipermail/r-sig-geo/2009-May/005781.html
owin2Polygons <- function(x, id="1") {
  stopifnot(is.owin(x))
  x <- as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

# From https://stackoverflow.com/questions/24496143/sample-points-on-polygon-edge-in-r
rPointOnPerimeter <- function(n, poly) {
  xy <- poly@coords
  dxy <- diff(xy)
  h <- gsl::hypot(dxy[,1], dxy[,2])
  e <- sample(nrow(dxy), n,replace=TRUE, prob=h)
  u <- runif(n)
  p <- xy[e,] + u * dxy[e,]
  return(p)
}

# Function to return a random point and species on the edge of a convex hull around the points
start_fun <- function(site_data, species = "none"){
  convex_hull <- spatstat::convexhull.xy(spatstat::as.ppp(cbind(site_data$X, site_data$Y), c(min(site_data$X),max(site_data$X), min(site_data$Y),max(site_data$Y))))
  chp <- owin2Polygons(convex_hull)
  startpt <- rPointOnPerimeter(1, chp@Polygons[[1]])
  return(list(X = startpt[1], Y = startpt[2], species = species))
}

###### Functions to run the simulation #####
source('/Users/jacobsocolar/Dropbox/Work/Code/beepaths/pathsim.R')  
# Define function transprobs() which computes transition probabilities between plants in a patch
# Define function runsim() which simulates paths on a patch given the transition probabilities

# Define wrapper to transprobs() and runsim() that handles the simulations
run_sims <- function(pd, mcs, colorlist = c('firebrick', 'cornflowerblue', 'magenta', 'goldenrod', 'olivedrab'), 
                     plantnames = c('Astrag', "Cast", "Lupine", "Unknown", "Penst")){
  # pd: one of Bethanne's data files providing locations, identities, and lambda (neighborhood plant density) for a given patch
  # mcs: mean coefficient estimates from RSF
  
  pd$X <- scale(pd$X, scale = F)
  pd$Y <- scale(pd$Y, scale = F)
  pd$uniqueID <- pd$Unique_ID
  
  pd$prop <- NA
  for(i in 1:nrow(pd)){
    pd$prop[i] <- sum(pd$Species == pd$Species[i])/nrow(pd)
  }
  
  pd_dist <- matrix(NA, nrow = nrow(pd), ncol = nrow(pd))
  dimnames(pd_dist) <- list(pd$uniqueID, pd$uniqueID)
  for(i in 1:nrow(pd)){
    for(j in 1:i){
      pd_dist[i,j] <- pd_dist[j,i] <- sqrt((pd$X[i]-pd$X[j])^2 + (pd$Y[i]-pd$Y[j])^2)
    }
  }
  
  trans_probs <- transprobs(fp = pd, fp_dist = pd_dist, mcs = mcs, msd_dist = msd_dist)
  paths <- pathsim(nsim = 1000, nstep = 10, fp = pd, fp_dist = pd_dist, mcs = mcs, pollenvector = pt_amts, start = start_fun(pd), msd_dist = msd_dist, trans_probs = trans_probs, pt_type = "amt")
  
  fp <- pd
  fp$visitcount <- 0
  fp$pollencount <- paths$pollen$success
  fp$color <- NA
  
  for(i in 1:nrow(fp)){
    fp$visitcount[i] <- length(which(paths$paths == i))
    fp$color[i] <- colorlist[which(plantnames == fp$Species[i])]
  }
  return(list(paths = paths, fp = fp, trans_probs = trans_probs, pd = pd, mcs = mcs))
}


##### Load data: fitted model coefficients and patch data files #####
load('/Users/jacobsocolar/Dropbox/Work/Bee_Path/new_data/bombus.352.Rdata') # fitted model object named 'bombus.352'
mcs <- bombus.352$coefficients

pt_amts <- .2*.8^(0:19) # pollen delivery as a function of number of intervening visits--taken from CITATION
msd_dist <- c(17.98465, 11.59513) # mean and sd of pairwise distances used to scale the distance for fitting bombus.352
                                  #  calculated from "ALL paths 2015 2016 w constancy lambda.csv"

# Load patch data files
setwd('/Users/jacobsocolar/Dropbox/Work/Bee_Path')
patch_data_files <- list.files('new_data')[grep('.csv', list.files('new_data'))]
patch_data <- list()
for(i in 1:length(patch_data_files)){
  patch_data[[i]] <- read.csv(paste0('patch_data/', patch_data_files[[i]]))
  # jitter plants with identical coords--none of these plants was ever actually visited.  Jitter necessary to avoid undefined turning angles.
  XY_handle <- paste(patch_data[[i]]$X, patch_data[[i]]$Y, sep = "__")
  dupes <- unique(c(which(duplicated(XY_handle)), which(duplicated(XY_handle, fromLast = T))))
  if(length(dupes > 0)){
    for(j in 1:length(dupes)){
      patch_data[[i]]$X[dupes[j]] <- patch_data[[i]]$X[dupes[j]] + runif(1,-.05, .05) 
      patch_data[[i]]$Y[dupes[j]] <- patch_data[[i]]$Y[dupes[j]] + runif(1,-.05, .05)
    }
  }
}
patch_data[[5]]$Species[patch_data[[5]]$Species == "purpleAstrag"] <- "Unknown" # Bethanne lumped purple Astrag in unknown category when fitting models

mcs_new <- as.data.frame(t(as.matrix(mcs))) # hack to get coefficients into named dataframe

# Zero the coefficients of non-interest
elements_to_zero <- list(c(1,9), 2, 3, 4, 5, 6, c(7,9), 8, 9, 3:6)
for(i in 1:length(elements_to_zero)){
  mcsn <- mcs_new[1, ]
  mcsn[, elements_to_zero[[i]]] <- 0
  mcs_new <- rbind(mcs_new, mcsn)
}

##### Do the simulations #####
rps_output <- list()    # "Real Patch Simulations" output
for(i in 1:length(patch_data_files)){
  pd <- patch_data[[i]]
  rps_output[[i]] <- list()
  for(j in 1:nrow(mcs_new)){
    print(c(i,j))
    mcs_j <- as.numeric(mcs_new[j, ])
    names(mcs_j) <- names(mcs_new)
    rps_output[[i]][[j]] <- run_sims(pd = pd, mcs = mcs_j)
  }
}



