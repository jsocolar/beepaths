library(spatstat)
library(spatstat.utils)
library(sp)
library(gsl)

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

# Function to run the simulation and produce a data object
run_sims <- function(pd, mcs, colorlist = c('firebrick', 'cornflowerblue', 'magenta', 'goldenrod'), 
                     plantnames = c('Astrag', "Cast", "Lupine", "Unknown")){
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
      pd_dist[i,j] <- pd_dist[j,i] <- sqrt((pd$X[i]-pd$X[j])^2+
                                             (pd$Y[i]-pd$Y[j])^2)
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


load('/Users/JacobSocolar/Dropbox/Work/Bee_Path/bombus.new.352.R')
pt_amts <- .2*.8^(0:19)
mcs <- bombus.new.352$coefficients
msd_dist <- c(17.98465, 11.59513)

source('/Users/jacobsocolar/Dropbox/Work/Code/beepaths/pathsim.R')

setwd('/Users/jacobsocolar/Dropbox/Work/Bee_Path')
patch_data_files <- list.files('patch_data')
patch_data <- list()
for(i in 1:length(patch_data_files)){
  patch_data[[i]] <- read.csv(paste0('patch_data/', patch_data_files[[i]]))
}
# Above, the long and short side calculation was to inform the smoothing standard deviation, which is implemented in Bethanne's code

#i <- 1

mcs_new <- as.data.frame(t(as.matrix(mcs)))
rps_output <- list()    # "Real Patch Simulations" output
for(i in 1:length(patch_data_files)){
  pd <- patch_data[[i]]
  rps_output[[i]] <- list()
  for(j in 1:nrow(mcs_new)){
    mcs_j <- as.numeric(mcs_new[j, ])
    names(mcs_j) <- names(mcs_new)
    rps_output[[i]][[j]] <- run_sims(pd = pd, mcs = mcs_j)
  }
}


# Example plotting code


fp <- test$fp
head(fp)

plot(fp$Y ~ fp$X, pch = 16, col = fp$color, main = 'plant locations', asp = 1, xlab = "", ylab = "")

fpa <- fp[fp$Species == "Astrag", ]
fpb <- fp[fp$Species != "Astrag", ]

plot(fpa$Y ~ fpa$X, pch = 16, col = fpa$color, cex = sqrt(fpa$visitcount)/5, main = 'number of visits (sqrt)', asp = 1, xlab = "", ylab = "")
points(fpb$Y ~ fpb$X, pch = 16, col = fpb$color, cex = sqrt(fpb$visitcount)/5)

# scaled by the square root of the total visit count to visualize better

plot(fpa$Y ~ fpa$X, pch = 16, col = fpa$color, cex = sqrt(fpa$pollencount)/2, main = 'conspecific pollen delivery (sqrt)', asp = 1, xlab = "", ylab = "")
points(fpb$Y ~ fpb$X, pch = 16, col = fpb$color, cex = sqrt(fpb$pollencount)/2)


