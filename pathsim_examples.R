load('/Users/JacobSocolar/Dropbox/Work/Bee_Path/bombus.new.352.R')

##### Example 1: Simulating 100 paths on a field of Astragalus, Castilleja, Lupinus, and Penstamon #####
pr <- 20 # patch radius
desired_values <- c(100, 100, 100, 100) # Desired numbers of plants on 50x50 meter field. There will be fewer!
nudiag <- c(3,1,3,.5) # Smoothness parameters for plant field
rho <- matrix(nc=4, c(1, .9, 0.2, .3,
                      .9, 1, 0.2, .5,
                      0.2, 0.2, 1, .5,
                      .3, .5, .5, 1)) # Covariance matrix for plant field that generates plant locations.
#To get something spatially not clumpy, can turn variance down towards zero, or smoothing down towards 0. If variance is
#near zero, then every "site" in the field has the same value. Increasing variance makes sites in the field more varied.
#Then manipulating covariance changes how likely specific plant types are to show up in the same sites together.
#If variance and covariance are both high, then multiple plants of multiple types will be trying to occupy the same site.
#Right now we have a limit of 1 plant per "site" (10 cm x 10 cm). This the EXPECTED covariance matrix, not necessarily the
#realized covariance matrix. Distribution of plants will deviate from distribution specified by field if >1 plant is
#supposed to be in a site.
#NU CANNOT = 0
pt_probs <- c(.9,.8,.7,.3,.15,.12,.1,.07,.03,.02) # Pollen transfer probs. This can be of
# arbitrary length, and will be cropped or padded with trailing zeros as necessary, based on
# the path length.
mcs <- bombus.new.352$coefficients
msd_dist <- c(17.98465, 11.59513)
plantnames <- c('Astrag', 'Castilleja', 'Penst', 'Lupinus')
# (If you want to use the species covariates (and species-by-distance covariates) from the 
# Bombus model, then plantnames needs to include at least two of these species. Any plantnames
# not in this list will be interpreted as zero-covariates, which makes them equivalent to 
# Astragalus for purposes of species ID and species-by-distance, but not for purposes of constancy
# or proportion. Note that you can assign plantnames = NULL if you don't want to bother specifying
# the names--that will give you names like plant1, plant2, etc, up to the length of nudiag.)
patch <- patchsim(nudiag = nudiag, rho = rho, desired_values = desired_values, radius = pr, plantnames = plantnames)
snow.start.time = proc.time()
trans_probs <- transprobs(fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs, msd_dist = msd_dist)
snow.end.time = proc.time()
snow.dtime = snow.end.time - snow.start.time

paths <- pathsim(nsim = 100, nstep = 10, fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs, pollenvector = pt_probs, start = list(patch$pr,"none"), msd_dist = msd_dist, trans_probs = trans_probs)
View(paths)
View(paths$paths) #numbers in cells correspond to row numbers in object fp
#paths$paths shows 100 x 10-step paths
#paths$pollen shows how many times a visit to each individual plant successfully delivered pollen



##### Example 2: Simulating 100 paths on a single-species field #####
pr <- 20 # patch radius
desired_values <- c(100) # Desired numbers of plants on 50x50 meter field
nudiag <- 1 # Smoothness parameters for plant field
rho <- as.matrix(1) # Covariance matrix for plant field. Even though it's a single-element,
# IT MUST STILL BE PASSED AS A MATRIX.
pt_probs <- c(.9,.8,.7,.3,.15,.12,.1,.07,.03,.02) # Pollen transfer probs. This can be of
# arbitrary length, and will be cropped or padded with trailing zeros as necessary, based on
# the path length.
mcs <- bombus.576$coefficients
plantnames <- NULL

patch <- patchsim(nudiag = nudiag, rho = rho, desired_values = desired_values, radius = pr, plantnames = plantnames)
# To re-use the code from the multi-plant system, patchsim does some internal shennanigans here.
# The visible result is that you will no longer get sequentially numbered uniqueIDs or even
# sequential row.names(fp). 
paths <- pathsim(nsim = 100, nstep = 10, fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs, pollenvector = pt_probs, start = list(patch$pr,"none"))
View(paths$paths)
View(paths$pollen)
View(paths$trans_probs[1:10,1:10])


##### Example 3: Resimulating the patch under a constant set of patch statistics, and assessing pollen limitaiton #####
pr <- 20 # patch radius
desired_values <- c(100, 100, 100, 100) # Desired numbers of plants on 50x50 meter field
nudiag <- c(3,1,3,.5) # Smoothness parameters for plant field
rho <- matrix(nc=4, c(1, .9, 0.2, .3,
                      .9, 1, 0.2, .5,
                      0.2, 0.2, 1, .5,
                      .3, .5, .5, 1)) # Covariance matrix for plant field
pt_probs <- c(.9,.8,.7,.3,.15,.12,.1,.07,.03,.02) # Pollen transfer probs. This can be of
# arbitrary length, and will be cropped or padded with trailing zeros as necessary, based on
# the path length.
mcs <- bombus.576$coefficients
plantnames <- c('Astrag', 'Castilleja', 'Penst', 'Lupinus')

pollen_freq <- list()
path_freq <- list()
for(i in 1:10){ # eventually this can be made much bigger
  print(i)
  patch <- patchsim(nudiag = nudiag, rho = rho, desired_values = desired_values, radius = pr, plantnames = plantnames)
  paths <- pathsim(nsim = 100, nstep = 10, fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs, pollenvector = pt_probs, start = list(patch$pr,"none"))
  pollen_freq[[i]] <- paths$pollen
  path_freq[[i]] <- paths$paths
}
pollen_all <- dplyr::bind_rows(pollen_freq, .id = "column_label")

vpi <- data.frame(Species = plantnames, visits_per_ind = 0) #how many visits each species gets per individual.
for(i in 1:nrow(vpi)){
  vpi$visits_per_ind[i] <- sum(pollen_all$success[pollen_all$Species == vpi$Species[i]])/sum(pollen_all$Species == vpi$Species[i])
}
vpi


##### Example 4: plotting results #####
pr <- 20
desired_values <- c(100, 100, 100, 100)
nudiag <- c(3,1,3,.5)
rho <- matrix(nc=4, c(1, .9, 0.2, .3,
                      .9, 1, 0.2, .5,
                      0.2, 0.2, 1, .5,
                      .3, .5, .5, 1))
pt_probs <- c(.9,.8,.7,.3,.15,.12,.1,.07,.03,.02)
mcs <- bombus.576$coefficients
plantnames <- c('Astrag', 'Castilleja', 'Penst', 'Lupinus')
patch <- patchsim(nudiag = nudiag, rho = rho, desired_values = desired_values, radius = pr, plantnames = plantnames)
paths <- pathsim(nsim = 500, nstep = 10, fp = patch$fp, fp_dist = patch$fp_dist, mcs = mcs, pollenvector = pt_probs, start = list(patch$pr,"none"))

fp <- patch$fp
fp$visitcount <- 0
fp$pollencount <- paths$pollen$success
fp$color <- NA
colorlist <- c('firebrick', 'cornflowerblue', 'orange1', 'magenta')
for(i in 1:nrow(fp)){
  fp$visitcount[i] <- length(which(paths$paths == i))
  fp$color[i] <- colorlist[which(plantnames == fp$Species[i])]
}
head(fp)

plot(fp$Y ~ fp$X, pch = 16, col = fp$color, main = 'plant locations')
plot(fp$Y ~ fp$X, pch = 16, col = fp$color, cex = sqrt(fp$visitcount)/5, main = 'number of visits (sqrt)')
# scaled by the square root of the total visit count to visualize better
plot(fp$Y ~ fp$X, pch = 16, col = fp$color, cex = sqrt(fp$pollencount)/5, main = 'number of pollen deliveries (sqrt)')
# scaled by 
