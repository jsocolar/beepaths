load('/Users/JacobSocolar/Dropbox/Work/Bee_Path/sims.Rdata')

str(simulations)
class(simulations)
length(simulations) # one list element for each row of the parameter set

class(simulations[[1]]) # each element is a list
length(simulations[[1]]) # one element for the output of each core

class(simulations[[1]][[1]]) 
length(simulations[[1]][[1]]) # one element for each replicate as controlled by the upper bound on j

class(simulations[[1]][[1]][[1]])
# This is the output from the pathsim function
names(simulations[[1]][[1]][[1]])

simulations[[1]][[1]][[1]]$args

