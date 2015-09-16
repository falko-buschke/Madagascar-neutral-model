###############################################################################
#####                                                                     #####
#####  Please consult README.txt file before continuing with this script  #####
#####                                                                     #####
###############################################################################

# This script outlines the functions used simulate bird diversity patterns in Madagascar

# set working directory and install required packages
setwd("C:/ADD YOUR OWN WORKING DIRECTORY")

install.packages("compiler")    # install the 'compiler' package to compile the R code to speed up calculations
library(compiler)

install.packages("vegan")       # install the vegan package for calculating the dipersion field heterogeneity
library(vegan)

install.packages("vegetarian")  # install the vegan package for calculating the dipersion field heterogeneity
library(vegetarian)


#--------------------------------------------------------------------------------
# 1) Open all the required datasets
#--------------------------------------------------------------------------------
 
# Open bird community dataset
birds <-  read.table("Mada_birds.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
birds <- birds[,-1]         # Remove column with species names
richness <- colSums(birds)  # Determine the empirical species richness

# Open grid information dataset
mad.grid <-   read.table("Mada_grid.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coord <- cbind(mad.grid$Long,mad.grid$Lat) # creates a data frame of geographical coordinates of quadrat centroids
energy <- mad.grid$NDVI                    # creates a vector of primary productivity as estimaeted from the NDVI


#--------------------------------------------------------------------------------
# 2) Set up the model requirements
#--------------------------------------------------------------------------------

##########################################################################

# Insert model parameters

para.a <- 0.095    # Parameter 'a' for adjusting short-distance dispersal
para.b <- 0.09     # Parameter 'b' for adjusting long-distance dispersal
K.par <- 350       # Parameter 'K' is the average habitat capacity

##########################################################################

# Create vector of habitat capacity that incorporates gradients of NPP

# Use a conversion factor to determine the habitat capacity of each quadrat relative to the average habitat capacity, K.
conversion <- length(energy)*(energy/sum(energy))   

K <- ceiling(K.par * conversion)     # Number of 'bird units' per quadrat (replace 'conversion' with 1 to remove gradients of NPP)

##########################################################################

# Create a matric of immigation probabilities with source-sink dynamics

# Create a distance matrix (measures distance in kilometers)
distance.1 <- as.matrix(dist(coord, method = "euclidean", diag = TRUE, upper = TRUE))
distance <- distance.1/50  # Convert geographic distance to "grid units"

# Create a dispersal kernel as a combination of exponential and Cauchy functions
dispersal <- para.a^(distance) + (para.b^2/(distance^2 + para.b^2))
diag(dispersal) <- 1      # Make sure to replace all the diagonal values with 1

ccapacity <- K            # Just create a 'holder' for habitat capacity for when you remove directional source-sink dispersal

# Use three steps for the probability of dispersal.
# First calculate the numerator, the the denominator and fianlly dived the former by the latter
numer <- t( apply(dispersal,1,function(dispersal) dispersal*ccapacity))
denom <- colSums(imm.null)
imm.prob <- t( apply(numer,1,function(numer) numer/denom))

##########################################################################

# Just read in all the information needed for the simulation

S <- dim(birds)[1]    # Number of species in assemblage
M <- dim(birds)[2]    # Number of sites/quadrats

species <-  (seq(1:S))    # Just a vector of species identitied
samples <- seq(1:M)       # A vector of quadrat identities

##########################################################################

#--------------------------------------------------------------------------------
# 3) All the functions needed for the nuetral model simulation
#--------------------------------------------------------------------------------

##############################################
#                                             #
#     Step 1(Set up starting conditions)     #
#                                             #
##############################################

# Create a blank species by site matrix for the simulated assemblage
Patches <- matrix(0, nrow=S,ncol=M)

# Randomly select the starting quadrates for each species and polulate these presenece in the matrix
start.cells <- sample(1:M,S,replace=T) 

for (k in 1:S) {
  Patches[k,start.cells[k]] <- 1
}

###################################
#                                 #
#     Step 2(death + replace)     #
#                                 #
###################################

# Creates the 'holder vectors for each quadrat:
ind.k <- rep (NA,M)         # The identitity of the species to die/be replaced
imm.cell <- rep (NA,M)      # The identitity of source quadrate for the colonist
ind.species <- rep (NA,M)   # The identitity of the species that colonises


# Create a function for the death and immigration function
death.imm.fun <- function (comm.t) {
  # Start with the selection of source cell and colinising species
  # But only do this for source quadrats that already contain bird units
  pres.ab <- ifelse(colSums(comm.t)== 0,0,1)    # Just a vector to ensure that the quadrate is occupied
  
  for (l in 1:M) {
    imm.cell[l] <- sample(samples,1,replace=TRUE, prob = (as.vector(imm.prob[l,])*pres.ab))  
  }
  for (k in 1:M) {
    ind.species[k] <- sample (species,1,replace=TRUE, prob = comm.t[,imm.cell[k]])
  }

  # Now choose an individual in each quadrate to be replaced, and then replace it.
  # However, only implement this 'death' process when the number of bird units reaches the habitat capacity 

  abun <- colSums(comm.t)  # A vector of the total abundance in each quadrate

  # If total abundance equals habitat capacity, implement death function, otherwise just colonisation
  for (j in 1:M) {
    if (abun[j] == K[j]) {
      ind.k[j] <- sample (species,1,replace=TRUE, prob=comm.t[,j])    # Choose species to kill
      comm.t[ind.k[j],j] <- comm.t[ind.k[j],j]-1                      # Remove one bird unit
      comm.t[ind.species[j],j] <- comm.t[ind.species[j],j]+1          # Add one bird unit of colonising species
    } else {
      comm.t[ind.species[j],j] <- comm.t[ind.species[j],j] + 1        # Just add colonising bird unit
    }
  }
comm.t
}  

# This uses the 'cmpfun' function fromthe 'compiler' package.
# This is a byte code compiler for R, which creates a compiled body expression that runs much faster.
death.imm <- cmpfun(death.imm.fun)



#--------------------------------------------------------------------------------
# 4) Actually run the whole model for 50 000 timesteps
#--------------------------------------------------------------------------------

ptm <- proc.time()  # Just records the time needed for the simulation
                        
for (i in 1:50000){
  Patches <- death.imm(Patches)   # Implements the function
  Patches
  if (i/1000 == floor(i/1000))    # This is just to reprot the progress of the simulation
  {print(i)                       # Reports every 1000 timesteps
  flush.console()}}
 
proc.time() - ptm    # Reports the time required for the simulation

################################################################
################################################################
############################


#--------------------------------------------------------------------------------
# 5) Make all the plots comparing the empirical and simulated data
#--------------------------------------------------------------------------------

par(mfrow=c(2,2))
#######################################################################################

# A-- Richness distribution
rich.null <- colSums(decostand(Patches,"pa"))     # Simulated richness
rich.obs <- colSums(birds)                        # Emperical richenss

bins <- seq(0,240,by=10)                          # Create bins for the histogram
# determine frequencies and plot the  histograms for empirical data (blue)
his <- hist((rich.obs),breaks=bins, plot=FALSE)
plot(his,main="",col="lightblue",ylim=c(50,200),xlab="Richness")

# Determine the frequencies of simualted data
sim <- hist(rich.null,plot=FALSE, breaks =bins)
# Plot simulated values as red points
points(sim$mids,sim$counts,pch=16,col="red",type="b")


#Calculate the error between observed and simulated frequencies for richness freqencies
E.rich <- (sum ((his$counts-sim$counts)^2)) / (sum((his$counts-(mean(his$counts)))^2))

######################################################################################

# Identify the geographical coordinates of 50km latitudinal bands
y.lat <- mad.grid$Degree
unique.lat <- unique(y.lat)

# B -- Latitude
obs.lat <- rep(NA,length(unique.lat))      # Blank holder for empirical richness in latitudinal bands
sim.lat <- rep(NA,length(unique.lat))      # Blank holder for simulated richness in latitudinal bands

# Use the 'd'  function in the 'vegetarian package to cumulate the number of species in each latitudinal band
for (i in 1:length(unique.lat)) {
  obs.lat[i] <- d(t(birds[,which(y.lat==unique.lat[i])]), lev = "gamma", q = 0)
  sim.lat[i] <- d(t(Patches[,which(y.lat==unique.lat[i])]), lev = "gamma", q = 0)
}

# Make the plot for latitudinal diversity
plot(unique.lat,obs.lat, ylab="Richness", xlab="Latitude",ylim=c(60,120),col="blue",pch=1)
points(unique.lat,sim.lat, col="red",pch=1)

#Calculate the error between observed and simulated frequencies for latitudinal diversity
E.latitude <- (sum ((obs.lat-sim.lat)^2)) / (sum((obs.lat-(mean(obs.lat)))^2))

#######################################################################################

# C-- rank occupancy curves

# Order the observed and simulated species from largest to smallest ranges 
Rank.obs <- sort(rowSums(decostand(birds,"pa")),decreasing=T)
Rank.null <- sort(rowSums(decostand(Patches,"pa")),decreasing=T)

# Make the plot
plot(1:S,Rank.obs,log="y",xlab="Rank",ylab="Occupancy",pch=16,col="blue",xlim=c(120,240))
points(1:S,Rank.null,pch=16,col="red")

#Calculate the error between observed and simulated frequencies for range size (occupancy)
E.occup <- (sum ((Rank.obs-Rank.null)^2)) / (sum((Rank.obs-(mean(Rank.obs)))^2))

#####################################################################################

# D-- beta dissimilarity

# Use the 'vegdist' function in the 'vegan' package to calculate jaccard dissimilarity
jac.obs <-   vegdist(t(birds), method="jaccard", binary=FALSE, diag=F, upper=F)
Neutral.mat <- decostand(t(Patches),"pa")
jac.null <- vegdist(Neutral.mat, method="jaccard", binary=FALSE, diag=F, upper=F)

# Divide similarities into distance classes of 25 km
obs.beta <- as.matrix (jac.obs)
null.beta <- as.matrix (jac.null)
brks <-  seq(0,1600, by=25)
grps.dist <-  findInterval(as.vector(distance.1), brks)

# Turn the matrices into vectors
beta.obs <- as.vector(obs.beta)
beta.null <- as.vector(null.beta)
dist.grps <- as.vector(grps.dist)

# Creat blank holders for mean and standard deviation
obs.mean <- rep(NA,length(brks))
null.mean <- rep(NA,length(brks))
obs.sd <- rep(NA,length(brks))
null.sd <- rep(NA,length(brks))

# Calculate the mean and standard deviation of similarity in 25 km distance classes
for (i in 1:length(brks)) {
  obs.mean[i] <- mean(beta.obs[which(dist.grps==i)])
  null.mean[i] <- mean(beta.null[which(dist.grps==i)])
  obs.sd[i] <- sd(beta.obs[which(dist.grps==i)])
  null.sd[i] <- sd(beta.null[which(dist.grps==i)])
}

# Create the upper and lower bounds of the standards deviations around the mean
obs.upper <- obs.mean + obs.sd
obs.lower <- obs.mean - obs.sd
null.upper <- null.mean + null.sd
null.lower <- null.mean - null.sd

# Make the plot for observed bird assemblages (blue)
plot(1-obs.mean~brks, pch=1,col="blue",ylab="Jaccard similarity", xlab="Geographic distance class", ylim=c(0,1))
arrows (brks,1-obs.mean,brks,1-obs.upper,code=2,len=0.05,angle=90,col="lightblue")
arrows (brks,1-obs.mean,brks,1-obs.lower,code=2,len=0.05,angle=90,col="lightblue")

# Add points for the simulated assemblage
points(1-null.mean~brks,pch=1,col="red")
arrows (brks,1-null.mean,brks,1-null.upper,code=2,len=0.05,angle=90,col="pink")
arrows (brks,1-null.mean,brks,1-null.lower,code=2,len=0.05,angle=90,col="pink")


#Calculate the error between observed and simulated Jaccard similarities
obs.jac <- na.omit(obs.mean)
null.jac <- na.omit(null.mean)
E.beta <- (sum ((obs.jac-null.jac)^2)) / (sum((obs.jac-(mean(obs.jac)))^2))
######################################################################################


# Report all the errors
E.rich
E.latitude
E.occup
E.beta
(E.tot <- E.rich + E.occup + E.beta + E.latitude)


# Report all the R2 values from the unity line regression
1 - E.rich
1 - E.latitude
1 - occup
1 - E.beta


#####################################################################################
