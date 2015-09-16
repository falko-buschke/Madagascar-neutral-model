###############################################################################
#####                                                                     #####
#####  Please consult README.txt file before continuing with this script  #####
#####                                                                     #####
###############################################################################

# This script outlines the functions used to plot the outputs from the paramaeter sensitivity analysis

# set working directory and install required packages
setwd("C:/ADD YOUR OWN WORKING DIRECTORY/Sensitivity_analyses")

####################################################################################
# Creates a plot for parameter K for all birds

# Read in data
birds_K <- read.table("Birds_K.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot
plot(birds_K[,1],birds_K[,2],xlab="Carrying capacity (K)", ylab="Error", ylim=c(0,12), type="o",col=1)
points(birds_K[,1],birds_K[,3],type="o",col=2)
points(birds_K[,1],birds_K[,4],type="o",col=3)
points(birds_K[,1],birds_K[,5],type="o",col=4)
points(birds_K[,1],birds_K[,6],type="o",col=5)
legend("topright", legend = c("alpha","latitude","occupancy", "beta","Total"), col=1:5, pch=16)
abline(v=700)

####################################################################################
# Creates a plot for parameter a for all birds

# Read in data
birds_a <- read.table("Birds_a.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot
plot(birds_a[,1],birds_a[,2],xlab="Short-distance dispersal (a)", ylab="Error", ylim=c(0,2.5), type="o",col=1)
points(birds_a[,1],birds_a[,3],type="o",col=2)
points(birds_a[,1],birds_a[,4],type="o",col=3)
points(birds_a[,1],birds_a[,5],type="o",col=4)
points(birds_a[,1],birds_a[,6],type="o",col=5)
legend("topright", legend = c("alpha","latitude","occupancy", "beta","Total"), col=1:5, pch=16)
abline(v=0.19)

####################################################################################
# Creates a plot for parameter b for all birds

# Read in data
birds_b <- read.table("Birds_b.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot
plot(birds_b[,1],birds_b[,2],xlab="Long-distance dispersal (b)", ylab="Error", ylim=c(0,6), type="o",col=1)
points(birds_b[,1],birds_b[,3],type="o",col=2)
points(birds_b[,1],birds_b[,4],type="o",col=3)
points(birds_b[,1],birds_b[,5],type="o",col=4)
points(birds_b[,1],birds_b[,6],type="o",col=5)
legend("topright", legend = c("alpha","latitude","occupancy", "beta","Total"), col=1:5, pch=16)
abline(v=0.18)

################################################################################
# Creates a plot for parameter K for endemic birds

# Read in data
end_K <- read.table("endemics_K.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot
plot(end_K[,1],end_K[,2],xlab="Carrying capacity (K)", ylab="Error", ylim=c(0,12), type="o",col=1)
points(end_K[,1],end_K[,3],type="o",col=2)
points(end_K[,1],end_K[,4],type="o",col=3)
points(end_K[,1],end_K[,5],type="o",col=4)
points(end_K[,1],end_K[,6],type="o",col=5)
legend("topright", legend = c("alpha","latitude","occupancy", "beta","Total"), col=1:5, pch=16)
abline(v=350)

####################################################################################
# Creates a plot for parameter a for endemic birds

# Read in data
end_a <- read.table("Endemics_a.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot
plot(end_a[,1],end_a[,2],xlab="Short-distance dispersal (a)", ylab="Error", ylim=c(0,5), type="o",col=1)
points(end_a[,1],end_a[,3],type="o",col=2)
points(end_a[,1],end_a[,4],type="o",col=3)
points(end_a[,1],end_a[,5],type="o",col=4)
points(end_a[,1],end_a[,6],type="o",col=5)
legend("topright", legend = c("alpha","latitude","occupancy", "beta","Total"), col=1:5, pch=16)
abline(v=0.095)

####################################################################################
# Creates a plot for parameter b for endemic birds

# Read in data
end_b <- read.table("Endemics_b.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot
plot(end_b[,1],end_b[,2],xlab="Long-distance dispersal (a)", ylab="Error", ylim=c(0,7), type="o",col=1)
points(end_b[,1],end_b[,3],type="o",col=2)
points(end_b[,1],end_b[,4],type="o",col=3)
points(end_b[,1],end_b[,5],type="o",col=4)
points(end_b[,1],end_b[,6],type="o",col=5)
legend("topright", legend = c("alpha","latitude","occupancy", "beta","Total"), col=1:5, pch=16)
abline(v=0.09)
