from metapopulationFunctionsAndClasses import classes
from matplotlib import pylab
import numpy as np
import random as random
from metapopulationFunctionsAndClasses import functions

cTimeTotal = 3000
cNumSteps = 100
cNumPatches = 200
cSize = 1000
cInitialPopulation = 1000
cMigrantProportion = 0.1
aBreedGamma = 5
bBreedGamma = 5
aFeedGamma = 3
bFeedGamma = 3
cRateParameter = 0.01
cInitialTime = 100 # Time to run ODEs before the migration starts
cTimeHegIn = 10 # Time to put the Heg in after the migration starts

# Landscape covariance terms
cCovarianceIndicator = 1 # Whether there is covariance among the housing locations
cNumDisturbances = 2
cKernelSigma = 100
cMeanNumPatchesPoisson = 50
vCovarianceParams = [cCovarianceIndicator,cNumDisturbances,cKernelSigma,cMeanNumPatchesPoisson]



# Thresholding of population size < this number set population to zero
cPoplationThresholdIndicator = 1
cPopulationZeroThreshold = 0

# Generate area and distances
aArea = classes.area(cNumPatches,cSize,cInitialPopulation,aBreedGamma,bBreedGamma,aFeedGamma,bFeedGamma,vCovarianceParams)
functions.getDistances(aArea)

# Evolve system
functions.evolveSystem(aArea,cTimeTotal,cNumSteps,cMigrantProportion,cRateParameter,cInitialTime,cTimeHegIn,cPoplationThresholdIndicator,cPopulationZeroThreshold)


