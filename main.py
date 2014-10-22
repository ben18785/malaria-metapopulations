from metapopulationFunctionsAndClasses import classes
from matplotlib import pylab
import numpy as np
import random as random

cTimeTotal = 3000
cNumSteps = 100
cNumPatches = 50
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
cCovarianceIndicator = 0 # Whether there is covariance among the housing locations
cNumDisturbances = 10
cKernelSigma = 5
vCovarianceParams = [cCovarianceIndicator,cNumDisturbances,cKernelSigma]



# Thresholding of population size < this number set population to zero
cPoplationThresholdIndicator = 1
cPopulationZeroThreshold = 500

aArea = classes.area(cNumPatches,cSize,cInitialPopulation,aBreedGamma,bBreedGamma,aFeedGamma,bFeedGamma,vCovarianceParams)
classes.evolveSystem(aArea,cTimeTotal,cNumSteps,cMigrantProportion,cRateParameter,cInitialTime,cTimeHegIn,cPoplationThresholdIndicator,cPopulationZeroThreshold)
# print(random.gammavariate(10,10))



