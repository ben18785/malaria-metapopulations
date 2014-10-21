from metapopulationFunctionsAndClasses import classes
from matplotlib import pylab
import numpy as np
import random as random

cTimeTotal = 100
cNumSteps = 100
cNumPatches = 200
cSize = 1000
cInitialPopulation = 100
cMigrantProportion = 0.1
aBreedGamma = 5
bBreedGamma = 6
aFeedGamma = 4
bFeedGamma = 4
cRateParameter = 0.001

aArea = classes.area(cNumPatches,cSize,cInitialPopulation,aBreedGamma,bBreedGamma,aFeedGamma,bFeedGamma)
classes.evolveSystem(aArea,cTimeTotal,cNumSteps,cMigrantProportion,cRateParameter)
# print(random.gammavariate(10,10))



