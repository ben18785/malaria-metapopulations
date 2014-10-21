from metapopulationFunctionsAndClasses import classes
from matplotlib import pylab
import numpy as np
import random as random

aArea = classes.area(100,1000)

print(aArea.getTimeIncrementAndSelectEvent())
Y=aArea.getPatches()[0].runODE(10,100)
print(aArea.getPatches()[0].runODEAndUpdateGroup(10,100))



