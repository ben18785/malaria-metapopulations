import numpy as np
import random as random
import scipy as sp
from scipy import integrate as integrate
import math as math
from matplotlib import pylab as plt
from metapopulationFunctionsAndClasses import classes as classes
import time as time


# Calculates rates, then gives time till next migration, then selects event, implements the migration and evolves
# the ODEs for that time period
def evolveSystem(aArea,cTimeTotal,cNumSteps,cMigrantProportion,cRateParameter,cInitialTime,cTimeHegIn,cPoplationThresholdIndicator,cPopulationZeroThreshold):
    cTimeCounter = 0
    cSteps = 0
    vTotalPopulation = []
    vHegPopulation = []
    vTimeEvents = []
    vTimeEvents.append(0)
    vLocations = aArea.getPatchesLocation()
    # Initial period of running to get individual cases to near equilibrium before starting migration
    if cPoplationThresholdIndicator == 0:
        aArea.evolveODE(cInitialTime,cNumSteps,0)
    else:
        aArea.evolveODE(cInitialTime,cNumSteps,cPopulationZeroThreshold)

    vTotalPopulation.append(aArea.getTotalPopulation())
    vHegPopulation.append(aArea.getHegPopulation())

    HegSwitch = 0
    PlotTimer = 0
    PlotSwitch = 0
    f, axarr = plt.subplots(3)
    while cTimeCounter < cTimeTotal:

        # Put Heg into a random patch if it is sufficiently big
        if cTimeCounter > cTimeHegIn and HegSwitch == 0:
            while HegSwitch == 0:
                aHegGroup = classes.mosquitoGroup(0,0,1000,0,0,0,0,0,0,0)
                cRandEntry = random.randint(0,aArea.getNumPatches()-1)
                if aArea.getPatches()[cRandEntry].getTotalPopulation() > 5000:
                    HegSwitch = 1
                    aArea.getPatches()[cRandEntry].immigrate(aHegGroup)
                    print(aArea.getPatches()[cRandEntry].getLocation())
                    print(aArea.getPatches()[cRandEntry].getNumBreedSites(),aArea.getPatches()[cRandEntry].getNumFeedSites())


        [cTimeIncrement,cPatchNumber,aOtherPatch] = aArea.getTimeIncrementAndSelectEvent(cRateParameter)
        if cPoplationThresholdIndicator == 0:
            aArea.evolveODE(cTimeIncrement,cNumSteps,0)
        else:
            aArea.evolveODE(cTimeIncrement,cNumSteps,cPopulationZeroThreshold)
        aGroup = classes.mosquitoGroup(cMigrantProportion*aArea.getPatches()[cPatchNumber].getPopulation())
        aArea.getPatches()[cPatchNumber].migration(aGroup,aOtherPatch)
        cTimeCounter += cTimeIncrement
        cSteps += 1
        vPatchPopulations = aArea.getPatchPopulation()
        vTotalPopulation.append(aArea.getTotalPopulation())
        vHegPopulation.append(aArea.getHegPopulation())
        vTimeEvents.append(cTimeCounter)
        print(cTimeCounter)
        PlotTimer += cTimeIncrement

        if PlotTimer > 5:
            PlotSwitch = 1

        if PlotSwitch == 1:
            axarr[0].scatter(vLocations[:,0],vLocations[:,1],s=0.03*vPatchPopulations,c=aArea.getHegIndicator(),vmin=0, vmax=1)
            axarr[0].set_xlim([0,aArea.getSize()])
            axarr[0].set_ylim([0,aArea.getSize()])
            axarr[0].hold(False)
            axarr[0].set_title("Village locations")
            plt.draw()

            p1, = axarr[1].plot(np.array(vTimeEvents),np.array(vTotalPopulation),c='b',linewidth=3)
            axarr[1].hold(True)
            p2, = axarr[1].plot(np.array(vTimeEvents),np.array(vHegPopulation),c='r',linewidth=3)
            axarr[1].legend([p1, p2],["Total","Heg"])
            axarr[1].hold(False)
            axarr[1].set_xlim([0,cTimeTotal])
            axarr[1].set_ylim([0,max(np.array(vTotalPopulation))+0.01e6])
            axarr[1].set_ylabel("Population size")
            axarr[1].set_xlabel("Time")
            axarr[1].set_title("Population dynamics")

            # Plot breeding and feeding site densities
            vColour = colourBreedFeed(np.array(aArea.getNumBreedSites()),np.array(aArea.getNumFeedSites()))
            axarr[2].set_title("Breeding and feeding site distributions")
            axarr[2].scatter(np.array(aArea.getNumBreedSites()),np.array(aArea.getNumFeedSites()),c=vColour)
            axarr[2].set_ylabel("Feeding sites")
            axarr[2].set_xlabel("Breeding sites")
            axarr[2].set_xlim([0,128])
            axarr[2].set_ylim([0,max(np.array(aArea.getNumFeedSites()))])

            f.show()
            PlotTimer = 0
            PlotSwitch = 0
            time.sleep(1e-6)

def getDistances(aArea):
    vPatches = aArea.getPatches()
    for patches in vPatches:
        for otherPatches in vPatches:
            cDistance = classes.getDistance(patches.getLocation(),otherPatches.getLocation())
            patches.setDistances(cDistance)

# Creates a colour vector to illustrate where the various points are on Ace's breed-feed ODE diagram. Have approximated
# the regions by y = -(1/8)x + 6 for WT non-survive boundary. For Heg survive boundary have assumed y = 362/x.
def colourBreedFeed(vNumBreed,vNumFeed):
    cLen = len(vNumBreed)
    vColour = np.zeros(cLen)
    for i in range(0,cLen):
        if vNumFeed[i] + (1/8)*vNumBreed[i] - 6 > 0: # WT non-survive test
            if vNumFeed[i] - (362/vNumBreed[i]) > 0: # Heg survive test
                vColour[i] = 2
            else:
                vColour[i] = 1
    return vColour
