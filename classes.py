import numpy as np
import random as random
import scipy as sp
from scipy import integrate as integrate
import math as math
from matplotlib import pylab as plt
import time as time
import matplotlib as mpl


class patch:
    def __init__(self,aLocation,anumBreedSites,anumFeedSites,aMosquitoGroup):
        self.location = aLocation
        self.numBreedSites = anumBreedSites
        self.numFeedSites = anumFeedSites
        self.mosquitoGroup = aMosquitoGroup

        cSO = 0.15
        cSH = 0.5
        cSM = 0.3
        cOmega = 40
        cKappa = 120

        cO = anumBreedSites*cSO*cSO*sp.pi*cOmega/cKappa

        cH = anumFeedSites*sp.pi*cSH*cSH


        cNu = 0.12*cO
        cE = 0.6
        cGammaJ = 0.1
        cMuJ = 0.1
        cAlpha = 0.05/anumBreedSites
        cMuM = 0.1
        cMuU = 0.1
        cM = 0.01*sp.pi*cSM*cSM
        cMuH = 0.1
        cGammaH = 0.015*cH
        cMuO = 0.1

        self.vParams = (cKappa,cNu,cE,cGammaJ,cMuJ,cAlpha,cMuM,cMuU,cM,cMuH,cGammaH,cMuO)


    def getNumBreedSites(self):
        return self.numBreedSites


    def getNumFeedSites(self):
        return self.numFeedSites

    def getTotalPopulation(self):
        return self.mosquitoGroup.getTotalPopulation()

    def getLocation(self):
        return self.location.getLocation()

    def getMosquitoGroup(self):
        return self.mosquitoGroup

    def getParams(self):
        return self.vParams

    def getPopulation(self):
        return self.mosquitoGroup.getPopulation()

    def runODE(self,cTime,cNumSteps):
        vTime = np.linspace(0,cTime,cNumSteps)
        yInit = self.mosquitoGroup.getPopulation()
        Y = integrate.odeint(derivativeODESystem,yInit,vTime,args=(self.vParams,))
        return Y

    def runODEAndUpdateGroup(self,cTime,cNumSteps):
        Y = self.runODE(cTime,cNumSteps)
        newGroup = mosquitoGroup(Y[-1,0],Y[-1,1],Y[-1,2],Y[-1,3],Y[-1,4],Y[-1,5],Y[-1,6],Y[-1,7],Y[-1,8],Y[-1,9])
        self.mosquitoGroup.removeMosquitoes(self.mosquitoGroup) # Clears the group
        self.immigrate(newGroup)


    # Once a migration to another patch is confirmed, this function removes the outgoing group numbers from the existing mosquitoes
    def emigrate(self,outgoingMosquitoGroup):
        self.mosquitoGroup.removeMosquitoes(outgoingMosquitoGroup)

    # Once a migration from another patch is confirmed, this function adds the incoming group to the existing mosquitoes
    def immigrate(self,incomingMosquitoGroup):
        self.mosquitoGroup.addMosquitoes(incomingMosquitoGroup)

    def migration(self,outgoingMosquitoGroup,otherPatch):
        self.emigrate(outgoingMosquitoGroup)
        otherPatch.immigrate(outgoingMosquitoGroup)

    # Gets the rate for a particular migration from this patch to another, which is a function of both distance and mosquito number
    def getRate(self,anotherPatch,cRateParameter):
        cdistance = getDistance(self.getLocation(),anotherPatch.getLocation())
        cRate = cRateParameter*self.mosquitoGroup.getTotalAdultPopulation()/(cdistance**2)
        return cRate

    def getAllRates(self,aArea,cRateParameter):
        vRates = []
        for patches in aArea.getPatches():
            cdistance = getDistance(self.getLocation(),patches.getLocation())
            if cdistance > 0:
                vRates.append([self,patches,self.getRate(patches,cRateParameter)])
        return vRates

    def setPopulationToZero(self):
        self.mosquitoGroup.setGroupSizeToZero()





class location:
    def __init__(self,ax,ay):
        self.x = ax
        self.y = ay

    def getLocation(self):
        return np.array([self.x,self.y])

class mosquitoGroup:
    def __init__(self,aJx,aJY0=None,aJY1=None,aMY0=None,aMY1=None,aU=None,aHY0=None,aHY1=None,aOY0=None,aOY1=None):
        if aJY0 is None: # For case of supplying arguments as a vector (therefore one input)
            self.Jx = aJx[0]
            self.JY0 = aJx[1]
            self.JY1 = aJx[2]
            self.MY0 = aJx[3]
            self.MY1 = aJx[4]
            self.U = aJx[5]
            self.HY0 = aJx[6]
            self.HY1 = aJx[7]
            self.OY0 = aJx[8]
            self.OY1 = aJx[9]
        else:
            self.Jx = aJx
            self.JY0 = aJY0
            self.JY1 = aJY1
            self.MY0 = aMY0
            self.MY1 = aMY1
            self.U = aU
            self.HY0 = aHY0
            self.HY1 = aHY1
            self.OY0 = aOY0
            self.OY1 = aOY1


    def getPopulation(self):
        return np.array([self.Jx,self.JY0,self.JY1,self.MY0,self.MY1,self.U,self.HY0,self.HY1,self.OY0,self.OY1])

    def getTotalPopulation(self):
        return np.sum(np.array([self.Jx,self.JY0,self.JY1,self.MY0,self.MY1,self.U,self.HY0,self.HY1,self.OY0,self.OY1]))

    def getTotalAdultPopulation(self):
        return np.sum(np.array([self.MY0,self.MY1,self.U,self.HY0,self.HY1,self.OY0,self.OY1]))


    def addMosquitoes(self,anotherMosquitoGroup):
        self.Jx += anotherMosquitoGroup.Jx
        self.JY0 += anotherMosquitoGroup.JY0
        self.JY1 += anotherMosquitoGroup.JY1
        self.MY0 += anotherMosquitoGroup.MY0
        self.MY1 += anotherMosquitoGroup.MY1
        self.U += anotherMosquitoGroup.U
        self.HY0 += anotherMosquitoGroup.HY0
        self.HY1 += anotherMosquitoGroup.HY1
        self.OY0 += anotherMosquitoGroup.OY0
        self.OY1 += anotherMosquitoGroup.OY1

    def removeMosquitoes(self,anotherMosquitoGroup):
        self.Jx -= anotherMosquitoGroup.Jx
        self.JY0 -= anotherMosquitoGroup.JY0
        self.JY1 -= anotherMosquitoGroup.JY1
        self.MY0 -= anotherMosquitoGroup.MY0
        self.MY1 -= anotherMosquitoGroup.MY1
        self.U -= anotherMosquitoGroup.U
        self.HY0 -= anotherMosquitoGroup.HY0
        self.HY1 -= anotherMosquitoGroup.HY1
        self.OY0 -= anotherMosquitoGroup.OY0
        self.OY1 -= anotherMosquitoGroup.OY1

    def setGroupSizeToZero(self):
        self.Jx = 0
        self.JY0 = 0
        self.JY1 = 0
        self.MY0 = 0
        self.MY1 = 0
        self.U = 0
        self.HY0 = 0
        self.HY1 = 0
        self.OY0 = 0
        self.OY1 = 0


class area:
    def __init__(self,anumPatches,aSize,aInitialPopulation,aBreedGamma,bBreedGamma,aFeedGamma,bFeedGamma,vCovarianceParams):
        self.vPatches = []
        for i in range(0,anumPatches):
            cNumBreedRand = random.gammavariate(aBreedGamma,bBreedGamma)
            cNumFeedRand = random.gammavariate(aFeedGamma,bFeedGamma)
            if vCovarianceParams[0] == 0:
                self.vPatches.append(patch(randomLocation(aSize),cNumBreedRand,cNumFeedRand,mosquitoGroup(aInitialPopulation,aInitialPopulation,0,0,0,0,0,0,0,0)))
            else:
                self.vPatches.append(patch(randomCovarianceLocation(aSize,vCovarianceParams),cNumBreedRand,cNumFeedRand,mosquitoGroup(aInitialPopulation,aInitialPopulation,0,0,0,0,0,0,0,0)))

        self.numPatches = anumPatches
        self.Size = aSize

    def getPatches(self):
        return self.vPatches

    def getNumPatches(self):
        return self.numPatches

    def getPatchesLocation(self):
        vLocations = []
        for patches in self.vPatches:
            vLocations.append(patches.getLocation())

        return np.array(vLocations)

    def getPatchPopulation(self):
        vPatchPopulation = []
        for patches in self.vPatches:
            vPatchPopulation.append(patches.getTotalPopulation())

        return np.array(vPatchPopulation)

    def getTotalPopulation(self):
        cCount = 0
        for patches in self.vPatches:
            cCount += patches.getTotalPopulation()

        return cCount

    def getAllRates(self,cRateParameter):
        mRates = []
        for patches in self.vPatches:
            mRates.append(patches.getAllRates(self,cRateParameter))
        return mRates

    def getTimeIncrement(self):
        cSum = sumAllRates(self.getAllRates())
        cRand = random.random()
        cTime = (1/cSum)*math.log(1/cRand)
        return cTime


    def getTimeIncrementAndSelectEvent(self,cRateParameter):
        mRates = self.getAllRates(cRateParameter)
        vSum = sumPatchRates(mRates)
        cSum = sum(vSum)
        cRand = random.random()
        cTime = (1/cSum)*math.log(1/cRand)

        # First select a patch based on relative rates
        cPatchNumber = pickPatch(vSum,cSum,self.numPatches)

        # Now select the destination amongst the other patches
        vRatesWithinPatch = sumWithinPatch(mRates,cPatchNumber,self.numPatches-1)
        # print(vRatesWithinPatch)
        # print(mRates[cPatchNumber])
        cSumWithinPatch = vSum[cPatchNumber]
        cOtherPatchNumber = pickPatch(vRatesWithinPatch,cSumWithinPatch,self.numPatches-1)
        aOtherPatch = mRates[cPatchNumber][cOtherPatchNumber][1]
        # print(cOtherPatchNumber)
        # print(aOtherPatch)
        # time.sleep(3)
        return [cTime,cPatchNumber,aOtherPatch]

    def getSize(self):
        return self.Size

    def getHegIndicator(self):
        vHegIndicator = np.zeros(self.numPatches)
        k = 0
        for patches in self.vPatches:
            vPopulationTemp = patches.getPopulation()
            if vPopulationTemp[2] > 0:
                vHegIndicator[k] = 1
            k += 1
        return vHegIndicator

    def evolveODE(self,cTime,cNumSteps,cPopulationZeroThreshold):
        for patches in self.vPatches:
            if patches.getTotalPopulation() > cPopulationZeroThreshold:
                patches.runODEAndUpdateGroup(cTime,cNumSteps)
            else:
                patches.setPopulationToZero()

# A construct that is used for creating spatially-correlated landscapes
class Disturbance:
    def __init__(self,aLocation):
        self.location = aLocation
        self.vPatches = []

    def addPatch(self,vCovarianceParams,anumBreedSites,anumFeedSites,aMosquitoGroup):
        aX = self.location[0]
        aY = self.location[1]

        # Generate x and y of patch
        cKernelSigma = vCovarianceParams[2]
        aPatchX = aX + random.normalvariate(0,cKernelSigma)
        aPatchY = aY + random.normalvariate(0,cKernelSigma)
        aLocation = location(aPatchX,aPatchY)

        # Create patch
        aPatch = patch(aLocation,anumBreedSites,anumFeedSites,aMosquitoGroup)
        return aPatch

    def createPatches(self,vCovarianceParams):
        pass


def randomLocation(cSize):
    ax = cSize*random.random()
    ay = cSize*random.random()
    aLocation = location(ax,ay)

    return aLocation

def getDistance(aLocation,bLocation):
    if type(aLocation) is np.ndarray:
        ax1 = aLocation[0]
        ay1 = aLocation[1]

        ax2 = bLocation[0]
        ay2 = bLocation[1]
    else:
        ax1 = aLocation.getLocation()[0]
        ay1 = aLocation.getLocation()[1]

        ax2 = bLocation.getLocation()[0]
        ay2 = bLocation.getLocation()[1]

    return math.sqrt((ax1-ax2)**2 + (ay1-ay2)**2)

def derivativeODESystem(Y,t,vParams):
    Jx = Y[0]; JY0 = Y[1]; JY1 = Y[2]; MY0 = Y[3]; MY1 = Y[4]; U = Y[5]; HY0 = Y[6]
    HY1 = Y[7]; OY0 = Y[8]; OY1 = Y[9]

    cKappa = vParams[0]
    cNu = vParams[1]
    cE = vParams[2]
    cGammaJ = vParams[3]
    cMuJ = vParams[4]
    cAlpha = vParams[5]
    cMuM = vParams[6]
    cMuU = vParams[7]
    cM = vParams[8]
    cMuH = vParams[9]
    cGammaH = vParams[10]
    cMuO = vParams[11]

    dJxdt = cKappa*cNu*(0.5*OY0 + 0.5*(1-cE)*OY1) - cGammaJ*Jx - cMuJ*Jx - cAlpha*Jx*(Jx + JY0 + JY1)
    dJY0dt = cKappa*cNu*0.5*OY0 - cGammaJ*JY0 - cMuJ*JY0 - cAlpha*JY0*(Jx + JY0 + JY1)
    dJY1dt = cKappa*cNu*0.5*(1+cE)*OY1 - cGammaJ*JY1 - cMuJ*JY1 - cAlpha*JY1*(Jx + JY0 + JY1)
    dM0dt = cGammaJ*JY0 - cMuM*MY0
    dM1dt = cGammaJ*JY1 - cMuM*MY1
    dUdt = cGammaJ*Jx - cMuU*U - cM*U*(MY0 + MY1)
    dHY0dt = cM*U*MY0 + cNu*OY0 - cMuH*HY0 - cGammaH*HY0
    dHY1dt = cM*U*MY1 + cNu*OY1 - cMuH*HY1 - cGammaH*HY1
    dOY0dt = -cNu*OY0 - cMuO*OY0 + cGammaH*HY0
    dOY1dt = -cNu*OY1 - cMuO*OY1 + cGammaH*HY1

    return np.array([dJxdt,dJY0dt,dJY1dt,dM0dt,dM1dt,dUdt,dHY0dt,dHY1dt,dOY0dt,dOY1dt])

def sumAllRates(mRates):
    cSum = 0
    for i in range(0,len(mRates)):
        for j in range(0,len(mRates)-1):
            cSum += mRates[i][j][2]
    return cSum

def sumPatchRates(mRates):
    vSum = []
    for i in range(0,len(mRates)):
        cSum = 0
        for j in range(0,len(mRates)-1):
            cSum += mRates[i][j][2]
        vSum.append(cSum)
    return vSum

def sumWithinPatch(mRates,cPatchNumber,cNumOtherPatches):
    vSum = []
    mRates1 = mRates[cPatchNumber]
    for i in range(0,cNumOtherPatches):
        vSum.append(mRates1[i][2])

    return vSum


def pickPatch(vSum,cSum,cNumPatches):
    test = 0
    while test == 0:
        cRandIndex = random.randint(0,cNumPatches-1)

        if random.random()<vSum[cRandIndex]/cSum:
            # print(vSum[cRandIndex]/cSum)
            test = 1
            return cRandIndex

# Calculates rates, then gives time till next migration, then selects event, implements the migration and evolves
# the ODEs for that time period
def evolveSystem(aArea,cTimeTotal,cNumSteps,cMigrantProportion,cRateParameter,cInitialTime,cTimeHegIn,cPoplationThresholdIndicator,cPopulationZeroThreshold):
    cTimeCounter = 0
    cSteps = 0
    vTotalPopulation = []
    vTimeEvents = []
    vTimeEvents.append(0)
    vLocations = aArea.getPatchesLocation()

    # Initial period of running to get individual cases to near equilibrium before starting migration
    if cPoplationThresholdIndicator == 0:
        aArea.evolveODE(cInitialTime,cNumSteps,0)
    else:
        aArea.evolveODE(cInitialTime,cNumSteps,cPopulationZeroThreshold)
    vTotalPopulation.append(aArea.getTotalPopulation())

    HegSwitch = 0

    f, axarr = plt.subplots(2)
    while cTimeCounter < cTimeTotal:

        # Put Heg into a random patch if it is sufficiently big
        if cTimeCounter > cTimeHegIn and HegSwitch == 0:
            while HegSwitch == 0:
                aHegGroup = mosquitoGroup(0,0,1000,0,0,0,0,0,0,0)
                cRandEntry = random.randint(0,aArea.getNumPatches()-1)
                if aArea.getPatches()[cRandEntry].getTotalPopulation() > 0:
                    HegSwitch = 1
                    aArea.getPatches()[cRandEntry].immigrate(aHegGroup)
                    print(aArea.getPatches()[cRandEntry].getLocation())


        [cTimeIncrement,cPatchNumber,aOtherPatch] = aArea.getTimeIncrementAndSelectEvent(cRateParameter)
        if cPoplationThresholdIndicator == 0:
            aArea.evolveODE(cTimeIncrement,cNumSteps,0)
        else:
            aArea.evolveODE(cTimeIncrement,cNumSteps,cPopulationZeroThreshold)
        aGroup = mosquitoGroup(cMigrantProportion*aArea.getPatches()[cPatchNumber].getPopulation())
        aArea.getPatches()[cPatchNumber].migration(aGroup,aOtherPatch)
        cTimeCounter += cTimeIncrement
        cSteps += 1
        vPatchPopulations = aArea.getPatchPopulation()

        axarr[0].scatter(vLocations[:,0],vLocations[:,1],s=0.03*vPatchPopulations,c=aArea.getHegIndicator(),vmin=0, vmax=1)
        axarr[0].set_xlim([0,aArea.getSize()])
        axarr[0].hold(False)
        plt.draw()
        # print(aArea.getHegIndicator())


        vTotalPopulation.append(aArea.getTotalPopulation())
        vTimeEvents.append(cTimeCounter)
        print(cTimeCounter)

        axarr[1].plot(np.array(vTimeEvents),np.array(vTotalPopulation),c='b')
        axarr[1].hold(False)
        axarr[1].set_xlim([0,cTimeTotal])
        axarr[1].set_ylim([0,max(np.array(vTotalPopulation))+0.1e6])

        f.show()


def randomCovarianceLocation(aSize,vCovarianceParams):
    for i in range(0,nNumDisturbances):
        pass
