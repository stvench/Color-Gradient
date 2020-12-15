from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


import inputFiles
import createStructs


def main(merge, waveband1, waveband2, galNum):
    minMaxRatio, minorAxis, majorAxis, axisRadians, inputCenterR, inputCenterC = inputFiles.read_arcsTSV(galNum)
    ############### 2.
    pixelLoc, rows, cols = inputFiles.read_clusMask(waveband=waveband1,galNum=galNum,group=True)
    ############### 3.
    armsPixels = createStructs.get_largestArm(galaxyArms=pixelLoc)
    ############### 4.
    arcsEllipse_Positions, overallMinTheta, overallMaxTheta = createStructs.arm_to_ArcsEllipse(majorAxis=majorAxis, minMaxRatio=minMaxRatio, axisRadians=axisRadians, armsPixels=armsPixels, center=(inputCenterR, inputCenterR))



    # LOOP 1 (CHANGES ELLIPSEINFO OBJECTS THETAS RELATIVE TO START OF 0 from FRONT) ###############################################################################################################################
    needSub360 = True if overallMaxTheta-overallMinTheta >= 360 else False
    newOverallMinTheta = None
    newOverallMaxTheta = None
    for ae_pObj in arcsEllipse_Positions:
        arcUpdatedThetas = []
        for theta,i,j in ae_pObj.arc:
            if (needSub360):
                if theta>360:
                    theta-=360
            if newOverallMinTheta is None or theta<newOverallMinTheta:
                newOverallMinTheta = theta
            if newOverallMaxTheta is None or theta>newOverallMaxTheta:
                newOverallMaxTheta = theta
            arcUpdatedThetas.append( (theta,i,j) )
        ae_pObj.arc = arcUpdatedThetas
    armsFrontTheta = arcsEllipse_Positions[0].arc[0][0]
    armsEndTheta   = arcsEllipse_Positions[-1].arc[-1][0]
    sWise = True if armsFrontTheta>armsEndTheta else False



    ####################### Figure out and write out reasoning for the +1, +2, -2 (They are definitely needed though) ###############################################################################################################################
    minMajAxLen = arcsEllipse_Positions[0].majorAxisLen
    maxMajAxLen = arcsEllipse_Positions[-1].majorAxisLen
    x = np.arange(0,newOverallMaxTheta-newOverallMinTheta+1)
    y = np.arange(minMajAxLen-2,maxMajAxLen+2)
    FINALPLOT = np.ones((len(y),len(x),3), dtype=float) # array of [1.0,1.0,1.0] (float needed as plt.imshow() for RBGs floats bounded by [0...1], with 1 as WHITE)
    fits1 = inputFiles.readFits(waveband1,galNum)
    fits2 = inputFiles.readFits(waveband2,galNum)
    # LOOP 2            ##################################       ADD A LOOP HERE FOR THE VALUES OF HOW MANY TO MERGE, DON'T NEED TO REDO EVERYTHING BEFORE
    for i in range(merge,len(arcsEllipse_Positions)-merge):
        current_aepObj = arcsEllipse_Positions[i]
        neighbor_aepObjs = [arcsEllipse_Positions[j] for j in np.arange(i-merge,i+merge+1) if (i!=j)]
        uniqueNeighborThetas = createStructs.remvSimThetas(middle=current_aepObj, neighbors=neighbor_aepObjs)


        # LOOP 1 (GETS THE RELATIVE MIN/MAX FOR THE CURRENT RADIUS)   ################# DO THIS DIRECTLY IN remvSIMTHETAS, ONE LESS LOOP
        curRadiusMinFlux = None
        curRadiusMaxFLux = None
        for theta,pixelList in uniqueNeighborThetas:
            for i,j in pixelList:
                avgFlux = createStructs.calcFlux(i,j,fits1,fits2)
                if (curRadiusMinFlux is None) or (avgFlux < curRadiusMinFlux):
                    curRadiusMinFlux = avgFlux
                if (curRadiusMaxFLux is None) or (avgFlux > curRadiusMaxFLux):
                    curRadiusMaxFLux = avgFlux
        
        
        # LOOP 2 (USES THE PREV LOOPS RESULTS OF MIN/MAX TO SCALE THE VALUES)
        majAxisIndex = int( (current_aepObj.majorAxisLen - minMajAxLen) + 2)
        for theta,pixelList in uniqueNeighborThetas:
            totalAvgFlux = 0
            for i,j in pixelList:
                avgFlux = createStructs.calcFlux(i,j,fits1,fits2)
                totalAvgFlux += avgFlux
            totalAvgFlux = totalAvgFlux/len(pixelList)
            thetaIndex = (newOverallMaxTheta - theta) if sWise else (theta - newOverallMinTheta)
            relScale = abs(avgFlux-curRadiusMinFlux) / abs(curRadiusMaxFLux-curRadiusMinFlux)
            if relScale == 1:            # If its equal to the max, tone down a little so can still see it
                relScale -= 0.01
            FINALPLOT[majAxisIndex,thetaIndex] = relScale



    # PLOTTING ###############################################################################################################################
    fig,ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [3, 1]})
    ### PLOT 1
    ax[0].imshow(FINALPLOT, origin="lower", extent = [0, newOverallMaxTheta-newOverallMinTheta+1, minMajAxLen-2,maxMajAxLen+2])
    ax[0].set_aspect(2)
    ax[0].set_xlabel("θ from front", size=13)
    ax[0].set_ylabel(f"Major Axis Length ({minMajAxLen}-{maxMajAxLen})", size=13)
    ### PLOT 2
    armReference = np.zeros((rows,cols))
    for i,j in armsPixels:
        armReference[i,j] = 1
    frontStartTheta = newOverallMaxTheta if sWise else newOverallMinTheta
    for majaxLen in range(1,int(maxMajAxLen/2)):
        i,j = createStructs.calcElpsPoint(majaxLen, majaxLen*minMaxRatio, axisRadians, frontStartTheta, (inputCenterR, inputCenterR))
        armReference[i,j] = 0.5
    ax[1].imshow(armReference)
    plt.suptitle(f"{galNum}_({waveband1}-{waveband2})_merge({merge})", size=17)
    

    plt.savefig(f"runRand200/{galNum}_({waveband1}-{waveband2})_merge({merge}).pdf")
    # plt.show()




if __name__ == "__main__":

    count = 0
    galFile = open('getItSp.txt','r')
    failedGalaxys = open("DEBUGLATER.txt","w")
    galaxy = galFile.readline()
    while galaxy:
        print("-------------------------------------------")
        print(galaxy)
        try:
            main(merge=0, waveband1='g',waveband2='i',galNum=galaxy.rstrip("\n"))
        except:
            print("NEED TO DEBUG THIS FAILURE")
            failedGalaxys.write(galaxy)
        galaxy = galFile.readline()
       if count == 200:
           break
       count += 1
    galFile.close()
    failedGalaxys.close()

    # main(merge=0, waveband1='g',waveband2='i',galNum="1237660635996291172")
    # main(merge=1, waveband1='g',waveband2='i',galNum="1237660635996291172")
    # main(merge=2, waveband1='g',waveband2='i',galNum="1237660635996291172")

    # main(merge=0, waveband1='g',waveband2='i',galNum="1237648705658486867")
    # main(merge=1, waveband1='g',waveband2='i',galNum="1237648705658486867")
    # main(merge=2, waveband1='g',waveband2='i',galNum="1237648705658486867")
    