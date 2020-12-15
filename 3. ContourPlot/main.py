from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
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



    # LOOP 1 (FIXES THE THETAS TO START AT 0 from FRONT)
    needSub360 = True if overallMaxTheta-overallMinTheta >= 360 else False
    newOverallMinTheta = None
    newOverallMaxTheta = None
    for ae_pObj in arcsEllipse_Positions:
        updatedThetas = []
        for theta,i,j in ae_pObj.arc:
            if (needSub360):
                if theta>360:
                    theta-=360
            if newOverallMinTheta is None or theta<newOverallMinTheta:
                newOverallMinTheta = theta
            if newOverallMaxTheta is None or theta>newOverallMaxTheta:
                newOverallMaxTheta = theta
            updatedThetas.append( (theta,i,j) )
        ae_pObj.arc = updatedThetas

    frontTheta = arcsEllipse_Positions[0].arc[0][0]
    endTheta   = arcsEllipse_Positions[-1].arc[-1][0]
    sWise = None 
    if frontTheta>endTheta: # Counter clockwise
        sWise = True
    else:
        sWise = False

    minMajAxLen = arcsEllipse_Positions[0].majorAxisLen
    maxMajAxLen = arcsEllipse_Positions[-1].majorAxisLen
    x = np.arange(0,newOverallMaxTheta-newOverallMinTheta+1)
    y = np.arange(minMajAxLen-2,maxMajAxLen+2)
    fig, ax = plt.subplots()
    FINALPLOT = np.ones((len(y),len(x),3), dtype=float)
    
    a = inputFiles.readFits(waveband1,galNum)
    b = inputFiles.readFits(waveband2,galNum)
    # LOOP 2
    for i in range(merge,len(arcsEllipse_Positions)-merge):
        current_aepObj = arcsEllipse_Positions[i]
        # print(np.arange(i-merge,i+merge+1))

        neighbor_aepObj = [arcsEllipse_Positions[j] for j in np.arange(i-merge,i+merge+1) if (i!=j)]
        uniqueNeighborThetas = createStructs.remvSimThetas(middle=current_aepObj, neighbors=neighbor_aepObj)
        # print(uniqueNeighborThetas)

        curRadiusMinFlux = None
        curRadiusMaxFLux = None
        # LOOP 1 (GETS THE RELATIVE MIN/MAX FOR THE CURRENT RADIUS)
        for theta,pixelList in uniqueNeighborThetas:
            for i,j in pixelList:
                flux = 0
                if i==0:
                    if j==0:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j+1] - b[i+1,j+1]
                        avgFlux = flux/4
                    elif j==len(a)-1:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j-1] - b[i+1,j-1]
                        avgFlux = flux/4
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j+1] - b[i+1,j+1]
                        flux += a[i+1,j-1] - b[i+1,j-1]
                        avgFlux = flux/6
                elif i==len(a)-1:
                    if j==0:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j+1] - b[i-1,j+1]
                        avgFlux = flux/4
                    elif j==len(a)-1:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j-1] - b[i-1,j-1]
                        avgFlux = flux/4
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j+1] - b[i-1,j+1]
                        flux += a[i-1,j-1] - b[i-1,j-1]
                        avgFlux = flux/6
                elif j==0:
                    if i==0:
                        # Already accounted for
                        pass
                    elif i==len(a)-1:
                        # Already accounted for
                        pass
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j+1] - b[i+1,j+1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j+1] - b[i-1,j+1]
                        avgFlux = flux/6
                elif j==len(a)-1:
                    if i==0:
                        # Already accounted for
                        pass
                    elif i==len(a)-1:
                        # Already accounted for
                        pass
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j-1] - b[i+1,j-1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j-1] - b[i-1,j-1]
                        avgFlux = flux/6
                else:
                    flux += a[i,j] - b[i,j]
                    flux += a[i,j+1] - b[i,j+1]
                    flux += a[i,j-1] - b[i,j-1]

                    flux += a[i+1,j] - b[i+1,j]
                    flux += a[i+1,j+1] - b[i+1,j+1]
                    flux += a[i+1,j-1] - b[i+1,j-1]

                    flux += a[i-1,j] - b[i-1,j]
                    flux += a[i-1,j+1] - b[i-1,j+1]
                    flux += a[i-1,j-1] - b[i-1,j-1]
                    avgFlux = flux/9
                # CURRENT INNER RADIUS Min/Max checks
                if (curRadiusMinFlux is None) or (avgFlux < curRadiusMinFlux):
                    curRadiusMinFlux = avgFlux
                if (curRadiusMaxFLux is None) or (avgFlux > curRadiusMaxFLux):
                    curRadiusMaxFLux = avgFlux
        majAxisIndex = int( (current_aepObj.majorAxisLen - minMajAxLen) + 2)

        # LOOP 2 (USES THE PREV LOOPS RESULTS OF MIN/MAX TO SCALE THE VALUES)
        for theta,pixelList in uniqueNeighborThetas:
            totalAvgFlux = 0
            for i,j in pixelList:
                flux = 0
                if i==0:
                    if j==0:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j+1] - b[i+1,j+1]
                        avgFlux = flux/4
                    elif j==len(a)-1:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j-1] - b[i+1,j-1]
                        avgFlux = flux/4
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j+1] - b[i+1,j+1]
                        flux += a[i+1,j-1] - b[i+1,j-1]
                        avgFlux = flux/6
                elif i==len(a)-1:
                    if j==0:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j+1] - b[i-1,j+1]
                        avgFlux = flux/4
                    elif j==len(a)-1:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j-1] - b[i-1,j-1]
                        avgFlux = flux/4
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j+1] - b[i-1,j+1]
                        flux += a[i-1,j-1] - b[i-1,j-1]
                        avgFlux = flux/6
                elif j==0:
                    if i==0:
                        # Already accounted for
                        pass
                    elif i==len(a)-1:
                        # Already accounted for
                        pass
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j+1] - b[i,j+1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j+1] - b[i+1,j+1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j+1] - b[i-1,j+1]
                        avgFlux = flux/6
                elif j==len(a)-1:
                    if i==0:
                        # Already accounted for
                        pass
                    elif i==len(a)-1:
                        # Already accounted for
                        pass
                    else:
                        flux += a[i,j] - b[i,j]
                        flux += a[i,j-1] - b[i,j-1]

                        flux += a[i+1,j] - b[i+1,j]
                        flux += a[i+1,j-1] - b[i+1,j-1]

                        flux += a[i-1,j] - b[i-1,j]
                        flux += a[i-1,j-1] - b[i-1,j-1]
                        avgFlux = flux/6
                else:
                    flux += a[i,j] - b[i,j]
                    flux += a[i,j+1] - b[i,j+1]
                    flux += a[i,j-1] - b[i,j-1]

                    flux += a[i+1,j] - b[i+1,j]
                    flux += a[i+1,j+1] - b[i+1,j+1]
                    flux += a[i+1,j-1] - b[i+1,j-1]

                    flux += a[i-1,j] - b[i-1,j]
                    flux += a[i-1,j+1] - b[i-1,j+1]
                    flux += a[i-1,j-1] - b[i-1,j-1]
                    avgFlux = flux/9
                totalAvgFlux += avgFlux
            totalAvgFlux = totalAvgFlux/len(pixelList)
            
            if sWise:
                thetaIndex = newOverallMaxTheta - theta
            else:
                thetaIndex = theta - newOverallMinTheta

            relScale = abs(avgFlux-curRadiusMinFlux) / abs(curRadiusMaxFLux-curRadiusMinFlux)
            if relScale == 1:            # If its equal to the max, tone down a little so can still see it
                relScale -= 0.01
            FINALPLOT[majAxisIndex,thetaIndex] = relScale


    ax.imshow(FINALPLOT, origin="lower", extent = [0, newOverallMaxTheta-newOverallMinTheta+1, minMajAxLen-2,maxMajAxLen+2])
    ax.set_aspect(2)
    plt.show()




if __name__ == "__main__":


    main(merge=1, waveband1='g',waveband2='i',galNum="1237660635996291172")

    # main(waveband1='g',waveband2='i',galNum="1237648705658486867")
    