from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


import inputFiles
import createStructs


def main(waveband1, waveband2, galNum):
    minMaxRatio, minorAxis, majorAxis, axisRadians, inputCenterR, inputCenterC = inputFiles.read_arcsTSV(galNum)
    ############### 2.
    pixelLoc, rows, cols = inputFiles.read_clusMask(waveband=waveband1,galNum=galNum,group=True)
    ############### 3.
    armsPixels = createStructs.get_largestArm(galaxyArms=pixelLoc)
    ############### 4.
    arcsEllipse_Positions, overallMinTheta, overallMaxTheta = createStructs.arm_to_ArcsEllipse(majorAxis=majorAxis, minMaxRatio=minMaxRatio, axisRadians=axisRadians, armsPixels=armsPixels, center=(inputCenterR, inputCenterR))



    # Loop 1
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
    # print(newOverallMinTheta,newOverallMaxTheta, minMajAxLen, maxMajAxLen)

    
    x = np.arange(0,newOverallMaxTheta-newOverallMinTheta+1)
    y = np.arange(minMajAxLen/2-1,maxMajAxLen/2+2)
    z = np.zeros((len(y),len(x)))
    a = inputFiles.readFits(waveband1,galNum)
    b = inputFiles.readFits(waveband2,galNum)
    minFlux = None
    maxFlux = None
    # Loop 2
    for aepObj in arcsEllipse_Positions:
        if (aepObj.arc):
            majAxisIndex = int( (aepObj.majorAxisLen - minMajAxLen)/2 +1)
            for theta,i,j in aepObj.arc:
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
                # Min/Max checks
                if (minFlux is None) or (avgFlux < minFlux):
                    minFlux = avgFlux
                if (maxFlux is None) or (avgFlux > maxFlux):
                    maxFlux = avgFlux
                if sWise:
                    thetaIndex = newOverallMaxTheta - theta
                else:
                    thetaIndex = theta - newOverallMinTheta
                z[majAxisIndex,thetaIndex] = avgFlux
    for i in range(len(y)):
        for j in range(len(x)):
            if z[i,j] == 0:
                z[i,j] = maxFlux

    plt.contourf(x,y,z,cmap=cm.gray)
    plt.colorbar()
    plt.xlabel("θ starting from min", size=15)
    plt.ylabel("Major Axis Length", size =15)
    plt.show()



if __name__ == "__main__":


    main(waveband1='g',waveband2='i',galNum="1237660635996291172")

    # main(waveband1='g',waveband2='i',galNum="1237648705658486867")
    