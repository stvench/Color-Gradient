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
    arcsEllipse_Positions = createStructs.arm_to_ArcsEllipse(majorAxis=majorAxis, minMaxRatio=minMaxRatio, axisRadians=axisRadians, armsPixels=armsPixels, center=(inputCenterR, inputCenterR))



    # NEED MINTHETA / MAXTHETA / MIN MAJAXLEN / MAX MAJAXLEN
        # COULD DO IT DIRECTLY IN CREATESTRUCTS-- MAYBE LATER
    # minTheta    = None
    # maxTheta    = None
    # minMajAxLen = None
    # maxMajAxLen = None
    # for aepObj in arcsEllipse_Positions:
    #     if (aepObj.arc):
    #         # len check
    #         if (minMajAxLen is None) or (aepObj.majorAxisLen < minMajAxLen):
    #             minMajAxLen = aepObj.majorAxisLen
    #         if (maxMajAxLen is None) or (aepObj.majorAxisLen > maxMajAxLen):
    #             maxMajAxLen = aepObj.majorAxisLen
    #         # theta check
    #         for theta, i, j in aepObj.arc:
    #             if (minTheta is None) or (theta < minTheta):
    #                 minTheta = theta
    #             if (maxTheta is None) or (theta > maxTheta):
    #                 maxTheta = theta
    # print(minTheta, maxTheta, minMajAxLen, maxMajAxLen)
    
    # y = np.arange(minMajAxLen-1,maxMajAxLen+2)
    # x = np.arange(minTheta,maxTheta+1)
    # z = np.zeros((len(y),len(x)))
    # filePath1 = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/1. Circle/Sample_Galaxies/{galNum}/{galNum}_{waveband1}.fits"
    # hdul = fits.open(filePath1)
    # a = hdul[0].data # 2-D Numpy Array
    # hdul.close()
    # filePath2 = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/1. Circle/Sample_Galaxies/{galNum}/{galNum}_{waveband2}.fits"
    # hdul = fits.open(filePath2)
    # b = hdul[0].data # 2-D Numpy Array
    # hdul.close()
    # minFlux = None
    # maxFlux = None
    # for aepObj in arcsEllipse_Positions:
    #     if (aepObj.arc):
    #         majAxisIndex = aepObj.majorAxisLen - minMajAxLen+1
    #         for theta,i,j in aepObj.arc:
    #                 flux = 0
    #                 if i==0:
    #                     if j==0:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j+1] - b[i,j+1]

    #                         flux += a[i+1,j] - b[i+1,j]
    #                         flux += a[i+1,j+1] - b[i+1,j+1]
    #                         avgFlux = flux/4
    #                     elif j==len(a)-1:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j-1] - b[i,j-1]

    #                         flux += a[i+1,j] - b[i+1,j]
    #                         flux += a[i+1,j-1] - b[i+1,j-1]
    #                         avgFlux = flux/4
    #                     else:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j+1] - b[i,j+1]
    #                         flux += a[i,j-1] - b[i,j-1]

    #                         flux += a[i+1,j] - b[i+1,j]
    #                         flux += a[i+1,j+1] - b[i+1,j+1]
    #                         flux += a[i+1,j-1] - b[i+1,j-1]
    #                         avgFlux = flux/6
    #                 elif i==len(a)-1:
    #                     if j==0:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j+1] - b[i,j+1]

    #                         flux += a[i-1,j] - b[i-1,j]
    #                         flux += a[i-1,j+1] - b[i-1,j+1]
    #                         avgFlux = flux/4
    #                     elif j==len(a)-1:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j-1] - b[i,j-1]

    #                         flux += a[i-1,j] - b[i-1,j]
    #                         flux += a[i-1,j-1] - b[i-1,j-1]
    #                         avgFlux = flux/4
    #                     else:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j+1] - b[i,j+1]
    #                         flux += a[i,j-1] - b[i,j-1]

    #                         flux += a[i-1,j] - b[i-1,j]
    #                         flux += a[i-1,j+1] - b[i-1,j+1]
    #                         flux += a[i-1,j-1] - b[i-1,j-1]
    #                         avgFlux = flux/6
    #                 elif j==0:
    #                     if i==0:
    #                         # Already accounted for
    #                         pass
    #                     elif i==len(a)-1:
    #                         # Already accounted for
    #                         pass
    #                     else:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j+1] - b[i,j+1]

    #                         flux += a[i+1,j] - b[i+1,j]
    #                         flux += a[i+1,j+1] - b[i+1,j+1]

    #                         flux += a[i-1,j] - b[i-1,j]
    #                         flux += a[i-1,j+1] - b[i-1,j+1]
    #                         avgFlux = flux/6
    #                 elif j==len(a)-1:
    #                     if i==0:
    #                         # Already accounted for
    #                         pass
    #                     elif i==len(a)-1:
    #                         # Already accounted for
    #                         pass
    #                     else:
    #                         flux += a[i,j] - b[i,j]
    #                         flux += a[i,j-1] - b[i,j-1]

    #                         flux += a[i+1,j] - b[i+1,j]
    #                         flux += a[i+1,j-1] - b[i+1,j-1]

    #                         flux += a[i-1,j] - b[i-1,j]
    #                         flux += a[i-1,j-1] - b[i-1,j-1]
    #                         avgFlux = flux/6
    #                 else:
    #                     flux += a[i,j] - b[i,j]
    #                     flux += a[i,j+1] - b[i,j+1]
    #                     flux += a[i,j-1] - b[i,j-1]

    #                     flux += a[i+1,j] - b[i+1,j]
    #                     flux += a[i+1,j+1] - b[i+1,j+1]
    #                     flux += a[i+1,j-1] - b[i+1,j-1]

    #                     flux += a[i-1,j] - b[i-1,j]
    #                     flux += a[i-1,j+1] - b[i-1,j+1]
    #                     flux += a[i-1,j-1] - b[i-1,j-1]
    #                     avgFlux = flux/9
    #                 # Min/Max checks
    #                 if (minFlux is None) or (avgFlux < minFlux):
    #                     minFlux = avgFlux
    #                 if (maxFlux is None) or (avgFlux > maxFlux):
    #                     maxFlux = avgFlux
    #                 thetaIndex = theta-minTheta
    #                 z[majAxisIndex,thetaIndex] = avgFlux
    # for i in range(len(y)):
    #     for j in range(len(x)):
    #         if z[i,j] == 0:
    #             z[i,j] = maxFlux

    # plt.contourf(x,y,z,cmap=cm.gray)
    # plt.colorbar()
    # plt.xlabel("θ starting from min", size=15)
    # plt.ylabel("Major Axis Length", size =15)
    # plt.show()



if __name__ == "__main__":


    main(waveband1='g',waveband2='i',galNum="1237660635996291172")