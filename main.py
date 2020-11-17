from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math
from statistics import mean


import inputFiles
import createStructs



class ThresholdError(Exception):
    pass



def main(waveband1, waveband2, galNum):
    GROUP = True

    minMaxRatio, minorAxis, majorAxis, axisRadians, inputCenterR, inputCenterC = inputFiles.read_arcsTSV(galNum)
    pixelLoc, rows, cols       = inputFiles.read_clusMask(waveband=waveband1,galNum=galNum,group=GROUP)
    armsPixels    = createStructs.get_largestArm(galaxyArms=pixelLoc)

    arcsEllipse_Positions = createStructs.arm_to_ArcsEllipse(majorAxis=majorAxis, minMaxRatio=minMaxRatio, axisRadians=axisRadians, armsPixels=armsPixels, center=(inputCenterR, inputCenterR))
    for i in range(1,len(arcsEllipse_Positions)-1):
        if (len(arcsEllipse_Positions[i].arc) > 30) and (len(arcsEllipse_Positions[i].arc) < 270):
            curEllipseInfo = arcsEllipse_Positions[i]
            thetaList = createStructs.simThetas(curEllipseInfo,arcsEllipse_Positions[i-1],arcsEllipse_Positions[i+1])### WHAT IF BEFORE/AFTER is 4 | 6 | 7, not next to each other            
            merged = inputFiles.mergedFits(thetaList, waveband1,waveband2,galNum)
            mergedTheta = []
            mergedDifs = []
            for theta, dif in merged: # adjust phi first, then consider overlap
                mergedTheta.append(theta)
                mergedDifs.append(dif)
            unmerged = inputFiles.unmergedFits(curEllipseInfo.arc,waveband1,waveband2,galNum)
            unmergedTheta = []
            unmergedDifs = []
            for theta, dif in unmerged:
                unmergedTheta.append(theta)
                unmergedDifs.append(dif)

            mng = plt.get_current_fig_manager()
            mng.window.state('zoomed')
            plt.subplot(121)
            a = np.zeros((rows,cols))
            for ii, ij in armsPixels:
                a[ii,ij]=2
            for ii, ij in curEllipseInfo.ellipse:
                a[ii,ij]=1
            for phi,ii, ij in curEllipseInfo.arc:
                a[ii,ij]=3
            plt.imshow(a) # origin='lower' for "normal" X/Y axis position
            plt.text(rows*0.6,cols*0.85,"Min θ: {}".format(curEllipseInfo.minTheta),color="white")
            plt.text(rows*0.6,cols*0.95,"Max θ: {}".format(curEllipseInfo.maxTheta),color="white")
            plt.title('Combination')

            plt.subplot(222)
            plt.plot(unmergedTheta,unmergedDifs)
            plt.title("Unmerged MajAx {}, {}-{}".format(curEllipseInfo.majorAxisLen,waveband1,waveband2))
            plt.subplot(224)

            plt.plot(mergedTheta,mergedDifs)
            plt.xlabel("θ (degrees)",size=15)
            plt.ylabel("Flux (nanomaggies)",size=15)
            plt.title("Merged MajAx ({},{},{}), {}-{}".format(curEllipseInfo.majorAxisLen-2,curEllipseInfo.majorAxisLen,curEllipseInfo.majorAxisLen+2,waveband1,waveband2))
            plt.subplots_adjust(hspace=0.3,wspace=0.4)
            plt.suptitle("MajorAxis: {}".format(curEllipseInfo.majorAxisLen), size=20)
            # plt.savefig("{}-_{}_{}-{}.pdf".format(galNum,curEllipseInfo.majorAxisLen,waveband1,waveband2))
            plt.show()
            plt.close()
            break


if __name__ == "__main__":



    galNum = "1237660635996291172"
    main("g", "r", galNum)