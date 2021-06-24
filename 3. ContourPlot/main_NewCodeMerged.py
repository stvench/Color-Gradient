# External Imports
from astropy.io import fits
import matplotlib.pyplot as plt


# Internal Imports
import inputFiles
import createStructs


def main(merge, waveband1, waveband2, galNum, onOpenlabs, makePDF):
    # Get all important info from arcs tsv for the galaxy
    minMaxRatio, minorAxis, majorAxis, axisRadians, inputCenterR, inputCenterC = inputFiles.read_arcsTSV(galNum,onOpenlabs)
    # Get the clustermasks for each waveband, confirm dimensions are equal
    pixelLoc1, rows1, cols1 = inputFiles.read_clusMask(waveband=waveband1,galNum=galNum,group=True,onOpenlabs=onOpenlabs)
    pixelLoc2, rows2, cols2 = inputFiles.read_clusMask(waveband=waveband2,galNum=galNum,group=True,onOpenlabs=onOpenlabs)
    rows = rows1 if rows1==rows2 else None
    cols = cols1 if cols1==cols2 else None

    # Use get_largestArm() for auto || get_PositionArm() for manual selection, looking at clusMask position of arms
    armsPixels1 = createStructs.get_largestArm(galaxyArms=pixelLoc1)
    ######## armsPixels1 = createStructs.get_PositionArm(galaxyArms=pixelLoc1,i=127,j=147)

    armsPixels  = createStructs.unionClosestArm(waveband1LargestArm=armsPixels1, waveband2AllArms=pixelLoc2)
    # Create the arcsEllipse_Positions & update thetas
    arcsEllipse_Positions, overallMinTheta, overallMaxTheta = createStructs.arm_to_ArcsEllipse(majorAxis=majorAxis, minMaxRatio=minMaxRatio, axisRadians=axisRadians, armsPixels=armsPixels, center=(inputCenterR, inputCenterR))
    newOverallMinTheta, newOverallMaxTheta, sWise = createStructs.updateThetaStarts(overallMinTheta=overallMinTheta, overallMaxTheta=overallMaxTheta, arcsEllipse_Positions=arcsEllipse_Positions)

    for i in range(1,len(arcsEllipse_Positions)-1):
        curEllipseInfo = arcsEllipse_Positions[i]
        filePath1 = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/3. ContourPlot/data/Sample_Galaxies/{galNum}/{galNum}_{waveband1}.fits"
        hdul = fits.open(filePath1)
        image_data1 = hdul[0].data # 2-D Numpy Array
        hdul.close()
        filePath2 = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/3. ContourPlot/data/Sample_Galaxies/{galNum}/{galNum}_{waveband2}.fits"
        hdul = fits.open(filePath2)
        image_data2 = hdul[0].data # 2-D Numpy Array
        hdul.close()
        dataVals = [(phi,image_data1[i,j]-image_data2[i,j]) for phi,i,j in (curEllipseInfo.arc)]

        theta = []
        difs = []
        for t,d in dataVals:
            theta.append(t)
            difs.append(d)
        # print(dataVals)
        plt.title(f"{galNum}_({waveband1}-{waveband2})_merge({merge})", size=17)
        plt.plot(theta,difs)


        if (makePDF):
            plt.savefig(f"runs/{galNum}_({waveband1}-{waveband2})_merge({merge}).pdf")
        else:
            plt.show()
        plt.close()




if __name__ == "__main__":

    onOpenlabs = True
    makePDF    = True

    if (onOpenlabs):
        ### # Openlabs Single Run of specific galaxy/wavebands
        main(merge=0, waveband1='g', waveband2='i', galNum="1237660635996291172", onOpenlabs=onOpenlabs, makePDF=makePDF)



        ### # Openlabs Batch Run of "data/Steven-Big-Galaxies.txt"
        # count    = 1
        # errCount = 0
        # galFile       = open('data/Steven-Big-Galaxies.txt','r')           # Each line is a galaxy to test
        # failedGalaxys = open("data/DEBUGLATER.txt","w")                    # Each line is a galaxy that failed
        # galaxy = galFile.readline()
        # while galaxy:
        #     print("-------------------------------------------")
        #     print(count,galaxy)
        #     try:
        #         main(merge=0, waveband1='u', waveband2='z', galNum=galaxy.rstrip("\n"), onOpenlabs=onOpenlabs, makePDF=makePDF)
        #     except KeyboardInterrupt:
        #         break
        #     except:
        #         failedGalaxys.write(galaxy)
        #         errCount += 1
        #         # raise 
        #     galaxy = galFile.readline()
        #     # if count == 200:
        #     #     break
        #     count += 1
        # galFile.close()
        # failedGalaxys.close()
        # print(f"{count-errCount}/{count} ({(count-errCount)/count*100:.2f}%) succeeded")

    else:
        ### Personal local runs
        #  1172 get_PositionArm Arm: Extra1(i=127,j=137)    Arm Extra2(i=127,j=147)
        main(merge=0, waveband1='g',waveband2='i',galNum="1237660635996291172",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=1, waveband1='g',waveband2='i',galNum="1237660635996291172",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=2, waveband1='g',waveband2='i',galNum="1237660635996291172",onOpenlabs=onOpenlabs,makePDF=makePDF)

        # main(merge=0, waveband1='g',waveband2='i',galNum="1237648705658486867",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=1, waveband1='g',waveband2='i',galNum="1237648705658486867",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=2, waveband1='g',waveband2='i',galNum="1237648705658486867",onOpenlabs=onOpenlabs,makePDF=makePDF)
    
