# External Imports
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull, convex_hull_plot_2d


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

    # Initialize variables for plotting
    minSemiMajAxLen = int(arcsEllipse_Positions[0].majorAxisLen/2)
    maxSemiMajAxLen = int(arcsEllipse_Positions[-1].majorAxisLen/2)
    xRange = newOverallMaxTheta-newOverallMinTheta+1     # Columns
    yRange = maxSemiMajAxLen-minSemiMajAxLen+1           # Rows
    FINALPLOT = np.zeros((yRange,xRange,3), dtype=float) # array of [1.0,1.0,1.0] (plt.imshow() for RBGs floats bounded by [0...1], 1==WHITE)
    fits1 = inputFiles.readFits(waveband1,galNum,onOpenlabs)
    fits2 = inputFiles.readFits(waveband2,galNum,onOpenlabs)

    # Make bins with range of min/max flux
    nBins = 100
    overallMinFlux, overallMaxFlux = createStructs.calcOverallMinMaxFlux(arcsEllipse_Positions=arcsEllipse_Positions,merge=merge,fits1=fits1,fits2=fits2)
    bins = np.linspace(overallMinFlux,overallMaxFlux,nBins)

    # Calculate fluxs into FINALPLOT
    for ii in range(merge,len(arcsEllipse_Positions)-merge):
        current_aepObj = arcsEllipse_Positions[ii]
        neighbor_aepObjs = [arcsEllipse_Positions[j] for j in np.arange(ii-merge,ii+merge+1) if (ii!=j)]
        groupedThetas = createStructs.groupNeighborThetas(middle=current_aepObj, neighbors=neighbor_aepObjs)
        semiMajAxisIndex = int(current_aepObj.majorAxisLen/2) - minSemiMajAxLen
        for theta,pixelList in groupedThetas:
            totalAvgFlux = 0
            for i,j in pixelList:
                avgFlux = createStructs.calcFlux(i,j,fits1,fits2)
                totalAvgFlux += avgFlux
            totalAvgFlux /= len(pixelList)
            thetaIndex = (newOverallMaxTheta - theta) if sWise else (theta - newOverallMinTheta)
            for i in range(nBins-1):
                # if the totalAvgFlux falls in between a bin range
                if (bins[i]<=totalAvgFlux and totalAvgFlux<=bins[i+1]) or (bins[i]>=totalAvgFlux and totalAvgFlux>=bins[i+1]):
                    break
            relScale = 1-i/nBins
            FINALPLOT[semiMajAxisIndex,thetaIndex] = relScale # Automatically converts relScale -> [relScale, relScale, relScale]

    # Get all points in FINALPLOT that are non-zero, put in np array, compute CONVEX HULL
    coords = []
    for row in range(0,FINALPLOT.shape[0]):
        for col in range(0,FINALPLOT.shape[1]):
            if np.all( FINALPLOT[row,col]!=0 ):
                coords.append( [col,row] )          # col==x, row ==y
            # GET POINTS THAT ARE ON THE OUTLINE OF THE PLOT plot lines as reference
            else:
                # TODO check border values
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                pass
    originalcoords = np.array(coords)
    originalHull   = ConvexHull(originalcoords)
    # Reverse compute the pixels within convex hull ( TODO: currently does not do the grouped thetas, so it doesn't takethe average of fluxs nearby, could maybe combine this loop with the previous one that utilizes groupedThetas)
    for row in range(0,FINALPLOT.shape[0]):
        for col in range(0,FINALPLOT.shape[1]):
            newcoords = coords.copy()
            newcoords.append( [col,row] )           # col==x, row ==y 
            newHull   = ConvexHull(np.array(newcoords))
            if (np.array_equal(originalHull.simplices,newHull.simplices)):
                # Reverse the adjustments made to the SemiMajAxisLen and Theta
                semiMajAxis  = row + minSemiMajAxLen
                theta        = (newOverallMaxTheta - col) if sWise else (col + newOverallMinTheta)
                fitsI, fitsJ = createStructs.calcElpsPoint(a=semiMajAxis, b=semiMajAxis*minMaxRatio, axisRadians=axisRadians, curTheta=theta, center=(inputCenterR, inputCenterR))
                avgFlux = createStructs.calcFlux(fitsI,fitsJ,fits1,fits2)
                for i in range(nBins-1):
                    if (bins[i]<=avgFlux and avgFlux<=bins[i+1]) or (bins[i]>=avgFlux and avgFlux>=bins[i+1]):
                        break
                relScale = 1-i/nBins
                FINALPLOT[row,col] = relScale # Automatically converts relScale -> [relScale, relScale, relScale]

   

    # PLOTTING  the WHITER it is, the larger WAVEBAND2 is.
               ### WHITE means the MINIMUM difference, meaning WAVEBAND2 is at its largest
               ### if negative range, waveband2>waveband1, if positive range, waveband2<waveband1
               ###     Take the range of flux, min(0) to max(100), and take current flux(5). To get i, iterate 
               ###     through the range until current is in between one of the ranges. Take that i, which in this 
               ###     case is very small(5/100 == 5) and put into relScale = 1-(5/100) to get 0.95, which is WHITE
    fig,ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [3, 1]})
    ### PLOT 1 (actual color difference)
    ax[0].imshow(FINALPLOT, origin="lower", extent = [0, newOverallMaxTheta-newOverallMinTheta+1, 0,maxSemiMajAxLen-minSemiMajAxLen+1])
    ax[0].set_aspect(2)
    ax[0].set_xlabel("θ from front", size=13)
    ax[0].set_ylabel(f"Semi-Major Axis Length ({0}-{maxSemiMajAxLen-minSemiMajAxLen})", size=13)
    for simplex in originalHull.simplices:
        ax[0].plot(originalcoords[simplex, 0], originalcoords[simplex, 1], '--',color ='red')
    ### PLOT 2
    imgAPng = inputFiles.read_imageAPng(waveband=waveband1,galNum=galNum,onOpenlabs=onOpenlabs)
    outlineColor = np.max(imgAPng)
    # Plot the arms outline
    armsPixelsOutline = createStructs.armOutline(armsPixels=armsPixels)
    for i,j in armsPixelsOutline:
        imgAPng[i,j] = outlineColor
    frontStartTheta = newOverallMaxTheta if sWise else newOverallMinTheta
    # Plot the "front" of the arm as a line
    for semiMajAxLen in range(1,maxSemiMajAxLen):
        i,j = createStructs.calcElpsPoint(semiMajAxLen, semiMajAxLen*minMaxRatio, axisRadians, frontStartTheta, (inputCenterR, inputCenterC))
        imgAPng[i,j] = outlineColor
    # Plot the overall shape of galaxy for ellipse reference
    for i,j in arcsEllipse_Positions[-1].ellipse:
        imgAPng[i,j] = outlineColor
    ax[1].imshow(imgAPng)

    plt.suptitle(f"{galNum}_({waveband1}-{waveband2})_merge({merge})", size=17)
    if (makePDF):
        plt.savefig(f"runs/{galNum}_({waveband1}-{waveband2})_merge({merge}).pdf")
    else:
        plt.show()
    plt.close()




if __name__ == "__main__":




    onOpenlabs = False
    makePDF    = False




    if (onOpenlabs):
        ### Openlabs runs

        ### 1. Should create a runs folder
        ### 2. Already needs the data folder
        count    = 1
        errCount = 0
        galFile       = open('data/getItSp.txt','r')           # Each line is a galaxy to test
        failedGalaxys = open("data/DEBUGLATER.txt","w")  # Each line is a galaxy that failed
        galaxy = galFile.readline()
        while galaxy:
            print("-------------------------------------------")
            print(count,galaxy)
            try:
                main(merge=0, waveband1='g', waveband2='i', galNum=galaxy.rstrip("\n"), onOpenlabs=onOpenlabs, makePDF=makePDF)
            except KeyboardInterrupt:
                break
            except:
                failedGalaxys.write(galaxy)
                errCount += 1
                # raise 
            galaxy = galFile.readline()
            if count == 200:
                break
            count += 1
        galFile.close()
        failedGalaxys.close()
        print(f"{count-errCount}/{count} ({(count-errCount)/count*100:.2f}%) succeeded")
    else:
        ### Local runs
        #  1172 get_PositionArm Arm: Extra1(i=127,j=137)    Arm Extra2(i=127,j=147)
        main(merge=0, waveband1='g',waveband2='i',galNum="1237660635996291172",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=1, waveband1='g',waveband2='i',galNum="1237660635996291172",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=2, waveband1='g',waveband2='i',galNum="1237660635996291172",onOpenlabs=onOpenlabs,makePDF=makePDF)

        # main(merge=0, waveband1='g',waveband2='i',galNum="1237648705658486867",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=1, waveband1='g',waveband2='i',galNum="1237648705658486867",onOpenlabs=onOpenlabs,makePDF=makePDF)
        # main(merge=2, waveband1='g',waveband2='i',galNum="1237648705658486867",onOpenlabs=onOpenlabs,makePDF=makePDF)
    
