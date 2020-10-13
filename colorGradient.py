from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math
from statistics import mean

def read_clusMask(waveband, galNum, group):
    # Get Image of Galaxy Arms
    imgClusMask = cv2.imread("1237660635996291172/{}/{}-K_clusMask-reprojected.png".format(waveband,galNum))
    dimensions = imgClusMask.shape
    rows = dimensions[0]
    cols = dimensions[1]
    # Store all Pixels of arms
    allPixelLoc = {}
    for i in range(rows):
        for j in range(cols):
            color = tuple(imgClusMask[i,j])
            if np.any(color != (0,0,0)):
                location = (i,j)
                if color in allPixelLoc:
                    allPixelLoc[color].add(location)
                else:
                    allPixelLoc[color] = {location}
    # All similar colors are added to the main color arms
    if group:
        for color, locs in allPixelLoc.items():
            if len(locs) <= 4:
                for i,j in locs:
                    borderColors = set()
                    D=imgClusMask[i+1,j]
                    if np.any(D != [0,0,0]):
                        borderColors.add(tuple(D))
                    U=imgClusMask[i-1,j]
                    if np.any(U != [0,0,0]):
                        borderColors.add(tuple(U))
                    R=imgClusMask[i,j+1]
                    if np.any(R != [0,0,0]):
                        borderColors.add(tuple(R))
                    L=imgClusMask[i,j-1]
                    if np.any(L != [0,0,0]):
                        borderColors.add(tuple(L))
                    DR=imgClusMask[i+1,j+1]
                    if np.any(DR != [0,0,0]):
                        borderColors.add(tuple(DR))
                    DL=imgClusMask[i+1,j-1]
                    if np.any(DL != [0,0,0]):
                        borderColors.add(tuple(DL))
                    UR=imgClusMask[i-1,j+1]
                    if np.any(UR != [0,0,0]):
                        borderColors.add(tuple(UR))
                    UL=imgClusMask[i-1,j-1]
                    if np.any(UL != [0,0,0]):
                        borderColors.add(tuple(UL))
                    # Remove any "small arm" colors
                    borderColors.discard(color)
                    for col in borderColors.copy():
                        if (len(allPixelLoc[col]) <= 4):
                            borderColors.remove(col)
                    # Calculate most similar, doesn't work if 1.(it is 2 pxiels away from nearest "large" arm) OR 2.(closest border color happens to be another border color)
                    if len(borderColors) != 0:
                        closestColor = None
                        colorDif = None
                        for c in borderColors:
                            curDif = sum([abs(int(color[index])-int(c[index])) for index in range(3)])
                            if (colorDif is None) or (curDif <= colorDif):
                                colorDif = curDif
                                closestColor = c
                        allPixelLoc[closestColor].add((i,j))
    # allPixelLoc(Dictionary) ->  key(color tuples) : value(set of (i,j location) tuples)
    # center(tuple) !math.ceil() should already return an INT, which is needed in arm_to_bestArc() range
    return allPixelLoc, (math.ceil(rows/2),math.ceil(cols/2)), rows, cols # {color : {(i,j) , (i,j)} } , (i,j) }

def get_singleArm(galaxyArms, position):
    ### Returns closest arm to the inputted position
    closestArm, minDist = None, None
    for arm in galaxyArms.values():
        for point in arm:
            curDist = calcDist(point,position)
            if (minDist is None) or (curDist<=minDist):
                minDist = curDist
                closestArm = arm
    return closestArm



def mergeArms(arms):
    ### arms should be a iterable of sets of tuples
    armPosition = set()
    for arm in arms:
        armPosition.update(arm)
    ### singleArm(single set of tuples)    
    return armPosition # {(i,j) , (i,j)}
    
def arm_to_ArcsCircles(armPosition,center):
    allDists = [calcDist(point,center) for point in armPosition]
    minRadius = math.ceil(min(allDists))
    maxRadius = math.floor(max(allDists))
    arcsCircles_Positions = []
    for radius in range(minRadius,maxRadius+1):
        arcsCircles_Positions.append( [radius,[],set(),None,None] )
        minPhi = None
        maxPhi = None
        # Obtain min/max angle of the arc overlap
        for phi in range(360):
            radian = phi/180*math.pi
            i = roundVal(radius*math.cos(radian)+center[0])
            j = roundVal(radius*math.sin(radian)+center[1])
            curDist = calcDist(center,(i,j))
            if (curDist >= radius-0.5) and (curDist <= radius+0.5) and ((i,j) in armPosition):
                if (minPhi is None):
                    minPhi = phi
                maxPhi = phi      
        arcsCircles_Positions[-1][3] = minPhi
        arcsCircles_Positions[-1][4] = maxPhi         
        # Stores the points for ARC and CIRCLE
        for phi in range(360):
            radian = phi/180*math.pi
            i = roundVal(radius*math.cos(radian)+center[0])
            j = roundVal(radius*math.sin(radian)+center[1])
            curDist = calcDist(center,(i,j))
            # If the current point is within 0.5 of the actual radius distance for the CIRCLE
            if (curDist >= radius-0.5) and (curDist <= radius+0.5):
                arcsCircles_Positions[-1][2].add((i,j))
                # If within the phi range for the ARC
                if (phi >= minPhi) and (phi <= maxPhi):
                    arcsCircles_Positions[-1][1].append((phi,i,j))
    # Returns a list of 5-tuples (radius, [arc], {circle}, minPhi, maxPhi)
    return arcsCircles_Positions



def t(middle, inner, outer):
    l = []
    for phi, i , j in middle[1]:
        l.append([(i,j)])
        for Iphi,Ii,Ij in inner[1]:
            if phi == Iphi:
                l[-1].append((Ii,Ij))
        for Ophi,Oi,Oj in outer[1]:
            if phi == Ophi:
                l[-1].append((Oi,Oj))
    return l

def unmergedFits(arcList,waveband1,waveband2,galNum):
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband1))
    image_data1 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband2))
    image_data2 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    return [image_data1[i,j]-image_data2[i,j] for phi,i,j in arcList]

def mergedFits(l, waveband1,waveband2,galNum):
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband1))
    image_data1 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband2))
    image_data2 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    return [mean([image_data1[i,j]-image_data2[i,j] for i,j in x]) for x in l]



def MAIN(waveband1, waveband2, galNum, position, group):
    pixelLoc1, center, rows, cols = read_clusMask(waveband=waveband1,galNum=galNum,group=group)
    waveband1Arm      = get_singleArm(galaxyArms=pixelLoc1, position=position)
    pixelLoc2, center, rows, cols = read_clusMask(waveband=waveband2,galNum=galNum,group=group)
    waveband2Arm      = get_singleArm(galaxyArms=pixelLoc2, position=position)

    armPosition                          = mergeArms(arms=[waveband1Arm,waveband2Arm])
    arcsCircles_Positions                = arm_to_ArcsCircles(armPosition=armPosition, center=center)

    for i in range(1,33):
        if len(arcsCircles_Positions[i][1]) > 50:
            unmerged = unmergedFits(arcsCircles_Positions[i][1],waveband1,waveband2,galNum)
            l = t(arcsCircles_Positions[i],arcsCircles_Positions[i-1],arcsCircles_Positions[i+1])
            merged = mergedFits(l, waveband1,waveband2,galNum)
            mng = plt.get_current_fig_manager()
            mng.window.state('zoomed')

            plt.subplot(121)
            a = np.zeros((rows,cols))
            for ii, ij in armPosition:
                a[ii,ij]=2
            for ii, ij in arcsCircles_Positions[i][2]:
                a[ii,ij]=1
            for phi,ii, ij in arcsCircles_Positions[i][1]:
                a[ii,ij]=3
            plt.imshow(a)
            plt.title('Combination')
            plt.subplot(222)
            plt.plot(unmerged)
            plt.title("Unmerged Radius {}, {}-{}".format(arcsCircles_Positions[i][0],waveband1,waveband2))
            plt.subplot(224)
            plt.plot(merged)
            plt.title("Merged Radius ({},{},{}), {}-{}".format(arcsCircles_Positions[i][0]-1,arcsCircles_Positions[i][0],arcsCircles_Positions[i][0]+1,waveband1,waveband2))
            plt.show()

    
### Helper Functions
def calcDist(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)  

def roundVal(val):
    if (val - int(val)) >= 0.5:
        return math.ceil(val)
    else:
        return math.floor(val)

### Debugging Function
def step_visualizer(rows, cols, waveband1, waveband2, galNum, position,pixelLoc1,pixelLoc2,waveband1Arm,waveband2Arm,armPosition,bestArcPosition,leftStartPoint,arcPixelOrdered,circle,subtracted):
        ### Plot waveband 1 
    plt.subplot(231)
    a = np.zeros((rows,cols))
    for color, pos in pixelLoc1.items():
        for i,j in pos:
            a[i,j]=2
    for i,j in waveband1Arm:
        a[i,j] = 1
    plt.imshow(a)
    plt.title(waveband1)
        ### Plot waveband 2 
    plt.subplot(232)
    a = np.zeros((rows,cols))
    for color, pos in pixelLoc2.items():
        for i,j in pos:
            a[i,j]=2
    for i,j in waveband2Arm:
        a[i,j] = 1
    plt.imshow(a)
    plt.title(waveband2)
        ### Plot Merged Arms
    plt.subplot(233)
    a = np.zeros((rows,cols))
    for i, j in armPosition:
        a[i,j]=2
    a[position[0],position[1]]=1
    plt.imshow(a)
    plt.title('{},{} arm merged at {}'.format(waveband1,waveband2,position))
        ### Plot waveband Best Arc and left startpoint
    plt.subplot(234)
    a = np.zeros((rows,cols))
    for i, j in armPosition:
        a[i,j]=2
    for i, j in circle:
        a[i,j]=1
    plt.imshow(a)
    plt.title('Circle + Arm')
        ### Plot Combination
    plt.subplot(235)
    a = np.zeros((rows,cols))
    for i, j in armPosition:
        a[i,j]=2
    for i, j in circle:
        a[i,j]=1
    for i, j in arcPixelOrdered:
        a[i,j]=3
    plt.imshow(a)
    plt.title('Combination')
        ### Plot wavebands color difference
    plt.subplot(236)
    plt.plot(subtracted)
    plt.title('{}-{} ({})'.format(waveband1.upper(),waveband2.upper(),galNum),size=20)
    plt.ylabel('Photon Count Dif',size=15)
    plt.xlabel('Position from Front',size=15)
    plt.suptitle("STEPS TAKEN")
    plt.subplots_adjust(left=0.05,right=0.95,top=0.92,wspace=0.30,hspace=0.35)
    mng = plt.get_current_fig_manager()
    mng.window.state('zoomed')
    plt.show()



if __name__ == "__main__":
    MAIN(waveband1='g',waveband2='i',galNum='1237660635996291172', position=(110,110), group=True)