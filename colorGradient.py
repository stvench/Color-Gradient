from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math


def read_clusMask(waveband, galNum, group):
    # Get Image of Galaxy Arms
    imgClusMask = cv2.imread("C:/Users/sc123/OneDrive/Desktop/gal/4. FITS/1237660635996291172/{}/{}-K_clusMask-reprojected.png".format(waveband,galNum))
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
    
def arm_Span(armPosition,center):
    allDists = [calcDist(point,center) for point in armPosition]
    minRadius = math.ceil(min(allDists))
    maxRadius = math.floor(max(allDists))
    #  Initialize minPhi, maxPhi
    minPhi = maxPhi = None
    for radius in range(minRadius,maxRadius+1):
        curRadiusMinPhi = curRadiusMaxPhi = None
        for phi in range(360): # Maybe directly range loop through radian?
            radian = phi/180*math.pi
            i = roundVal(radius*math.cos(radian)+center[0])
            j = roundVal(radius*math.sin(radian)+center[1])
            if (i,j) in armPosition:
                if curRadiusMinPhi is None:
                    curRadiusMinPhi = phi
                curRadiusMaxPhi = phi
        if (minPhi is None) or (curRadiusMinPhi <= minPhi):
            minPhi = curRadiusMinPhi
        if (maxPhi is None) or (curRadiusMaxPhi >= maxPhi):
            maxPhi = curRadiusMaxPhi
    ### Returns the min/max DEGREE value in range of [0,359]
    return minPhi, maxPhi, minRadius, maxRadius

def arm_to_ArcsCircles(armPosition,center,minPhi,maxPhi,minRadius,maxRadius):
    arcsCircles_Positions = []
    for radius in range(minRadius,maxRadius+1):
        arcsCircles_Positions.append( (radius,set(),set()) )
        for phi in range(360): ##################### Maybe directly range loop through radian?
            radian = phi/180*math.pi
            i = roundVal(radius*math.cos(radian)+center[0])
            j = roundVal(radius*math.sin(radian)+center[1])
            curDist = calcDist(center,(i,j))
            # If the current point is within 0.5 of the actual radius distance for the CIRCLE
            if (curDist >= radius-0.5) and (curDist <= radius+0.5):
                arcsCircles_Positions[-1][2].add((i,j))
                # If within the phi range for the ARC
                if (phi >= minPhi) and (phi <= maxPhi):
                    arcsCircles_Positions[-1][1].add((i,j))
    # Returns a list of 3-tuples (radius, {arc}, {circle})
    return arcsCircles_Positions


def arcStart(arcPositions):
    startPoints = []
    for i,j in arcPositions:
        # Get all neighbors of current pixel 
        neighbors = []
        right, left, up, down, mr, ml ,mu, md = False, False, False, False, False, False, False, False
        if (i,j+1) in arcPositions:
            neighbors.append((i,j+1))
            right = True
            mr = True
        if (i,j-1) in arcPositions:
            neighbors.append((i,j-1))
            left = True
            ml = True
        if (i+1,j) in arcPositions:
            neighbors.append((i+1,j))
            down = True
            md = True
        if (i-1,j) in arcPositions:
            neighbors.append((i-1,j))
            up = True
            mu = True
        if (i+1,j+1) in arcPositions:
            neighbors.append((i+1,j+1))
            down = True
            right = True
        if (i-1,j+1) in arcPositions:
            neighbors.append((i-1,j+1))
            up = True
            right = True
        if (i+1,j-1) in arcPositions:
            neighbors.append((i+1,j-1))
            down = True
            left = True
        if (i-1,j-1) in arcPositions:
            neighbors.append((i-1,j-1))
            up = True
            left = True
        # Check for "C" pattern
        if (right) and ( (not mu) and (not md) and (not left) ):
            startPoints.append((i,j))
        if (left)  and ( (not mu) and (not md) and (not right) ):
            startPoints.append((i,j))
        if (up)    and ( (not mr) and (not ml) and (not down) ):
            startPoints.append((i,j))
        if (down)  and ( (not mr) and (not ml) and (not up) ):
            startPoints.append((i,j))
    leftmostStartPoint = sorted(startPoints, key=lambda x: x[1])[0]
    return leftmostStartPoint # tuple (i,j)### Helpers

def arcOrder(arcPositions, startPoint):
    # Starting from {startPoint} calculate the closest pixel within {arcPositions}, add to new list to be returned, then remove it from arcPositions
    ordered = []
    copyArcPositions = arcPositions.copy()
    prevPoint = startPoint
    while len(copyArcPositions) != 0:
        minDist = None
        closestPoint = None      
        for curPoint in copyArcPositions:
            curDist = calcDist(curPoint,prevPoint)
            if (minDist is None) or (curDist<=minDist):
                minDist = curDist
                closestPoint = curPoint
        copyArcPositions.remove(closestPoint)
        ordered.append(closestPoint)
        prevPoint = closestPoint
    return ordered

def read_fits(arcPixelOrdered,waveband,galNum):
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband))
    image_data = hdul[0].data # 2-D Numpy Array
    hdul.close()
    fitsOrdered = [image_data[i,j] for i,j in arcPixelOrdered]
    return fitsOrdered



def MAIN(waveband1, waveband2, galNum, position, group):
    pixelLoc1, center, rows, cols = read_clusMask(waveband=waveband1,galNum=galNum,group=group)
    waveband1Arm      = get_singleArm(galaxyArms=pixelLoc1, position=position)
    pixelLoc2, center, rows, cols = read_clusMask(waveband=waveband2,galNum=galNum,group=group)
    waveband2Arm      = get_singleArm(galaxyArms=pixelLoc2, position=position)

    armPosition                          = mergeArms(arms=[waveband1Arm,waveband2Arm])
    minPhi, maxPhi, minRadius, maxRadius = arm_Span(armPosition=armPosition,center=center)
    arcsCircles_Positions                = arm_to_ArcsCircles(armPosition=armPosition, center=center, minPhi=minPhi, maxPhi=maxPhi, minRadius=minRadius, maxRadius=maxRadius)

### FOR ONLY THE INNER RADIUS
    radius, arc, circle = arcsCircles_Positions[int(len(arcsCircles_Positions)/2)]
    leftStartPoint  = arcStart(arcPositions=arc)
    arcPixelOrdered = arcOrder(arcPositions=arc, startPoint=leftStartPoint)
    fitsOrdered1    = read_fits(arcPixelOrdered=arcPixelOrdered,waveband=waveband1,galNum=galNum)
    fitsOrdered2    = read_fits(arcPixelOrdered=arcPixelOrdered,waveband=waveband2,galNum=galNum)
    subtracted      = [fitsOrdered1[index]-fitsOrdered2[index] for index in range(len(fitsOrdered1))]
    step_visualizer(rows, cols, waveband1,waveband2,galNum,position,pixelLoc1,pixelLoc2,waveband1Arm,waveband2Arm,armPosition,arc,leftStartPoint,arcPixelOrdered,circle,subtracted)
### FOR ALL POSSIBLE RADIUS
    # for radius, arc, circle in arcsCircles_Positions:
    #     leftStartPoint  = arcStart(arcPositions=arc)
    #     arcPixelOrdered = arcOrder(arcPositions=arc, startPoint=leftStartPoint)
    #     fitsOrdered1    = read_fits(arcPixelOrdered=arcPixelOrdered,waveband=waveband1,galNum=galNum)
    #     fitsOrdered2    = read_fits(arcPixelOrdered=arcPixelOrdered,waveband=waveband2,galNum=galNum)
    #     subtracted      = [fitsOrdered1[index]-fitsOrdered2[index] for index in range(len(fitsOrdered1))]
    #     step_visualizer(rows, cols, waveband1,waveband2,galNum,position,pixelLoc1,pixelLoc2,waveband1Arm,waveband2Arm,armPosition,arc,leftStartPoint,arcPixelOrdered,circle,subtracted)
    #     break
    
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
    plt.xlabel('Pixels from Leftmost Endpoint',size=15)
    plt.suptitle("STEPS TAKEN")
    plt.subplots_adjust(left=0.05,right=0.95,top=0.92,wspace=0.30,hspace=0.35)
    mng = plt.get_current_fig_manager()
    mng.window.state('zoomed')
    plt.show()



if __name__ == "__main__":
    MAIN(waveband1='g',waveband2='r',galNum='1237660635996291172', position=(110,110), group=True)