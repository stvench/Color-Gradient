from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math
from statistics import mean


def read_clusMask(waveband, galNum, group):
    """
    Returns a dict containing {color:set of position pairs}
        - color:    color
        - position: (i,j) pixel locations
    """
    # Get Image of Galaxy Arms
    imgClusMask = cv2.imread("{}/{}/{}-K_clusMask-reprojected.png".format(galNum,waveband,galNum))
    dimensions = imgClusMask.shape
    rows = dimensions[0]
    cols = dimensions[1]
    # Store all Pixels of arms
    allPixelLoc = dict()
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
            if len(locs) <= 10:
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
    return allPixelLoc, (math.ceil(rows/2),math.ceil(cols/2)), rows, cols # {color : {(i,j) , (i,j)} } , (i,j) }

def get_closestArm(galaxyArms, position):
    """Returns closest arm to the inputted position (row,col) | (y,x) 
        -Ignores "small" pixels merged into similar color"""
    closestArm, minDist = None, None
    for arm in galaxyArms.values():
        if (len(arm) > 10):
            for point in arm:
                curDist = calcDist(point,position)
                if (minDist is None) or (curDist<=minDist):
                    minDist = curDist
                    closestArm = arm
    return closestArm

def get_largestArm(galaxyArms):
    """Returns the largest arm in galaxy | mean postion
            -most pixels, not length"""
    largestArmPosition = set()
    largestArmSize = 0
    for positions in galaxyArms.values():
        if len(positions)>=largestArmSize:
            largestArmSize = len(positions)
            largestArmPosition = positions
    return largestArmPosition




def mergeArms(waveband1Arm,waveband2Arm):
    ### Check for overlap (needs at least 50% in one)
    overlapThreshold = 0.5
    overlap = waveband1Arm.intersection(waveband2Arm)
    arm1Overlap = len(overlap)/len(waveband1Arm)
    arm2Overlap = len(overlap)/len(waveband2Arm)
    if (arm1Overlap > overlapThreshold) or (arm2Overlap > overlapThreshold):
        return waveband1Arm.union(waveband2Arm)
    else:
        print("{} or {} <= {} Threshold".format(arm1Overlap,arm2Overlap,overlapThreshold))
        raise TypeError
    
def arm_to_ArcsCircles(armPosition,center):
    allDists = [math.dist(point,center) for point in armPosition]
    minRadius = math.ceil(min(allDists))
    maxRadius = math.floor(max(allDists))
    arcsCircles_Positions = []
    for radius in range(minRadius,maxRadius+1):
        arcsCircles_Positions.append( radiusInfo(radius,[],set(),None,None) )
        minPhi = None
        maxPhi = None
        # Obtain startingPhi (eliminates 360-0 overlap)
        startingPhi=None
        for phi in range(40,360,40):
            emptyGap = True
            for curPhi in range(phi-40,phi+40):
                radian = curPhi/180*math.pi
                i = round(radius*math.cos(radian)+center[0])
                j = round(radius*math.sin(radian)+center[1])
                curDist = calcDist(center,(i,j))
                if (curDist >= radius-0.5) and (curDist <= radius+0.5) and ((i,j) in armPosition):
                    emptyGap = False      
                    break
            if (emptyGap):
                startingPhi=phi
                break
        # Obtain min/max angle of the arc overlap (relative to startingPhi)
        for phi in range(startingPhi,startingPhi+360):
            radian = phi/180*math.pi
            i = round(radius*math.cos(radian)+center[0])
            j = round(radius*math.sin(radian)+center[1])
            curDist = calcDist(center,(i,j))
            if (curDist >= radius-0.5) and (curDist <= radius+0.5) and ((i,j) in armPosition):
                if (minPhi is None):
                    minPhi = phi
                maxPhi = phi      
        arcsCircles_Positions[-1].minPhi = adjustPhi(minPhi)
        arcsCircles_Positions[-1].maxPhi = adjustPhi(maxPhi)+360 if adjustPhi(maxPhi)<adjustPhi(minPhi) else adjustPhi(maxPhi)
        # Stores the points for ARC and CIRCLE
        for phi in range(startingPhi,startingPhi+360):
            radian = phi/180*math.pi
            i = round(radius*math.cos(radian)+center[0])
            j = round(radius*math.sin(radian)+center[1])
            curDist = calcDist(center,(i,j))
            # If the current point is within 0.5 of the actual radius distance for the CIRCLE
            if (curDist >= radius-0.5) and (curDist <= radius+0.5):
                arcsCircles_Positions[-1].circle.add((i,j))
                # If within the phi range for the ARC
                if (phi >= minPhi) and (phi <= maxPhi):
                    arcsCircles_Positions[-1].arc.append((phi,i,j))
    # Returns a list of radiusInfo objects
    return arcsCircles_Positions




def simPhis(middle, inner, outer):
    phiList = []
    for phi, i , j in middle.arc:
        phiList.append((phi,[(i,j)]))
        for Iphi,Ii,Ij in inner.arc:
            if (phi == Iphi):
                phiList[-1][1].append((Ii,Ij))
                break
        for Ophi,Oi,Oj in outer.arc:
            if (phi == Ophi):
                phiList[-1][1].append((Oi,Oj))
                break
    ### REMOVES MOST DUPLICATE PHIs
    uniquePhiList = []
    prev = None
    same = []
    for ele in phiList:
        if ele[1] == prev:
            same.append(ele)
        else:
            if len(same)>0:
                uniquePhiList.append(same[int(len(same)/2)])
            same = [ele]
            prev = ele[1]
    if len(same)>0:
        uniquePhiList.append(same[int(len(same)/2)])

    return uniquePhiList

def unmergedFits(arcList,waveband1,waveband2,galNum):
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband1))
    image_data1 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband2))
    image_data2 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    return [(phi,image_data1[i,j]-image_data2[i,j]) for phi,i,j in (arcList)]

def mergedFits(l, waveband1,waveband2,galNum):
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband1))
    image_data1 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    hdul = fits.open("{}/{}_{}.fits".format(galNum,galNum,waveband2))
    image_data2 = hdul[0].data # 2-D Numpy Array
    hdul.close()
    return [(phi,mean([image_data1[i,j]-image_data2[i,j] for i,j in x])) for phi,x in (l)]





def manualTest(waveband1, waveband2, galNum, position, group):
    pixelLoc1, center, rows, cols = read_clusMask(waveband=waveband1,galNum=galNum,group=group)
    waveband1Arm      = get_closestArm(galaxyArms=pixelLoc1, position=position)
    
    pixelLoc2, center, rows, cols = read_clusMask(waveband=waveband2,galNum=galNum,group=group)
    waveband2Arm      = get_closestArm(galaxyArms=pixelLoc2, position=position)
    
    
    armPosition                          = mergeArms(arms=[waveband1Arm,waveband2Arm])
    arcsCircles_Positions                = arm_to_ArcsCircles(armPosition=armPosition, center=center)
    for i in range(1,len(arcsCircles_Positions)-1):
        if (len(arcsCircles_Positions[i].arc) > 30): #and (len(arcsCircles_Positions[i].arc) < 270)
            curRadiusInfo = arcsCircles_Positions[i]
            phiList = simPhis(curRadiusInfo,arcsCircles_Positions[i-1],arcsCircles_Positions[i+1])
            merged = mergedFits(phiList, waveband1,waveband2,galNum)
            mergedPhi = []
            mergedDifs = []
            for phi, dif in merged: # adjust phi first, then consider overlap
                mergedPhi.append(adjustPhi(phi)+360 if adjustPhi(phi)<curRadiusInfo.minPhi else adjustPhi(phi))
                mergedDifs.append(dif)
            unmerged = unmergedFits(curRadiusInfo.arc,waveband1,waveband2,galNum)
            unmergedPhi = []
            unmergedDifs = []
            for phi, dif in unmerged:
                unmergedPhi.append(adjustPhi(phi)+360 if adjustPhi(phi)<curRadiusInfo.minPhi else adjustPhi(phi))
                unmergedDifs.append(dif)
            mng = plt.get_current_fig_manager()
            mng.window.state('zoomed')

            plt.subplot(121)
            a = np.zeros((rows,cols))
            for ii, ij in armPosition:
                a[ii,ij]=2
            for ii, ij in curRadiusInfo.circle:
                a[ii,ij]=1
            for phi,ii, ij in curRadiusInfo.arc:
                a[ii,ij]=3
            plt.imshow(a) # origin='lower' for "normal" X/Y axis position
            plt.text(rows*0.6,cols*0.85,"Min θ: {}".format(curRadiusInfo.minPhi),color="white")
            plt.text(rows*0.6,cols*0.95,"Max θ: {}".format(curRadiusInfo.maxPhi),color="white")
            plt.title('Combination')

            plt.subplot(222)
            plt.plot(unmergedPhi,unmergedDifs)
            plt.title("Unmerged Radius {}, {}-{}".format(curRadiusInfo.radius,waveband1,waveband2))
            plt.subplot(224)

            plt.plot(mergedPhi,mergedDifs)
            plt.xlabel("θ (degrees)",size=15)
            plt.ylabel("Flux (nanomaggies)",size=15)
            plt.title("Merged Radius ({},{},{}), {}-{}".format(curRadiusInfo.radius-1,curRadiusInfo.radius,curRadiusInfo.radius+1,waveband1,waveband2))
            plt.subplots_adjust(hspace=0.3,wspace=0.4)
            plt.suptitle("Radius: {}".format(curRadiusInfo.radius), size=20)
            # plt.savefig("{}/{}_{}_{}-{}.pdf".format(galNum,position,curRadiusInfo.radius,waveband1,waveband2))
            plt.show()
            plt.close()
            # break

def automateTest(waveband1, waveband2, galNum):
    """
    WAVEBANDS ["u","g","r","i","z"]     
    Parameters:
        waveband1 (str) : From WAVEBANDS
        waveband2 (str) : From WAVEBANDS
        galNum    (int) :
    Returns:
        None (but plots graph)
    """
    # try:
    group = True
    pixelLoc1, center, rows, cols = read_clusMask(waveband=waveband1,galNum=galNum,group=group)
    waveband1Arm      = get_largestArm(galaxyArms=pixelLoc1)
    pixelLoc2, center, rows, cols = read_clusMask(waveband=waveband2,galNum=galNum,group=group)
    waveband2Arm      = get_largestArm(galaxyArms=pixelLoc2)
    
    armPosition                          = mergeArms(waveband1Arm=waveband1Arm,waveband2Arm=waveband2Arm)
    arcsCircles_Positions                = arm_to_ArcsCircles(armPosition=armPosition, center=center)
    for i in range(1,len(arcsCircles_Positions)-1):
        if (len(arcsCircles_Positions[i].arc) > 30): #and (len(arcsCircles_Positions[i].arc) < 270)
            curRadiusInfo = arcsCircles_Positions[i]
            phiList = simPhis(curRadiusInfo,arcsCircles_Positions[i-1],arcsCircles_Positions[i+1])
            merged = mergedFits(phiList, waveband1,waveband2,galNum)
            mergedPhi = []
            mergedDifs = []
            for phi, dif in merged: # adjust phi first, then consider overlap
                mergedPhi.append(adjustPhi(phi)+360 if adjustPhi(phi)<curRadiusInfo.minPhi else adjustPhi(phi))
                mergedDifs.append(dif)
            unmerged = unmergedFits(curRadiusInfo.arc,waveband1,waveband2,galNum)
            unmergedPhi = []
            unmergedDifs = []
            for phi, dif in unmerged:
                unmergedPhi.append(adjustPhi(phi)+360 if adjustPhi(phi)<curRadiusInfo.minPhi else adjustPhi(phi))
                unmergedDifs.append(dif)
            mng = plt.get_current_fig_manager()
            mng.window.state('zoomed')

            plt.subplot(121)
            a = np.zeros((rows,cols))
            for ii, ij in armPosition:
                a[ii,ij]=2
            for ii, ij in curRadiusInfo.circle:
                a[ii,ij]=1
            for phi,ii, ij in curRadiusInfo.arc:
                a[ii,ij]=3
            plt.imshow(a) # origin='lower' for "normal" X/Y axis position
            plt.text(rows*0.6,cols*0.85,"Min θ: {}".format(curRadiusInfo.minPhi),color="white")
            plt.text(rows*0.6,cols*0.95,"Max θ: {}".format(curRadiusInfo.maxPhi),color="white")
            plt.title('Combination')

            plt.subplot(222)
            plt.plot(unmergedPhi,unmergedDifs)
            plt.title("Unmerged Radius {}, {}-{}".format(curRadiusInfo.radius,waveband1,waveband2))
            plt.subplot(224)

            plt.plot(mergedPhi,mergedDifs)
            plt.xlabel("θ (degrees)",size=15)
            plt.ylabel("Flux (nanomaggies)",size=15)
            plt.title("Merged Radius ({},{},{}), {}-{}".format(curRadiusInfo.radius-1,curRadiusInfo.radius,curRadiusInfo.radius+1,waveband1,waveband2))
            plt.subplots_adjust(hspace=0.3,wspace=0.4)
            plt.suptitle("Radius: {}".format(curRadiusInfo.radius), size=20)
            # plt.savefig("{}/{}_{}_{}-{}.pdf".format(galNum,position,curRadiusInfo.radius,waveband1,waveband2))
            plt.show()
            plt.close()
            break
    # except:
    #     print("{} {}-{} Failed".format(galNum,waveband1,waveband2))
    #     print("-------------------------------------------")



class radiusInfo:
    def __init__(self,radius,arc,circle,minPhi,maxPhi):
        self.radius = radius
        self.arc    = arc
        self.circle = circle
        self.minPhi = minPhi
        self.maxPhi = maxPhi

### Helper Functions
def adjustPhi(phi):
    ''' Adjust Phi 90 degrees counterclockwise (Easier than changing whole code) '''
    phi-=90
    if (phi<0):
        phi = 360+phi
    return phi

def calcDist(point,center):
    return math.sqrt((point[0]-center[0])**2+(point[1]-center[1])**2)

if __name__ == "__main__":
    """
                        --------------REMINDER--------------
    When setting position, uses (ROW,COL) which is equal to (Y,X) in the plotted axis
    """
    ############## MANUALLY CHECK POSITIONING
    ### 1237648702986125622
    #
    # manualTest(waveband1='g',waveband2='i',galNum='1237648702986125622', position=(80,70), group=True)      #G Dark Green
    # manualTest(waveband1='g',waveband2='i',galNum='1237648702986125622', position=(100,100), group=True)    #G Light Green


    ### 1237648704586514654
    #
    # manualTest(waveband1='g',waveband2='i',galNum='1237648704586514654', position=(90,75), group=True)      #G Top Red
    # manualTest(waveband1='g',waveband2='i',galNum='1237648704586514654', position=(100,120), group=True)    #G Bottom Orange


    ### 1237648705658486867
    #
    # manualTest(waveband1='g',waveband2='i',galNum='1237648705658486867', position=(45,45), group=True)      #G Red


    ### 1237660635996291172
    #
    # manualTest(waveband1='g',waveband2='i',galNum='1237660635996291172', position=(110,110), group=True)    #G GREEN
    # manualTest(waveband1='g',waveband2='i',galNum='1237660635996291172', position=(120,140), group=True)    #G PURPLE
    # manualTest(waveband1='g',waveband2='i',galNum='1237660635996291172', position=(140,120), group=True)    #G ORANGE



    ############## AUTOMATE POSITIONING
    # automateTest(waveband1='g',waveband2='i',galNum='1237660635996291172')
    # automateTest(waveband1='g',waveband2='i',galNum='1237648705658486867')
    # automateTest(waveband1='g',waveband2='i',galNum='1237648704586514654')
    # automateTest(waveband1='g',waveband2='i',galNum='1237648702986125622')
    print("OLD VERSION - NOT COMPATIBLE WITH OPENLAB")


    # imgClusMask = cv2.imread("1237648705658486867/i/1237648705658486867-K_clusMask-reprojected.png")
    # plt.imshow(imgClusMask)
    # plt.show()