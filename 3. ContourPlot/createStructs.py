import math
import numpy as np


def get_largestArm(galaxyArms):
    """Returns SET of the largest arm in galaxy
            -most pixels, not length"""
    largestArmPosition = set()
    largestArmSize = 0
    for positions in galaxyArms.values():
        if len(positions)>=largestArmSize:
            largestArmSize = len(positions)
            largestArmPosition = positions
    return largestArmPosition



def get_PositionArm(galaxyArms,i,j):
    """ 
        Returns arm closest to the (i,j) position <-- (x,y) or (y,x)

        CHANGE THIS TO LOOK FOR CLOSEST MATCH NOT EXACT MATCH?

        MAKE SECOND FUNCTION TO GET SECOND LARGEST ARM?
    """
    armGood = set()
    for arm in galaxyArms.values():
        for iI,jJ in arm:
            if iI==i and jJ==j:
                armGood = arm
    return armGood



def unionClosestArm(waveband1LargestArm, waveband2AllArms):
    """ Returns the arm from waveband2 that has
        the highest overall overlap with waveband1s largest arm.
            -Overlap here is looking at the size of the intersection of the 
             2 arms over the total size of waveband1's arm
             (waveband1 && waveband2) / (waveband1)

        rtype: set()
    """
    overlapAmt = 0
    closestArm = None
    for armPos in waveband2AllArms.values():
        curOverlap = len(waveband1LargestArm.intersection(armPos)) / len(waveband1LargestArm)
        if (curOverlap>=overlapAmt):
            overlapAmt = curOverlap
            closestArm = armPos
    return waveband1LargestArm.union(closestArm) 



def arm_to_ArcsEllipse(majorAxis, minMaxRatio, axisRadians, armsPixels, center):
    arcsEllipse_Positions = []
    overallMinTheta = None
    overallMaxTheta = None
    for semiMajLen in range(5, int(majorAxis/2)):# Initial is 5 as any smaller just looks weird, probably not useful
        currentRadius = ellipseInfo(semiMajLen*2,[],set(),None,None)
        # Skip all small majorAxis's without little actual
        absoluteOverlapCount = 0 
        for curTheta in range(360):
            i,j = calcElpsPoint(semiMajLen, semiMajLen*minMaxRatio, axisRadians, curTheta, center)
            if ((i,j) in armsPixels):
                absoluteOverlapCount += 1
        if (absoluteOverlapCount >= 5):

            # Obtain startingTheta (eliminates 360-0 overlap)
            startingTheta=calcStartingTheta(semiMajLen, minMaxRatio, axisRadians, center, armsPixels)

            # Obtain min/max angle of the arc overlap (relative to startingTheta)
            for curTheta in range(startingTheta,startingTheta+360):
                i,j = calcElpsPoint(semiMajLen, semiMajLen*minMaxRatio, axisRadians, curTheta, center)
                if ((i,j) in armsPixels):
                    if (currentRadius.minTheta is None):
                        currentRadius.minTheta = curTheta
                    currentRadius.maxTheta = curTheta    
            if (overallMinTheta is None) or (currentRadius.minTheta < overallMinTheta):
                overallMinTheta = currentRadius.minTheta
            if (overallMaxTheta is None) or (currentRadius.maxTheta > overallMaxTheta):
                overallMaxTheta = currentRadius.maxTheta 
            
            # Stores the points for ARC and ELLIPSE
            for curTheta in range(startingTheta,startingTheta+360):
                i,j = calcElpsPoint(semiMajLen, semiMajLen*minMaxRatio, axisRadians, curTheta, center)
                currentRadius.ellipse.add( (i,j) )
                # If within the theta range for the ARC
                if (curTheta >= currentRadius.minTheta) and (curTheta <= currentRadius.maxTheta):
                    currentRadius.arc.append( (curTheta,i,j) )
            arcsEllipse_Positions.append(currentRadius)
    # Return a list of ellipseInfo objects
    return arcsEllipse_Positions, overallMinTheta, overallMaxTheta



def groupNeighborThetas(middle, neighbors):
    """
        Takes a single(middle) aep object and a list of neighboring aep objects.
        
        For each theta in the middle aep, get all same thetas in the neighboring aep's 
        and group together, returning a new list of thetas with repsective (i,j) lists
    """
    groupedThetas = []
    for theta, i , j in middle.arc:
        groupedThetas.append((theta,[(i,j)]))
        for aep_Obj in neighbors:
            for Itheta, Ii, Ij in aep_Obj.arc:
                if (theta == Itheta):
                    groupedThetas[-1][1].append((Ii,Ij))
                    break
                if (Itheta > theta):
                    break
    return groupedThetas



def calcFlux(i,j,fits1,fits2):
    """
        Given an ( i,j ) position, returns the average (float) of all 8 values surrounding it in a square
    """
    flux = 0
    if i==0:
        if j==0:
            flux += fits1[i,j]      - fits2[i,j]
            flux += fits1[i,j+1]    - fits2[i,j+1]

            flux += fits1[i+1,j]    - fits2[i+1,j]
            flux += fits1[i+1,j+1]  - fits2[i+1,j+1]
            avgFlux = flux/4
        elif j==len(fits1)-1:
            flux += fits1[i,j]      - fits2[i,j]
            flux += fits1[i,j-1]    - fits2[i,j-1]

            flux += fits1[i+1,j]    - fits2[i+1,j]
            flux += fits1[i+1,j-1]  - fits2[i+1,j-1]
            avgFlux = flux/4
        else:
            flux += fits1[i,j]      - fits2[i,j]
            flux += afits1[i,j+1]   - fits2[i,j+1]
            flux += fits1[i,j-1]    - fits2[i,j-1]

            flux += fits1[i+1,j]    - fits2[i+1,j]
            flux += fits1[i+1,j+1]  - fits2[i+1,j+1]
            flux += fits1[i+1,j-1]  - fits2[i+1,j-1]
            avgFlux = flux/6
    elif i==len(fits1)-1:
        if j==0:
            flux += fits1[i,j]      - fits2[i,j]
            flux += fits1[i,j+1]    - fits2[i,j+1]

            flux += fits1[i-1,j]    - fits2[i-1,j]
            flux += fits1[i-1,j+1]  - fits2[i-1,j+1]
            avgFlux = flux/4
        elif j==len(fits1)-1:
            flux += fits1[i,j]      - fits2[i,j]
            flux += fits1[i,j-1]    - fits2[i,j-1]

            flux += fits1[i-1,j]    - fits2[i-1,j]
            flux += fits1[i-1,j-1]  - fits2[i-1,j-1]
            avgFlux = flux/4
        else:
            flux += fits1[i,j]      - fits2[i,j]
            flux += fits1[i,j+1]    - fits2[i,j+1]
            flux += fits1[i,j-1]    - fits2[i,j-1]

            flux += fits1[i-1,j]    - fits2[i-1,j]
            flux += fits1[i-1,j+1]  - fits2[i-1,j+1]
            flux += fits1[i-1,j-1]  - fits2[i-1,j-1]
            avgFlux = flux/6
    elif j==0:
        if i==0:
            # Already accounted for
            pass
        elif i==len(fits1)-1:
            # Already accounted for
            pass
        else:
            flux += a[i,j]          - fits2[i,j]
            flux += a[i,j+1]        - fits2[i,j+1]

            flux += a[i+1,j]        - fits2[i+1,j]
            flux += a[i+1,j+1]      - fits2[i+1,j+1]

            flux += a[i-1,j]        - fits2[i-1,j]
            flux += a[i-1,j+1]      - fits2[i-1,j+1]
            avgFlux = flux/6
    elif j==len(fits1)-1:
        if i==0:
            # Already accounted for
            pass
        elif i==len(fits1)-1:
            # Already accounted for
            pass
        else:
            flux += fits1[i,j]      - fits2[i,j]
            flux += fits1[i,j-1]    - fits2[i,j-1]

            flux += fits1[i+1,j]    - fits2[i+1,j]
            flux += fits1[i+1,j-1]  - fits2[i+1,j-1]

            flux += fits1[i-1,j]    - fits2[i-1,j]
            flux += fits1[i-1,j-1]  - fits2[i-1,j-1]
            avgFlux = flux/6
    else:
        flux += fits1[i,j]          - fits2[i,j]
        flux += fits1[i,j+1]        - fits2[i,j+1]
        flux += fits1[i,j-1]        - fits2[i,j-1]

        flux += fits1[i+1,j]        - fits2[i+1,j]
        flux += fits1[i+1,j+1]      - fits2[i+1,j+1]
        flux += fits1[i+1,j-1]      - fits2[i+1,j-1]

        flux += fits1[i-1,j]        - fits2[i-1,j]
        flux += fits1[i-1,j+1]      - fits2[i-1,j+1]
        flux += fits1[i-1,j-1]      - fits2[i-1,j-1]
        avgFlux = flux/9
    return avgFlux



def getMinMaxFlux(groupedThetas, fits1, fits2):
    """
    Compute the min/max flux of this groupedTheta
    """
    curRadiusMinFlux = None
    curRadiusMaxFLux = None
    for theta,pixelList in groupedThetas:
        for i,j in pixelList:
            avgFlux = calcFlux(i,j,fits1,fits2)
            if (curRadiusMinFlux is None) or (avgFlux < curRadiusMinFlux):
                curRadiusMinFlux = avgFlux
            if (curRadiusMaxFLux is None) or (avgFlux > curRadiusMaxFLux):
                curRadiusMaxFLux = avgFlux
    return curRadiusMinFlux, curRadiusMaxFLux



def calcOverallMinMaxFlux(arcsEllipse_Positions,merge,fits1,fits2):
    """
    Find the overall min/max flux for every aepObj's grouped thetas
    """
    overallMinFlux = None
    overallMaxFlux = None
    for i in range(merge,len(arcsEllipse_Positions)-merge): # SHOULD CALCULATE MIN/MAX FLUX FROM CURRENT AEPOBJ DIRECTLY, dont need neighbors-extra work
        current_aepObj = arcsEllipse_Positions[i]
        neighbor_aepObjs = [arcsEllipse_Positions[j] for j in np.arange(i-merge,i+merge+1) if (i!=j)]
        groupedThetas = groupNeighborThetas(middle=current_aepObj, neighbors=neighbor_aepObjs)
        curRadiusMinFlux, curRadiusMaxFlux = getMinMaxFlux(groupedThetas=groupedThetas, fits1=fits1, fits2=fits2)
        if (overallMinFlux is None) or (curRadiusMinFlux < overallMinFlux):
            overallMinFlux = curRadiusMinFlux
        if (overallMaxFlux is None) or (curRadiusMaxFlux > overallMaxFlux):
            overallMaxFlux = curRadiusMaxFlux
    return overallMinFlux, overallMaxFlux



def calcElpsPoint(a, b, axisRadians, curTheta, center):
    # https://math.stackexchange.com/questions/315386/ellipse-in-polar-coordinates
    radians = curTheta/180*math.pi

    bottomB = (b*math.cos(radians -axisRadians -math.pi/2))**2     # PI/2 must be here to correctly rotate (BASICALLY DOING adjust Theta?)
    bottomA = (a*math.sin(radians -axisRadians -math.pi/2))**2     # diskMajAxAngleRadians rotates COUNTER-CLOCKWISE(both a>b & a<b)
    top = (a*b)
    bottom = math.sqrt(bottomB + bottomA)
    r = top/bottom

    i = round(r*math.cos(radians) + center[0])
    j = round(r*math.sin(radians) + center[1])
    return i,j # Determine whether this is rows,col or col,row



def calcStartingTheta(semiMajLen, minMaxRatio, axisRadians, center, armsPixels):
    for curTheta in range(40,360,40):
        emptyGap = True
        for innerTheta in range(curTheta-40,curTheta+40):
            i,j = calcElpsPoint(semiMajLen, semiMajLen*minMaxRatio, axisRadians, innerTheta, center)
            if ((i,j) in armsPixels):
                emptyGap = False      
                break
        if (emptyGap):
            startingTheta=curTheta
            break
    return startingTheta



def armOutline(armsPixels):
    """
        Given pixels of an arm, gets only pixels that have >= threshold missing in
        a square around it, returing an outline of the arm

        Increase threshold to get thicker outline, Smaller threshold has some points missing from outline
    """
    threshold = 3
    outline=set()
    for i,j in armsPixels:
        missing = 0
        if ((i,j+1) not in armsPixels):
            missing += 1
        if ((i,j-1) not in armsPixels):
            missing += 1
        if ((i+1,j) not in armsPixels):
            missing += 1
        if ((i-1,j) not in armsPixels):
            missing += 1
        if ((i+1,j+1) not in armsPixels):
            missing += 1
        if ((i+1,j-1) not in armsPixels):
            missing += 1
        if ((i-1,j+1) not in armsPixels):
            missing += 1
        if ((i-1,j-1) not in armsPixels):
            missing += 1
        if (missing >= threshold):
            outline.add( (i,j) )
    return outline



def updateThetaStarts(overallMinTheta, overallMaxTheta, arcsEllipse_Positions):
    """
        Change all ellipseinfo objects thetas to be relative to start, 0 from FRONT

        Returns new overall min/max thetas, sWise bool
    """
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
    return newOverallMinTheta, newOverallMaxTheta, sWise



class ellipseInfo:
    """
        majorAxisLen  
        arc  
        ellipse  
        minTheta  
        maxTheta  
    """
    def __init__(self,majorAxisLen,arc,ellipse,minTheta,maxTheta):
        self.majorAxisLen   = majorAxisLen
        self.arc            = arc
        self.ellipse        = ellipse
        self.minTheta       = minTheta
        self.maxTheta       = maxTheta

