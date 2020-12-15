import math



def get_largestArm(galaxyArms):
    """Returns the largest arm in galaxy
            -most pixels, not length"""
    largestArmPosition = set()
    largestArmSize = 0
    for positions in galaxyArms.values():
        if len(positions)>=largestArmSize:
            largestArmSize = len(positions)
            largestArmPosition = positions
    return largestArmPosition



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



def remvSimThetas(middle, neighbors):
    ### Change function name, doesnt remove similar thetas across neighbors, just groups thetas together
    thetaList = []
    for theta, i , j in middle.arc:
        thetaList.append((theta,[(i,j)]))
        for aep_Obj in neighbors:
            for Itheta, Ii, Ij in aep_Obj.arc:
                if (theta == Itheta):
                    thetaList[-1][1].append((Ii,Ij))
                    break
                if (Itheta > theta):
                    break
    return thetaList



def calcFlux(i,j,fits1,fits2):
    # Given an (i,j) position, returns the average of all 8 values surronding it, unless on border
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
    return i,j



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



class ellipseInfo:
    def __init__(self,majorAxisLen,arc,ellipse,minTheta,maxTheta):
        self.majorAxisLen   = majorAxisLen
        self.arc            = arc
        self.ellipse        = ellipse
        self.minTheta       = minTheta
        self.maxTheta       = maxTheta

