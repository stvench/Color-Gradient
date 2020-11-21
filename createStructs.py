from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math
from statistics import mean



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


   
def arm_to_ArcsEllipse(majorAxis, minMaxRatio, axisRadians, armsPixels, center):
    arcsEllipse_Positions = []
    for a in range(5, int(majorAxis/2)):
        arcsEllipse_Positions.append( ellipseInfo(a*2,[],set(),None,None) )
        minTheta = None
        maxTheta = None
        # Skip all small majorAxis's without little actual
        absoluteOverlapCount = 0 
        for curTheta in range(40,360):
            i,j = calcElpsPoint(a, a*minMaxRatio, axisRadians, curTheta, center)
            if ((i,j) in armsPixels):
                absoluteOverlapCount += 1
        if (absoluteOverlapCount >= 5):
            # Obtain startingTheta (eliminates 360-0 overlap)
            startingTheta=None
            for curTheta in range(40,360,40):
                emptyGap = True
                for innerTheta in range(curTheta-40,curTheta+40):
                    i,j = calcElpsPoint(a, a*minMaxRatio, axisRadians, innerTheta, center)
                    if ((i,j) in armsPixels):
                        emptyGap = False      
                        break
                if (emptyGap):
                    startingTheta=curTheta
                    break
            # Obtain min/max angle of the arc overlap (relative to startingTheta)
            for curTheta in range(startingTheta,startingTheta+360):
                i,j = calcElpsPoint(a, a*minMaxRatio, axisRadians, curTheta, center)
                if ((i,j) in armsPixels):
                    if (minTheta is None):
                        minTheta = curTheta
                    maxTheta = curTheta      
            arcsEllipse_Positions[-1].minTheta = adjustTheta(minTheta)
            arcsEllipse_Positions[-1].maxTheta = adjustTheta(maxTheta)+360 if adjustTheta(maxTheta)<adjustTheta(minTheta) else adjustTheta(maxTheta)
            # Stores the points for ARC and ELLIPSE
            for curTheta in range(startingTheta,startingTheta+360):
                i,j = calcElpsPoint(a, a*minMaxRatio, axisRadians, curTheta, center)
                arcsEllipse_Positions[-1].ellipse.add((i,j))
                # If within the theta range for the ARC
                if (curTheta >= minTheta) and (curTheta <= maxTheta):
                    arcsEllipse_Positions[-1].arc.append((curTheta,i,j))
    # Return a list of ellipseInfo objects
    return arcsEllipse_Positions



def remvSimThetas(middle, inner, outer):
    thetaList = []
    for theta, i , j in middle.arc:
        thetaList.append((theta,[(i,j)]))
        for Itheta,Ii,Ij in inner.arc:
            if (theta == Itheta):
                thetaList[-1][1].append((Ii,Ij))
                break
        for Otheta,Oi,Oj in outer.arc:
            if (theta == Otheta):
                thetaList[-1][1].append((Oi,Oj))
                break
    ### REMOVES MOST DUPLICATE PHIs (ones that point to the same i,j position)
    uniqueThetaList = []
    prev = None
    same = []
    for ele in thetaList:
        if ele[1] == prev:
            same.append(ele)
        else:
            if len(same)>0:
                uniqueThetaList.append(same[int(len(same)/2)])
            same = [ele]
            prev = ele[1]
    if len(same)>0:
        uniqueThetaList.append(same[int(len(same)/2)])
    return uniqueThetaList



def calcElpsPoint(a, b, axisRadians, curTheta, center):
    # https://math.stackexchange.com/questions/315386/ellipse-in-polar-coordinates
    radians = curTheta/180*math.pi

    bottomB = (b*math.cos(radians -axisRadians -math.pi/2))**2     # PI/2 must be here to correctly rotate
    bottomA = (a*math.sin(radians -axisRadians -math.pi/2))**2     # diskMajAxAngleRadians rotates COUNTER-CLOCKWISE(both a>b & a<b)
    top = (a*b)
    bottom = math.sqrt(bottomB + bottomA)
    r = top/bottom

    i = round(r*math.cos(radians) + center[0])
    j = round(r*math.sin(radians) + center[1])
    return i,j



class ellipseInfo:
    def __init__(self,majorAxisLen,arc,ellipse,minTheta,maxTheta):
        self.majorAxisLen   = majorAxisLen
        self.arc            = arc
        self.ellipse        = ellipse
        self.minTheta       = minTheta
        self.maxTheta       = maxTheta


def adjustTheta(theta):
    ''' Adjust Phi 90 degrees counterclockwise '''
    theta-=90
    if (theta<0):
        theta = 360+theta
    return theta