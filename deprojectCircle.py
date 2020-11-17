import matplotlib.pyplot as plt
import numpy as np
import math
import cv2





def structureInfo(galNum): # P_CS >=
    filePath = f"C:/Users/sc123/Desktop/Galaxys/r/{galNum}/{galNum}.tsv"
    with open(filePath, 'r') as data:
        foundHeaders = 0
        ##### FIND info
        for dataIndex, name in enumerate(data.readline().split('\t')):
            if name == 'diskAxisRatio':
                diskAxisRatioIndex = dataIndex
                foundHeaders += 1
            if name == 'diskMinAxsLen':
                diskMinAxsLenIndex = dataIndex
                foundHeaders += 1
            if name == 'diskMajAxsLen':
                diskMajAxsLenIndex = dataIndex
                foundHeaders += 1
            if name == 'diskMajAxsAngleRadians':
                diskMajAxsAngleRadiansIndex = dataIndex
                foundHeaders += 1
            if name == 'inputCenterR':
                inputCenterRIndex = dataIndex
                foundHeaders += 1
            if name == 'inputCenterC':
                inputCenterCIndex = dataIndex
                foundHeaders += 1
            if name == 'iptSz':
                iptSzIndex = dataIndex
                foundHeaders += 1
        if (foundHeaders != 7):
            print("Missing Header Index")
            return
        ##### GET info
        dataInfo = data.readline().split('\t')
        diskAxisRatio           = float(dataInfo[diskAxisRatioIndex])           #   diskAxisRatio = diskMinAxsLen / diskMajAxsLen
        diskMinAxsLen           = float(dataInfo[diskMinAxsLenIndex])
        diskMajAxsLen           = float(dataInfo[diskMajAxsLenIndex])
        diskMajAxsAngleRadians  = float(dataInfo[diskMajAxsAngleRadiansIndex])
        inputCenterR            = float(dataInfo[inputCenterRIndex])
        inputCenterC            = float(dataInfo[inputCenterCIndex])
        iptSz                   = dataInfo[iptSzIndex]
        #####
        print("MINOR / MAJOR:",diskAxisRatio)
        print("MINOR:",diskMinAxsLen)           ### MAX?
        print("MAJOR:",diskMajAxsLen)           ### MAX?
        print("ANGLE:",diskMajAxsAngleRadians)
        print("SIZE:", iptSz)
        print("CENTER:", inputCenterC, inputCenterR)  
        # https://math.stackexchange.com/questions/315386/ellipse-in-polar-coordinates
        ### A is horizontal x (diskMajAxsLen)
        ### B is vertical y   (diskMinAxsLen)
        for a in range(1,int(diskMajAxsLen),1):
            b = a*diskAxisRatio
            ellipse = []
            for i in range(90,450):                                                  # Changing to (90,450) moves start THETA to intuitive point
                radians = i/180*math.pi

                bottomB = (b*math.cos(radians -diskMajAxsAngleRadians -math.pi/2))**2     # PI/2 must be here to correctly rotate
                bottomA = (a*math.sin(radians -diskMajAxsAngleRadians -math.pi/2))**2     # diskMajAxAngleRadians rotates COUNTER-CLOCKWISE(both a>b & a<b)
                top = (a*b)
                bottom = math.sqrt(bottomB + bottomA)
                r = top/bottom

                x = r*math.cos(radians) + 70
                y = r*math.sin(radians) + 70
                ellipse.append((round(x),round(y)))





if __name__ == "__main__":
    

    structureInfo("1237648703510478979")
