from astropy.io import fits
import numpy as np
import cv2
from statistics import mean



def read_arcsTSV(galNum, onOpenlabs):
    if (onOpenlabs):
        FILEPATH = f"/extra/wayne1/research/drdavis/SDSS/SpArcFiRe/2016-09/r/{galNum[-3:]}/{galNum}/{galNum}.tsv"
    else:
        FILEPATH = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/3. ContourPlot/data/Sample_Galaxies/{galNum}/r/{galNum}.tsv"

    with open(FILEPATH, 'r') as data:
        foundHeaders = 0
        ##### Get the Header, indexes of info needed
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
        minMaxRatio             = float(dataInfo[diskAxisRatioIndex])           #   diskAxisRatio = diskMinAxsLen / diskMajAxsLen
        minorAxis               = float(dataInfo[diskMinAxsLenIndex])
        majorAxis               = float(dataInfo[diskMajAxsLenIndex])
        axisRadians             = float(dataInfo[diskMajAxsAngleRadiansIndex])
        inputCenterR            = float(dataInfo[inputCenterRIndex])
        inputCenterC            = float(dataInfo[inputCenterCIndex])
        iptSz                   = dataInfo[iptSzIndex]
        #####
        # print("MINOR / MAJOR:",minMaxRatio)
        # print("MINOR:",minorAxis)           ### MAX?
        # print("MAJOR:",majorAxis)           ### MAX?
        # print("ANGLE:",axisRadians)
        # print("SIZE:", iptSz)
        # print("CENTER:", inputCenterC, inputCenterR)  
        return minMaxRatio, minorAxis, majorAxis, axisRadians, inputCenterR, inputCenterC


def read_imageAPng(waveband,galNum,onOpenlabs):
    """
    Returns the image as an np array, for use as referene
    """
    if (onOpenlabs):
        FILEPATH = f"/extra/wayne1/research/drdavis/SDSS/SpArcFiRe/2016-09/{waveband}/{galNum[-3:]}/{galNum}/{galNum}-A_input.png"
    else:
        FILEPATH = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/3. ContourPlot/data/Sample_Galaxies/{galNum}/{waveband}/{galNum}-A_input.png"
    imgPng = cv2.imread(FILEPATH)
    return imgPng


def read_clusMask(waveband, galNum, group, onOpenlabs):
    """
    Returns a dict containing {color:set of position pairs}
        - color:    color
        - position: (i,j) pixel locations
    """
    if (onOpenlabs):
        FILEPATH = f"/extra/wayne1/research/drdavis/SDSS/SpArcFiRe/2016-09/{waveband}/{galNum[-3:]}/{galNum}/{galNum}-K_clusMask-reprojected.png"
    else:
        FILEPATH = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/3. ContourPlot/data/Sample_Galaxies/{galNum}/{waveband}/{galNum}-K_clusMask-reprojected.png"

    imgClusMask = cv2.imread(FILEPATH)
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
                    # Ignore border pixels (QUIZ FIX BUT NOT THE BEST: GETS RID OF IndexErrors)
                    if (i != 0) and (i != rows-1) and ( j!= 0) and (i != cols-1):
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
    return allPixelLoc, rows, cols # {color : {(i,j) , (i,j)} } , (i,j) }


def readFits(waveband,galNum,onOpenlabs):
    if (onOpenlabs):
        FILEPATH = f"/extra/wayne1/research/drdavis/SDSS/FITS/color/{waveband}/{galNum[-3:]}/{galNum}.fits.gz"
    else:
        FILEPATH = f"C:/Users/sc123/Desktop/gal/5. NonCircular FITS/3. ContourPlot/data/Sample_Galaxies/{galNum}/{galNum}_{waveband}.fits"
    hdul = fits.open(FILEPATH)
    fitsData = hdul[0].data # 2-D Numpy Array
    hdul.close()
    return fitsData