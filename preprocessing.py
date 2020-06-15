from accessData import getHDFData
from voronoi import voronoi
import inSituCluster as clusterLib
import matplotlib.pyplot as plt
from IPython import embed
import numpy as np
#import cv2
import os
import scipy.spatial as kd

kdt = None
arr = None
gradx = None
grady = None
grad = None
imgCntr = 0
yPts = None
xPts = None


def createHist(inputData, lowerBound=None, upperBound=None, ignore=None):
    """ Create Histogram with matplotlib
        Parameters:
            data - The input data to be processed
            lowerBound - The smallest value to consider for the histogram
            upperBound - the largest value to consider for the histogram
            ignore - A list of values to ignore from the histogram
    """

    data = inputData.copy()
    if ignore is not None:
        for val in ignore:  # Consecutively remove the values to be ignored
            data = data[np.where(data != val)]
    if lowerBound is None:
        lowerBound = data.min()
    if upperBound is None:
        upperBound = data.max()
    sampleData = data[np.where(data < upperBound)]
    sampleData = sampleData[np.where(sampleData >= lowerBound)]
    sampleSize = sampleData.size
# Scott's (1979) bin width. Based on the standard deviation
# and the data size
    n = np.power(sampleSize, 1.0 / 3)
    std = np.std(sampleData)
    bin_width = (3.49 * std) / n
    histSize = int((upperBound - lowerBound) / bin_width)
    ranges = (lowerBound, upperBound)
    hist = plt.hist(sampleData, histSize, ranges)
    plt.axes().set_aspect('equal')
    plt.show()

    return hist


def quantize(data, cutoff, numRegions=8, tolerance=0.0001):

    newData = data[::4, ::4]  # Downsample data to speed up by a factor of 16
    # Next remove clouds poorly represented data
    cleanData = newData[np.logical_and(newData > 0, newData < cutoff)]

    tess = voronoi(cleanData, numRegions, tolerance, data)

    thresh = []
    for val in tess.values():
        thresh.append(np.min(val))
    thresh = np.sort(thresh)

    (rows, cols) = data.shape
    qData = np.ndarray((rows, cols), dtype='uint8')
    qData.fill(255)

    for val, t in enumerate(thresh):
        qData[data > t] = val + 1

    qData[data == 0] = 0  # Put back the land values
    qData[data == -1] = 255  # Put back the cloud values

    return qData


def colorMap(data):

    colors = [(123, 0, 255), (0, 0, 255), (0, 255, 255), (0, 255, 0),
              (123, 255, 0), (255, 255, 0), (255, 123, 0), (255, 0, 0)]

    (rows, cols) = data.shape
    color = np.ndarray((rows, cols, 3), dtype='uint8')
    color.fill(255)

    # Fill land with black
    color[data == 0] = (0, 0, 0)

    for i, c in enumerate(colors):
        color[data == i + 1] = c
    #plt.imshow(color)
    #plt.show()

    return color


def writeClouds(data):

    # Only clouds should be left as -1
    locs = np.where(data == -1)
    numPts = locs[0].size
    with open("data/clouds.dat", 'w') as fd:
        for i in xrange(numPts):
            if i % 100000 == 0:
                print i
            x = locs[0][i]
            y = locs[1][i]
            fd.write("%d %d\n" % (x, y))


def findSurroundingVal(qData, cluster, clusterVal=255):
    """ Summary:
            Check and see if the cluster is surrounded by a solid value, if it
            is, return the value which surrounds it. This one checks for border
            conditions
        Parameters:
            qData: The quantized satellite data.
            cluster: A 2D array of shape (2,n) where n is the number of points
                composing the cluster. The two rows are parallel arrays of x
                then y coordinates respectively.
        Returns:
            theValue: The surrounding value found.
    """

    theValue = None

    lk = set([(0, 1), (1, 1), (1, 0), (-1, 0), (-1, 1)])  # Left Kernel
    rk = set([(1, -1), (1, 0), (-1, 0), (-1, -1), (0, -1)])  # Right Kernel
    tk = set([(1, -1), (0, 1), (1, 1), (0, -1), (1, 0)])  # Top Kernel
    bk = set([(-1, -1), (0, 1), (-1, 1), (0, -1), (-1, 0)])  # Bottom Kernel
    mk = np.array(tuple(lk.union(rk)))  # Middle Kernel
    tlk = np.array(tuple(tk.intersection(lk)))  # Top Left Kernel
    trk = np.array(tuple(tk.intersection(rk)))  # Top Right Kernel
    blk = np.array(tuple(bk.intersection(lk)))  # Bottom Left Kernel
    brk = np.array(tuple(bk.intersection(rk)))  # Bottom Right Kernel
    lk = np.array(tuple(lk))
    rk = np.array(tuple(rk))
    tk = np.array(tuple(tk))
    bk = np.array(tuple(bk))

    rows, cols = qData.shape
    neighbors = {}

    for pt in cluster.T:
        row, col = pt
        k = None
        kused = None
        if row > 0:
            if col > 0:
                if row < rows - 1:
                    if col < cols - 1:
                        k = (mk + pt).T
                        kused = "middle"
                    else:
                        k = (rk + pt).T
                        kused = "right"
                elif col < cols - 1:
                    k = (bk + pt).T
                    kused = "bottom"
                else:
                    k = (brk + pt).T
                    kused = "bottom right"
            elif row < rows - 1:
                k = (lk + pt).T
                kused = "left"
            else:
                k = (blk + pt).T
                kused = "bottom left"
        elif col > 0:
            if col < cols - 1:
                k = (tk + pt).T
                kused = "top"
            else:
                k = (trk + pt).T
                kused = "top right"
        else:
            k = (tlk + pt).T
            kused = "top left"

        try:
            values = qData[tuple(k)]
        except:
            print "row %d out of %d" % (row, rows)
            print "col %d out of %d" % (col, cols)
            print "pt = ", pt
            print "k used = ", kused
            print "resulting k = ", k
            qData.fill(0)
            qData[tuple(pt)] = 255
            plt.imshow(qData)
            plt.axes().set_aspect('equal')
            plt.show()

        if clusterVal == 255:
            for value in values:
                if value == clusterVal or value == 0 or value == 255:
                    continue

                if theValue is None:
                    theValue = value
                elif value != theValue:
                    return None
        else:
            for value in values:
                if value == clusterVal or value == 0:
                    continue

                if value in neighbors.keys():
                    neighbors[value] += 1
                else:
                    neighbors[value] = 1

            maximum = 0
            for k, v in neighbors.iteritems():
                if v > maximum:
                    theValue = k

    if theValue is None:
        return clusterVal
    else:
        return theValue


def filter(img):
    """
    Summary: A filter to remove all points with less than 3 neighbors.
    Parameters:
        img - The image to be filtered
    Returns:
        output - The filtered dataset. Same dimensions as img
    """

    # Add the border mask
    cpy = np.ndarray((img.shape[0] + 2, img.shape[1] + 2))
    cpy[:, 0] = 0
    cpy[0, :] = 0
    cpy[:, -1] = 0
    cpy[-1, :] = 0
    cpy[1:-1, 1:-1] = img.copy()
    h, w = cpy.shape
    for y in xrange(1, h - 2):
        for x in xrange(1, w - 2):
            val = img[y, x]
            # Extract a 3x3 sub matrix
            submat = cpy[y - 1:y + 2, x - 1:x + 2]
            #print "submat \n", submat
            #print "kernel \n", kernel
            numOccurence = submat[submat == val].size

            if numOccurence < 3:
                cpy[y, x] = 0.0  # So remove it!

    output = np.where(cpy == 1)

    return output


def findMaximama(array):
    """ Find the maxima of a 1D array """

    binSizes = array[0]
    binVals = array[1][:-1]
    c = (np.diff(np.sign(np.diff(binSizes))) < 0).nonzero()[0] + 1  # local max

    # Set to True to plot the maxima
    if False:
        plt.plot(binVals, binSizes)
        plt.plot(binVals[c], binSizes[c], "o", label="max")
        plt.legend()
        plt.axes().set_aspect('equal')
        plt.show()

    return binVals[c]


def extractStepEdges(img):
    """
    Summary: A filter to find all step edges in the dataset:
        img - The image to be filtered
    Returns:
        output - The filtered dataset. Same dimensions as img
    """

    subImg = img[1:-1, 1:-1]
    pts = np.array(np.where(subImg == 1))

    yVals = []
    xVals = []
    for y, x in pts.T:
        y += 1  # Pts were found with borders extracted
        x += 1  # Pts were found with borders extracted
        val = img[y,x]
        # Extract a 3x3 sub matrix
        submat = img[y - 1:y + 2, x - 1:x + 2]
        targetPts = np.logical_and(submat != val, submat != 255)
        numOccurence = submat[targetPts].size

        if numOccurence > 0:
            yVals.append(y)
            xVals.append(x)

    edgePts = np.array((yVals, xVals))

    #img[(yVals, xVals)] = 8
    #plt.imshow(colorMap(img))
    #plt.show()

    return edgePts


def buildKDTree(data):
    global kdt, arr

    arr = np.asarray(zip(data[0], data[1]))
    kdt = kd.cKDTree(arr)


def getGrads(img, kernelSize):
    global gradx, grady, grad

    print "Computing xderiv"
    gradx = xDeriv(img, kernelSize)
    print "Computing yderiv"
    grady = yDeriv(img, kernelSize)


def getSobelKernel(size, direction):
    if size % 2 == 0:
        print "\nKernel size must be odd!\n"
        exit()

    magic = size / 2
    kernel = list()
    for i in xrange(-magic, magic + 1):
        row = list()
        if i < 0:
            pivot = i - magic
            for j in xrange(-magic, magic + 1):
                e = pivot + np.abs(j)
                row.append(e)
            kernel.append(row)
        elif i == 0:
            row = [0] * size
            kernel.append(row)
        else:
            pivot = i + magic
            for j in xrange(-magic, magic + 1):
                e = pivot - np.abs(j)
                row.append(e)
            kernel.append(row)

    if direction == 'y':
        return np.asarray(kernel)
    elif direction == 'x':
        return np.asarray(kernel).T


def xDeriv(img, kernelSize):

    s = kernelSize
    if s % 2 == 0:
        sys.exit("Kernel size must be odd!\n")
    kx = getSobelKernel(s, 'x')

    h, w = img.shape
    xgrad = np.zeros(img.shape)
    sigma = s / 2
    #blurImg = cv2.GaussianBlur(img, (s, s), sigma)
    for j in xrange(s / 2, h - s / 2):
        for i in xrange(s / 2, w - s / 2):
            submat = blurImg[j - s / 2: j + s / 2 + 1, i - s / 2: i + s / 2 + 1]
            val = (kx * submat).sum()
            xgrad[j, i] = val

    return xgrad


def yDeriv(img, kernelSize):

    s = kernelSize
    if s % 2 == 0:
        sys.exit("Kernel size must be odd!\n")
    ky = getSobelKernel(s, 'y')

    h, w = img.shape
    ygrad = np.zeros(img.shape)
    sigma = s / 2
    #blurImg = cv2.GaussianBlur(img, (s, s), sigma)
    for j in xrange(s / 2, h - s / 2):
        for i in xrange(s / 2, w - s / 2):
            submat = blurImg[j - s / 2: j + s / 2 + 1, i - s / 2: i + s / 2 + 1]
            val = (ky * submat).sum()
            ygrad[j, i] = val

    return ygrad


def preprocess(data, cutoff):

    global kdt

    numRegions = 8
    clustDist = 1

    print "Quantizing the Data...\n"
    qData = quantize(data, cutoff, numRegions)

    #plt.figure(0)
    #plt.title("After Quantization")
    #plt.imshow(colorMap(qData))

    # Cluster the Clouds using efficient KDT
    print "Clustering Clouds...\n"
    cloudsB4 = np.array(np.where(qData == 255))
    C = clusterLib.Cluster()
    cloudClusters = C.clusterGrid(cloudsB4, clustDist)

    # Remove the clouds once and prepare for filtering
    print "Removing clouds...\n"
    cloudsB4 = np.array(np.where(qData == 255))
    for cluster in cloudClusters:
        c = type(cluster)
        if c is int:
            c = [cluster]
        else:
            c = list(cluster)

        cloud = cloudsB4[:, c]
        val = findSurroundingVal(qData, cloud)
        if val is not None:
            qData[tuple(cloud)] = val

    #plt.figure(1)
    #plt.title("After Removing Clouds")
    #plt.imshow(colorMap(qData))

    # Filter the data
    print "Filtering out Small Clusters...\n"
    for region in xrange(numRegions):
        clusterVal = numRegions - region
        if not clusterVal in qData:
            continue
        dataPts = np.array(np.where(qData == clusterVal))
        print "\tFiltering %d Data Points in Region %d..." % \
            (dataPts[0].size, clusterVal)
        C = clusterLib.Cluster()
        regionClusters = C.clusterGrid(dataPts, clustDist)
        for cluster in regionClusters:
            if type(cluster) is int or len(cluster) > 250:
                continue
            r = dataPts[:, tuple(cluster)]
            val = findSurroundingVal(qData, r, clusterVal)
            qData[tuple(r)] = val

    #plt.figure(2)
    #plt.title("After Filtering Data")
    #plt.imshow(colorMap(qData))

    # Remove the clouds again
    print "Removing Clouds Again..."
    cloudsB4 = np.array(np.where(qData == 255))
    C = clusterLib.Cluster()
    cloudClusters = C.clusterGrid(cloudsB4, clustDist)

    for cluster in cloudClusters:
        c = type(cluster)
        if c is int:
            c = [cluster]
        else:
            c = list(cluster)
        try:
            cloud = cloudsB4[:, c]
        except:
            embed()
        val = findSurroundingVal(qData, cloud)
        if val is not None:
            qData[tuple(cloud)] = val

    #plt.figure(3)
    #plt.title("After Removing Clouds After Filtering Data")
    #plt.imshow(colorMap(qData))
    #plt.axes().set_aspect('equal')
    #plt.show()

    #print "Getting Gradient Values...\n"
    #if gradx is None or grady is None:
        #getGrads(qData, ptsBtwnKnts)

    print "Extracting Edges...\n"
    edgePts = extractStepEdges(qData)

    #colorData = colorMap(qData)
    #plt.imshow(colorData)
    #plt.axes().set_aspect('equal')
    #plt.show()

    return edgePts
