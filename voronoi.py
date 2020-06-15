import numpy as np
import sys
from matplotlib import colors as mplColors
import matplotlib.pyplot as plt
from IPython import embed
from accessData import getHDFData


def findGens(data, numRegions):

    if len(data) == 0:
        print "data is empty"
        exit()
    partitions = numRegions + 2

    # Scott's (1979) bin width. Based on the standard deviation
    # and the data size
    data_std = np.std(data)
    n = np.power(len(data), 1.0 / 3)
    bin_width = (3.49 * data_std) / n
    histMin = np.min(data)
    histMax = np.max(data)
    bins = np.arange(histMin, histMax + bin_width, bin_width)

    cdf, edges, patches = plt.hist(data, bins=bins,
                                   cumulative=True, normed=True)

    minVal = np.min(cdf)
    maxVal = np.max(cdf)
    yVals = np.linspace(minVal, maxVal, partitions)

    gens = []
    for val in yVals[1:-1]:
        loc = np.argmin(np.absolute(cdf - val))
        if cdf[loc] > val:
            pt1 = cdf[loc - 1]
            pt2 = cdf[loc]
            if pt1 == pt2:
                gens.append(edges[loc])
                continue
            ratio = (val - pt1) / (pt2 - pt1)
            x = edges[loc - 1] + ratio * (edges[loc] - edges[loc - 1])
            gens.append(x)
        elif cdf[loc] < val:
            pt1 = cdf[loc]
            pt2 = cdf[loc + 1]
            if pt1 == pt2:
                gens.append(edges[loc + 1])
                continue
            ratio = (val - pt1) / (pt2 - pt1)
            x = edges[loc] + ratio * (edges[loc + 1] - edges[loc])
            gens.append(x)
        else:
            gens.append(edges[loc])

    if False:
        retVal = np.array((yVals[1:-1], gens))
        for b in retVal.T:
            print "b = ", b
            plt.plot((edges[0], b[1]), (b[0], b[0]), 'g')  # horizontal
            plt.plot((b[1], b[1]), (minVal, b[0]), 'r')  # vertical
        plt.title("Optimal Starting Generators", size=18)
        plt.xlabel("Chlorophyll Intensity (mg/m^3)", size=12)
        plt.ylabel("Cumulative Probability", size=12)
        plt.figure(1)
        d2 = data[data<1]
        data_std = np.std(d2)
        n = np.power(len(d2), 1.0 / 3)
        bin_width = (3.49 * data_std) / n
        histMin = np.min(d2)
        histMax = np.max(d2)
        bins = np.arange(histMin, histMax + bin_width, bin_width)
        plt.title("Chlorophyll Data Histogram", size=18)
        plt.xlabel("Chlorophyll Intensity (mg/m^3)", size=12)
        plt.ylabel("Relative Frequency", size=12)
        print "bins = ", bins
        hist, edges, patches = plt.hist(d2, bins=bins)
        plt.show()
        exit()
    else:
        plt.clf()

    return np.asarray(gens)


def voronoi(inputSet, numRegions=8, tolerance=0.0001, data=None):

# ----------------------------------------------------------------
    frame = 0
    #filename = "figures/cvts/%dreg=%d_tol=%f.png" % (frame, numRegions, tolerance)
    filename = "/home/nathan/code/theses/masters_thesis/figures/cvts/jpeg/%04d.jpeg" % (frame)
# ----------------------------------------------------------------

    gens = findGens(inputSet, numRegions)
    #print "\tInitial %d Generators\n\t" % (numRegions)
    #print gens
    while len(set(gens)) != numRegions:  # Ensure generator uniqueness
        print "Duplicate Generators Found. Finding new Generators"
        gens = findGens(inputSet, numRegions)
        print "\tNew Voronoi Region Generators\n\t", gens

    regions = {gen: [] for gen in gens}

    for s in inputSet:
        distance = np.absolute(gens - s)  # 1D Euclidean distance from point s
        # to all the generators
        closestGen = gens[np.argmin(distance)]
        regions[closestGen].append(s)  # Add the point s to the region around
        # the generator, closestGen

    oldGens = gens.copy()
    for i in xrange(gens.size):
        #gens[i] = massCentroid(regions[gens[i]])  # Update the
        # generators
        gens[i] = np.average(regions[gens[i]])

    #norm = np.linalg.norm(gens - oldGens)
    norm = np.max(np.absolute(gens - oldGens))
    print "Tolerance = ", tolerance
    print "norm = ", norm
    while norm > tolerance:
    #for i in xrange(22):
# ----------------------------------------------------------------
        thresh = []
        colors = ['black', 'blueviolet', 'blue', 'aqua', 'limegreen', 'greenyellow',
                  'yellow', 'orange', 'red', 'white']
        data[data == -1] = 254
        cmap = mplColors.ListedColormap(colors)
        for val in regions.values():
            thresh.append(np.min(val))
        thresh = np.sort(thresh)
        bounds = [0]
        bounds.extend(thresh)
        bounds.extend([50,255])
        bnorm = mplColors.BoundaryNorm(bounds, cmap.N)

        #(rows, cols) = data.shape
        #qData = np.ndarray((rows, cols), dtype='uint8')
        #qData.fill(255)

        #for val, t in enumerate(thresh):
            #qData[data > t] = val + 1

        #qData[data == 0] = (0, 0, 0)  # Put back the land values
        #qData[data == -1] = (255, 255, 255)  # Put back the cloud values

        plt.clf()
        plt.title("Rho - %d Gens, L2 Norm = %1.04f" % (numRegions, norm))
        img = plt.imshow(data, cmap=cmap, norm=bnorm)
        plt.colorbar(img, cmap=cmap, norm=bnorm, boundaries=bounds, ticks=bounds)
        #plt.show()
        #exit()
        plt.savefig(filename)
        frame += 1
        #filename = "figures/cvts/%03dreg=%d_tol=%1.02f.png" % \
                   #(frame + 1, numRegions, tolerance)
        filename = "/home/nathan/code/theses/masters_thesis/figures/cvts/jpeg/%04d.jpeg" % (frame)

        data[data == 254] = -1
# ----------------------------------------------------------------

        regions = {gen: [] for gen in gens}
        for s in inputSet:
            distance = np.absolute(gens - s)  # 1D Euclidean distance from
            # point s to all the generators
            closestGen = gens[np.argmin(distance)]
            regions[closestGen].append(s)  # Add the point s to the region
            # around the generator, closestGen

        oldGens = gens.copy()
        for i in xrange(gens.size):
            #gens[i] = massCentroid(regions[gens[i]])  # Update the
            gens[i] = np.average(regions[gens[i]])
            # generators

        #norm = np.linalg.norm(gens - oldGens)
        norm = np.max(np.absolute(gens - oldGens))
        print "norm = ", norm

    return regions


def massCentroid(data):

    if len(data) == 0:
        print "data is empty"
        exit()

    partitions = 75

    # Scott's (1979) bin width. Based on the standard deviation
    # and the data size
    data_std = np.std(data)
    n = np.power(len(data), 1.0 / 3)
    bin_width = (3.49 * data_std) / n
    histMin = np.min(data)
    histMax = np.max(data)
    bins = np.arange(histMin, histMax + bin_width, bin_width)
    if histMax - histMin == 0:
        print "len(data) = ", len(data)
        print "histMin = ", histMin
        print "histMax = ", histMax
        print "bins = ", bins
        print "bin_width = ", bin_width
    cdf, edges, patches = plt.hist(data, bins=bins, color="blue",
                                   cumulative=True, normed=True)

    # Extrapolate to find the next value in the CDF
    x1 = edges[-2] - edges[-3]
    x2 = edges[-1] - edges[-3]
    y1 = cdf[-1] - cdf[-2]
    y2 = (y1 * x2) / x1
    lst = cdf.tolist()
    lst.append(cdf[-2] + y2)
    cdf = np.asarray(lst)

    minVal = np.min(cdf)
    maxVal = np.max(cdf)
    yVals = np.linspace(minVal, maxVal, partitions)

    optimalBins = []
    for val in yVals:
        loc = np.argmin(np.absolute(cdf - val))
        if cdf[loc] > val:
            pt1 = cdf[loc - 1]
            pt2 = cdf[loc]
            if pt1 == pt2:
                optimalBins.append(edges[loc])
                continue
            ratio = (val - pt1) / (pt2 - pt1)
            x = edges[loc - 1] + ratio * (edges[loc] - edges[loc - 1])
            optimalBins.append(x)
        elif cdf[loc] < val:
            pt1 = cdf[loc]
            pt2 = cdf[loc + 1]
            if pt1 == pt2:
                optimalBins.append(edges[loc + 1])
                continue
            ratio = (val - pt1) / (pt2 - pt1)
            x = edges[loc] + ratio * (edges[loc + 1] - edges[loc])
            optimalBins.append(x)
        else:
            optimalBins.append(edges[loc])

    hist, edges, patches = plt.hist(data, optimalBins)
    plt.figure(0)
    if False:
        cdf, edges, patches = plt.hist(data, bins=bins,
                                    cumulative=True, normed=True)
        retVal = np.array((yVals, optimalBins))
        plt.figure(0)
        for b in retVal.T:
            plt.plot((edges[0], b[1]), (b[0], b[0]), 'g')  # horizontal
            plt.plot((b[1], b[1]), (minVal, b[0]), 'r')  # vertical
        plt.title("Optimal Bins for Density Histogram", size=18)
        plt.xlabel("Chlorophyll Intensity (mg/m^3)", size=12)
        plt.ylabel("Cumulative Probability", size=12)

        plt.figure(1)
        plt.title("Histogram Approximating Density", size=18)
        plt.xlabel("Chlorophyll Intensity (mg/m^3)", size=12)
        plt.ylabel("Relative Frequency", size=12)
        plt.show()

    num = 0.0  # Numerator
    den = 0.0  # Denominator
    for val in data:
        distance = np.absolute(val - edges[:-1])  # Find the distance of val
        # from every bin
        loc = np.argmin(distance)
        edge = edges[loc]
        if val < edge:
            loc -= 1
        density = hist[loc]
        #num += val * density * density * density
        #den += density * density * density
        #num += val * density * density
        #den += density * density
        num += val * density
        den += density

    centroid = num / den

    return centroid


if __name__ == "__main__":

    #FILENAME = "../extract_splines/data/crop_square.png"
    #print "Opening ", FILENAME
    #img = np.array(Image.open( FILENAME )
    print "Getting HDF Data..."
    img = getHDFData()

    data = img.flatten()
    print "Finding the Voronoi Regions of the flattened HDF data..."
    r = voronoi(data[np.logical_and(data > 0, data < 10)])
    print "Final Thresholds = ", r.keys()
    #exit()
    downSampled = img.copy()
    for key, vals in r.iteritems():
        minVal = np.min(vals)
        maxVal = np.max(vals)
        downSampled[np.logical_and(img > minVal, img < maxVal)] = key

    plt.clf()
    plt.figure(1)
    plt.imshow(img, cmap="gray")
    plt.figure(2)
    plt.imshow(downSampled, cmap="gray")
    plt.show()
