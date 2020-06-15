import os
import preprocessing as pp  # I've rather ungracefully bundled all the
                            # preprocessing routines into one file
from accessData import getHDFData  # Routine to get the hdf data

# A handy library to embed an interpreter in the middle of a script for debugging
# simply call embed() anywhere in the program and it will pause there and launch
# an interactive interpreter
from IPython import embed

# The below packages are all for plotting and color mapping
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as m

# Will operate on every file in this directory
rootDir = "/home/nathan/code/theses/masters_thesis/data/chuanminsData/"
# The resulting images will be saved here
rsltsDir = "/home/nathan/code/theses/masters_thesis/results/testing_cutoffs_gray/"
# All land values for SECOORA region
landMask = pp.np.load("/home/nathan/code/theses/masters_thesis/data/landMask.npy")

landMask = landMask[750:, :1350]  # cropped the SECOORA to be mostly GOM
landMask = landMask.astype(pp.np.bool)  # Convert ints to bools for masked array

colors = ['blueviolet', 'blue', 'aqua', 'limegreen', 'greenyellow',
          'yellow', 'orange', 'red']

#for f in os.listdir(rootDir):
#files = ["A20102032010209.1KM.SECOORA.7DAY.L3D.OCI.h5",
#         "A20102912010297.1KM.SECOORA.7DAY.L3D.OCI.h5",
#         "A20101442010150.1KM.SECOORA.7DAY.L3D.OCI.h5"]  # Three data sets from
                                                         # summer, fall, and spring

files = ["A20102032010209.1KM.SECOORA.7DAY.L3D.OCI.h5"]

cutoffs = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
           0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]  # List of all cutoff points to test

c = 0
for f in files:
    for cutoff in cutoffs:

        print "Loading Datafile %s..." % (f)
        data = getHDFData(rootDir + f)
        data = data[750:, :1350]  # Crop the SECOORA to the GOM
        cols, rows = data.shape
        masked = pp.np.ma.array(data, mask=~landMask)  # For easy visualization

        # Create custom colormap...
        minv = data[data > 0].min()
        maxv = data.max()
        rng = maxv - minv
        thresh = [0.0, 0.05 / rng, 0.1 / rng, 0.25 / rng, 0.5 / rng, 1.3 / rng,
                    5.0 / rng, 1.0]  # All colors are scaled between 0 and 1
        clr_indx = zip(thresh, colors)
        cmap = m.colors.LinearSegmentedColormap.from_list('custom', clr_indx, 256)
        if True:  # Gordon requested we examine the fronts in gray scale mapping
            cmap = cm.get_cmap('gray')  # overwrite the previos color mapping
            minv = 0.01
            maxv = 0.2
        cmap.set_bad('black')  # Set the land to black (from masked)
        cmap.set_under('white')  # Set the clouds to white. set_under
                                 # means set all the values under minval to
                                 # white. (clouds are -1 and the minval is > 0.
                                 # Roughly 0.011)

        cmap.set_over('white')  # set_over does the opposite of set_under

        edgePts = pp.preprocess(data, cutoff)  # This is the actual routine which
                                          # finds the 'front'. The edge pixels
                                          # around each voronoi region are the
                                          # approximated front

        masked.data[edgePts.tolist()] = 255  # Here we set all the edge pixels
                                             # to 255, or white thanks to the
                                             # cmap.set_over routine.

        plt.imshow(masked, cmap=cmap, vmin=minv, vmax=maxv)  # create plot
        fig = plt.gcf()  # get current figure

        fig.set_size_inches((2 * cols) / 100.0, (2 * rows) / 100.0)  # set save
                                                                     # parameters
        plt.savefig(rsltsDir+f[:-3]+"_%f.png" % (c), dpi=100)
        # By combining set_size and the dpi in savefig, you have complete
        # control over the size and resolution of the saved figure. If the
        # denominator of the fraction in set_size is the same as dpi in savefig
        # then the numerator of the set_size function denotes the resolution of
        # the saved image.

        plt.clf()  # clear current figure
        c += 1
