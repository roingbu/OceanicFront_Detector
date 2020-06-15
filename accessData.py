"""
netcdf sschlor.ncom.GOMh0.04.012.nest1.2009101000_00000000 {
dimensions:
    Time = UNLIMITED ; // (1 currently)
    Latitude = 349 ;
    Longitude = 517 ;
variables:
    float sschlor(Time, Latitude, Longitude) ;
        sschlor:_FillValue = -1.e+34f ;
        sschlor:long_name = "Sea Surface Chlorophyll" ;
        sschlor:units = "mg/m^3" ;
        sschlor:FORTRAN_format = "f7.3" ;
    int Time(Time) ;
        Time:long_name = "Time" ;
        Time:units = "seconds since 2009-10-10 00:00:00.00 0:00" ;
    float Latitude(Latitude) ;
        Latitude:long_name = "Latitude" ;
        Latitude:units = "degrees_north" ;
        Latitude:FORTRAN_format = "f8.2" ;
    float Longitude(Longitude) ;
        Longitude:long_name = "Longitude" ;
        Longitude:units = "degrees_east" ;
        Longitude:FORTRAN_format = "f8.2" ;

// global attributes:
        :Conventions = "COARDS" ;
        :title = "NCOM1 GOMh0.04 01.2" ;
        :institution = "Naval Research Laboratory, Stennis Space Center" ;
        :source = "NCOM native output" ;
        :experiment = "01.2" ;
        :history = "ncom2nc" ;
}
"""

#import netCDF4 as nc
import h5py as h5
import numpy as np

def getNCData(filename=None):

    data = nc.Dataset(filename, 'r')
    ncVars = data.variables
    sschlor = ncVars["sschlor"]
    maskNarray = sschlor[0]
    allVals = maskNarray.data
    mask = maskNarray.mask
    goodVals = allVals[~mask]

    return (allVals, goodVals)


def getHDFData(filename=None):
    """ The HDF data has all the nans and land values set to 0. This
    is changed below for improved visualization.
    land = 0, clouds = -1"""

    dataSet = h5.File(filename, mode="r")
    data = dataSet.get("data").value
    landMask = np.load("/home/nathan/code/theses/masters_thesis/data/landMask.npy")
    # A precalculated land mask.
    data[np.where(data == 0)] = -1  # Make all fills -1
    data[np.where(landMask == 0)] = 0  # Then set land back to 0!

    return data
    # The below numbers are the range to isolate the gulf of mexico
    #return vals[1350:1750, 1950:2500]

