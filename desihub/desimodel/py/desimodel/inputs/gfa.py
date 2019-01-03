# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

#import desimodel.io
#import scipy.interpolate
import numpy as np
from ..focalplane import get_radius_deg, get_radius_mm

def build_gfa_table(outfile = 'gfa.ecsv'):
    '''
    Builds the GFA table given the data from DESI-0530-v13 Excel spreadsheet
    and writes a .ecsv file using the astropy table library. The data is pulled from
    the "GFALocation" tab on the spreadsheet and from rows 16-23 and columns A-I.
    Parameters
    ----------
    outfile: a default parameter that represents the desired filename which is returned by 
    this function. The filename defaults to "gfa.ecsv" if no parameters are given. 
    '''
    # Uses the reference projection of active area to create data table of GFAs
    from astropy.table import Table
    # Initial x and y coordinates for the GFAs
    """ Uses the x and y from the petal indexed at 9 so the first petal 
        added to the table is indexed at 0
    [[-125.10482863 -370.01790486]
     [-129.83038525 -384.56223777]
     [-159.04283509 -375.05643893]
     [-154.31646944 -360.51151824]]   
    """
    # Data obtained from DESI-0530-v13 Excel spreadsheet
    x = [-125.10482863, -129.83038525, -159.04283509, -154.31646944]
    y = [-370.01790486, -384.56223777, -375.05643893, -360.51151824]
    z = [-17.053, -18.487, -18.631, -17.198]
    rotatemat = np.zeros(shape=(2,2))
    rotatemat[0] = [np.cos(36*np.pi/180), -np.sin(36*np.pi/180)]
    rotatemat[1] = [np.sin(36*np.pi/180), np.cos(36*np.pi/180)]
    
    # Note: the corners are 0 indexed 
    gfatable = Table(names = ('PETAL', 'CORNER', 'X', 'Y', 'Z', 'Q', 'RADIUS_DEG', 'RADIUS_MM'), 
                 dtype = ('int', 'int', 'float', 'float', 'float', 'float', 'float', 'float'))
    # Sets the units for the GFA table
    gfatable['X'].unit = 'mm'
    gfatable['Y'].unit = 'mm'
    gfatable['Z'].unit = 'mm'
    gfatable['Q'].unit = 'degrees'
    gfatable['RADIUS_DEG'] = 'degrees'
    gfatable['RADIUS_MM'] = 'mm'
    find_gfa_coordinates(x, y, z, gfatable, rotatemat)
    # Saves the table of data as an ecsv file
    gfatable.write(outfile, format='ascii.ecsv')
    
# Function that obtains the x and y coordinates for each corner of the GFAs
def find_gfa_coordinates(x, y, z, gfatable, rotatemat):
    '''
    Finds all the GFA coordinates by rotating the initial coordinates and adding
    the respective coordinates to the gfaTable
    Parameters
    ----------
    x : Array of four x initial coordinates of the GFA
    y : Array of four y initial coordinates of the GFA
    z: Array of four z initial coordinates of the GFA
    gfaTable: Astropy Table object which stores the petal number, corner number, 
    and x, y, and z coordinates in mm within each row
    rotateMat: Rotation matrix to rotate the coordinates 36 degrees counterclockwise
    '''
    coord = np.zeros(shape=(2,1))
    gfacoord = np.zeros(shape=(4, 2))
    oldxcoord = x
    oldycoord = y
    for j in range(10):
        for i in range(4):
            coord[0] = oldxcoord[i]
            coord[1] = oldycoord[i]
            newcoord = np.matmul(rotatemat, coord)
            oldxcoord[i] = newcoord[0]
            oldycoord[i] = newcoord[1]
            gfacoord[i] = [newcoord[0], newcoord[1]]
            
            theta = np.degrees(np.arctan2(newcoord[0], newcoord[1]))
            # radius is the radius in mm
            radius = np.sqrt(newcoord[0]**2 + newcoord[1]**2)
            # degree is the radius in degrees
            degree = get_radius_deg(newcoord[0], newcoord[1])
            # Could be building the table in O(N^2), which is notably inefficient
            gfatable.add_row([j, i, newcoord[0], newcoord[1], z[i], theta, degree, radius])


            
