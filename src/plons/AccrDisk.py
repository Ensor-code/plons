import os
import math                     as math
import numpy                    as np

# import plons scripts
import plons.ConversionFactors_cgs    as cgs
import plons.GeometricalFunctions     as gf
import plons.Tools                    as tl
import plons.PhysicalQuantities       as pq

# from plons.SmoothingKernelScript import rotatePixCoordAroundZ

# ignore warnings
import warnings
warnings.filterwarnings("ignore")




'''
To Do
- change coordinate system to (0,0,0) at location companion sink 
- only select data with r <... ? probably just do calculations up to certain r? (rotate both sinks on x-axis? not nescessary?)

- estimate radius of accretion disk, by cut-off on radial density profile (with r=0 at companion sink):
    use smoothingKernelScript !things that have to be updated: position, AGBcoord, compcoord, posAGB!
    !In smoothinkernelscript centred around AGB star?
    which direction? multiple directions? mean in all directions?

- at that radius, use similar method to estimate scale height in z-direction 
- estimate mass, with a fit for the disk geometry 
    For a vertically isothermal disc the self-consistent solution is T ∝ R^−1/2 , which gives H/R ∝ R^(5/4)
     https://rdalexander.github.io/planets_2022/lecture2_notes.pdf 

- Evolution in time of disk radius, height and mass? But you need dumpdata to calculate, and we only read in one dump

- Inflowing mass in disk (to compare to mass accretion onto sink)?

Lee et al:

- estimate disk height at several r: for each r 10 vertical lines and find on each line drop in density 
    (no longer exponential function) and take mean of the 10
- estimate disk mass at several radii (with their corresponding scale heights), until mass increment is negligible 
--> to find M and radius
- do this for each timestep :p

'''

def change_origin_coordinates(x_coords,y_coords,z_coords, new_origin):
    # Get the current origin
    current_origin = (0, 0, 0)

    # Calculate the translation vector
    translation_vector = (new_origin[0] - current_origin[0],
                          new_origin[1] - current_origin[1],
                          new_origin[2] - current_origin[2])

    # Apply the translation to each coordinate array
    new_x_coords = x_coords + translation_vector[0]
    new_y_coords = y_coords + translation_vector[1]
    new_z_coords = z_coords + translation_vector[2]

    return new_x_coords, new_y_coords, new_z_coords


'''
For rotation:
    theta = pq.getPolarAngleCompanion(dumpData['posComp'][0], dumpData['posComp'][1])

def getPixels(shape, n, r, data, bound):
Get a 2D array containing the [x,y,z]-components of all pixels 
          shape == 'r' --> spherical splice
          shape == 'z' --> planar slice: orbital plane (xy-plane)
          shape == 'x' --> yz-plane
          shape == 'y' --> xz-plane

pixCoord    = getPixels(shape, n, r, data, bound)

# Rotate the plane/line to align with companion
if shape != 'line_z':
    pixCoord = rotatePixCoordAroundZ(theta, pixCoord)

'''

'''
main definition

INPUT:
    - dumpData      data of dump file   [dictionary]
    - setup         setup data          [dictionary]
    - run           run number          [string]    
    - loc           output directory    [string]
'''

def accrDiskAnalysis(run, loc, dumpData, setup):

    dumpData['new_position']    = np.array(change_origin_coordinates(dumpData['position'].transpose()[0],dumpData['position'].transpose()[1],dumpData['position'].transpose()[2], dumpData['posComp'])).transpose()
    new_r, new_phi, new_theta   = gf.TransformToSpherical(dumpData['new_position'].transpose()[0],dumpData['new_position'].transpose()[1],dumpData['new_position'].transpose()[2])              # sperical coordinates





