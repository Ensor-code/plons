import numpy    as np
import math     as math


'''
README:

This script contains geometrical function for defining planes, shells,
distances to planes, normals of planes,... needed for the first script
ever to plot PHANTOM data. Many of these functions might not be
relevant anymore, but might ever come in handy.

If you ever need to define planes analytically with python, this is where to look :)

Useful functions:
    - coordinate transformation from Cartesian to spherical
    - get radius from a certain point, not the origin
    - get radial and tangential velocity
'''


# ------------------------------------------------------------



'''
Perform a coordinate transformation (just a rotation) of two dimensions (x,y) with a given angle alpha. ---> input is array OR floats
  Can also input (x,z) or (y,z) if this is needed.
'''
def coordTransf(x,y,alpha):
    if x is list:
        x_transf = []
        y_transf = []
        for i in range(len(x)):
            x_transf.append(x[i]*np.cos(alpha)+y[i]*np.sin(alpha))
            y_transf.append(-x[i]*np.sin(alpha) +y[i]*np.cos(alpha))
    else:
        x_transf = x*math.cos(alpha)+y*math.sin(alpha)
        y_transf = -x*math.sin(alpha) +y*math.cos(alpha)
    return x_transf, y_transf


'''
Get the radial and tangential velocity
  x, y, and phi are lists
'''
def getRadTanVelocity(x,y,v_x,v_y):
    phi = np.arctan2(y,x)
    v_rad = (coordTransf(v_x,v_y,phi)[0])
    v_tan = (coordTransf(v_x,v_y,phi)[1])
    return v_rad,v_tan

'''
Get ratio between the radial and tangential velocity
'''
def getRatioRadTan(v_r,v_t):
    ratio = []
    for i in range(len(v_r)):
        ratio.append(v_r[i]/v_t[i])
    return np.abs(ratio)

'''
Return the position of all particles as a function of radius with the AGB star as centre
'''
def getRadiusCoordinate(position,AGBcoord):
    position = position.transpose()
    x = position[0]
    y = position[1]
    z = position[2]
    r = np.sqrt((x-AGBcoord[0])**2+(y-AGBcoord[1])**2+(z-AGBcoord[2])**2)
    return r

'''
Calculate radius from the centre of mass
'''
def calc_r(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    return r

'''
Calculate radius from the centre of mass
'''
def calc_r_2D(x,y):
    r = np.sqrt(x**2+y**2)
    return r

'''
Calculate azimuthal angle theta
      Returns np.array
'''
def calcTheta(x,y,z):
    return np.arctan2(np.sqrt(x * x + y * y), z)


'''
Calculate azimuthal angle theta
      Returns np.array
'''
def calcPhi(x,y):
    return np.arctan2(y, x)


######-----    SHELL   -----######

'''
Get a the shell you want to model:
  Returns a list of indices which belong to the chosen shell
      r = [int] a certain radius where you want to take the shell
      dr = [int] thickness of the shell
      o = centre of your shell
      x,y,z [arrays] = cartesian coordinates
'''
def TakeShell(r, dr, o, x, y, z):
    index = []
    for i in range(len(x)):
        R = np.sqrt((x[i]-o[0])**2+(y[i]-o[1])**2+(z[i]-o[2])**2)
        if R <= r+dr and R >= r-dr:
            index.append(i)
    return index

'''
Select the data which is positioned within a given region ('slice'):
  Given an index list, this function returns the data which corresponds to the indices,
  and thus belongs to the slice.
      boolean == True for sperical selections
              == False for cartesian slices
'''
def getDataInSlice(index, data, boolean, sphericalCoords):
    if len(index) != 0:
        dataInSlice = []
        for i in range(len(index)-2):
            dataInSlice.append(data[:,index[i]])
        if boolean == True:
            rInSlice = []
            for i in range(len(index)-2):
                rInSlice.append(sphericalCoords[:,index[i]])
            return np.matrix.transpose(np.array(dataInSlice)), np.matrix.transpose(np.array(rInSlice))
        else:
            return np.matrix.transpose(np.array(dataInSlice))
    else:
        print('No data in this shell!')

'''
Make X and Y grids in order to plot the plane
'''
def makeGrids(linspace_x,linspace_y):
    X, Y = np.meshgrid(linspace_x, linspace_y)
    return X,Y

'''
Get z values to that the chosen shells can be plotted
      x,y are cartesian coordinates
      r, dr chosen specifics for the shell
      o is centre of the shells
'''
def makeShell(x,y,r,dr,o):
    z_innerShell = np.sqrt((r-dr)**2-(x-o[0])**2-(y-o[1])**2)+o[2]
    z_outerShell = np.sqrt((r+dr)**2-(x-o[0])**2-(y-o[1])**2)+o[2]
    return z_innerShell, z_outerShell

'''
Get the y-values for a 2D circle with given radius (r), x-coordinates (x, array) and centre (c = (x0,y0)).
'''
def getCircle(x,centre,r):
    cmin = []       # positive root of circle
    cplus = []      # negative root of circle
    for i in range(len(x)):
        cmin.append(-np.sqrt(r**2-(x[i]-centre[0])**2)+centre[1])
        cplus.append(np.sqrt(r**2-(x[i]-centre[0])**2)+centre[1])
    return cmin,cplus


######-----    PLANE    -----#######


### GENERAL PLANE STUFF

'''
Function to define a plane (orbital plane: xy-plane).
      p1, p2, p3 are arrays with three elements, which are the coordinates in 3D space
      general equation of a plane:  ax + by + cz + d = 0
      returns [a, b, c, d]
'''
def makePlane(p1,p2,p3):
    vect1 = []
    vect2 = []
    for i in range(3):
        vect1.append(p1[i]-p2[i])
        vect2.append(p1[i]-p3[i])
    # the elements of this 'norm'vector are a, b, c
    norm = [vect1[1]*vect2[2]-vect1[2]*vect2[1], -vect1[0]*vect2[2]+vect1[2]*vect2[0], vect1[0]*vect2[1]-p1[1]*vect2[0]]
    d = -norm[0]*p1[0]-norm[1]*p1[1]-norm[2]*p1[2]
    return [norm[0],norm[1],norm[2],d]

'''
Use the formula E = ax+by+cz+d and see what the result is. Fill in a point,
      if E = 0, then this point is located on the plane.
'''
def pointOnPlane(p,plane):
    result = plane[0]*p[0]+plane[1]*p[1]+plane[2]*p[2]+plane[3]
    return result

'''
Get the new value for d (when a plane is described by ax+by+cz+d=0) by using the formula
      D = |ax_1+by_1+cz_1+d|/sqrt(a^2+b^2+c^3), when D is the input [AU].
      a, b and c don't change when you consider perpendicular planes.
      Use for 'point' a point in the perpendicular plane (like eg. the AGB star)
  NOTE: there is also a second solution to this equation: -d
'''
def getPlaneAtDistance(D, plane, point):
    d = np.sqrt(plane[0]**2+plane[1]**2+plane[2]**2)*D-plane[0]*point[0]-plane[1]*point[1]-plane[2]*point[2]
    return d


def getDistancePlane(plane, point):
    D = (plane[0]*point[0]+plane[1]*point[1]+plane[2]*point[2]+plane[3])/(np.sqrt(plane[0]**2+plane[1]**2+plane[2]**2))
    return D


### ORBITAL PLANE

'''
Make two different planes, which define the range of the slice: z +/- dz for the orbital plane
'''
def makeOrbitalPlaneRange(p1,p2,p3,dz):
    # upper range: z_new = z+dz
    p1_plus = [p1[0], p1[1],p1[2]+dz]
    p2_plus = [p2[0], p2[1],p2[2]+dz]
    p3_plus = [p3[0], p3[1],p3[2]+dz]
    upperPlane = makePlane(p1_plus,p2_plus,p3_plus)
    # lower range: z_new = z-dz
    p1_min = [p1[0], p1[1],p1[2]-dz]
    p2_min = [p2[0], p2[1],p2[2]-dz]
    p3_min = [p3[0], p3[1],p3[2]-dz]
    lowerPlane = makePlane(p1_min,p2_min,p3_min)
    return upperPlane, lowerPlane

'''
The orbital plane is in the xy plane, so this means that when we take a slice in the orbital plane,
      we only need to have a look at the z-values of all coordinates.
  --> this function returns all rows in the data set that are in between a slice
      {only applicable on the orbital plane slice}
'''
def getRowsInSlice(z, dz):
    inSlice = []
    for i in range(len(z)):
        if z[i] <= dz and z[i] >= -dz:
            inSlice.append(i)
    return inSlice


#### PLANE PERPENDICULAR TO ORBITAL PLANE

'''
Get normal of the plane perpendicular to the orbital plane
      vector l is perpendicular with the normal of the orbital
      plane (n) and with the vector from the positions of the two stars (a)
      --> vect(n)*vect(l) = 0
          vect(a)*vect(l) = 0
'''
def makePerpPlane(p1,p2,p3):
    orbitalPlane = makePlane(p1,p2,p3)  # returns [normal, d]
    n = [orbitalPlane[0],orbitalPlane[1],orbitalPlane[2]]
    a = []
    for i in range(3):
        a.append(p1[i]-p2[i])
    l3 = 1
    l1 = ((n[2]*a[2]*l3)/(a[1]*n[0])-(n[2]*l3)/(n[0]))/(1-(n[1]*a[0])/(a[1]*n[0]))
    l2 = (-a[0]*l1-a[2]*l3)/(a[1])
    d = -l1*p1[0]-l2*p1[1]-l3*p1[2]
    return [l1,l2,l3,d]

'''
Get all rows in the data set of the points that lie between two planes perpendicular to the orbital plane.
      plane1 = orbitalPerpPlane with d-K
      plane2 = orbitalPerpPlane with d+K
'''
def getRowsInSlicePerpPlane(x,y,z,plane1,plane2):
    inSlice = []
    for i in range(len(x)):
        if pointOnPlane([x[i],y[i],z[i]],plane1) <= 0 and pointOnPlane([x[i],y[i],z[i]],plane2) >= 0:
            inSlice.append(i)
    return inSlice

'''
Function to define a plane (orbital plane: xy-plane).
      p1, p2, p3 are arrays with three elements, which are the coordinates in 3D space
      general equation of a plane:  ax + by + cz + d = 0
      returns [a, b, c], which is the normal of the plane given by the above formula
'''
def getNormalPlaneThrough(p1,p2,p3):
    vect1 = []
    vect2 = []
    for i in range(3):
        vect1.append(p1[i]-p2[i])
        vect2.append(p1[i]-p3[i])
    # the elements of this 'norm'vector are a, b, c
    norm = [vect1[1]*vect2[2]-vect1[2]*vect2[1], -vect1[0]*vect2[2]+vect1[2]*vect2[0], vect1[0]*vect2[1]-p1[1]*vect2[0]]
    return [norm[0],norm[1],norm[2]]

'''
Get normal of the plane perpendicular to the orbital plane
      vector l is perpendicular with the normal of the orbital
      plane (n) and with the vector from the positions of the two stars (a)
      --> vect(n)*vect(l) = 0
          vect(a)*vect(l) = 0
'''
def getNormalPerpTo(p1,p2,p3):
    n = getNormalPlaneThrough(p1,p2,p3)
    a = []
    for i in range(3):
        a.append(p1[i]-p2[i])
    l3 = 1
    l1 = ((n[2]*a[2]*l3)/(a[1]*n[0])-(n[2]*l3)/(n[0]))/(1-(n[1]*a[0])/(a[1]*n[0]))
    l2 = (-a[0]*l1-a[2]*l3)/(a[1])
    return [l1,l2,l3]


### FUNCTIONS TO PLOT THE PLANES

'''
Get the z component of the plane
'''
def makeZ(x,y,plane):
    z = (-plane[3]-plane[0]*x-plane[1]*y)/plane[2]
    return z

'''
Make a Z grid with the X,Y grid and the plane formula:    z = (-d-ax-by)/c
      X and Y are grids, thus Z too
'''
def makeZgrid(X,Y,plane):
    Z = (-plane[3]-plane[0]*X-plane[1]*Y)/plane[2]
    return Z

'''
Make a Y grid with the X,Z grid and the plane formula:    y = (-d-ax-cz)/b
      X and Y are grids, thus Z too
'''
def makeYgrid(X,Z,plane):
    Y = (-plane[3]-plane[0]*X-plane[2]*Z)/plane[1]
    return Y

'''
Get the distance between two points in 3D
'''
def getDistanceTwoPoints(p1,p2):
    x1=p1[0]
    y1=p1[1]
    z1=p1[2]
    x2=p2[0]
    y2=p2[1]
    z2=p2[2]
    return np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)


'''
Get the angle alpha, which is the angle between the perpendicular orbital plane and
the xz-plane (y=0) in order to do a coordinate transformation in order to plot this plane.
      - choose 3 points (p1, p2 and p3):
          p1 = AGB star
          p2 = (0,0,0)
          p3 = point in the plane: y3=0, z3=0, so choose a x3
'''
def getAlpha(p1,p2,p3):
    a = getDistanceTwoPoints(p1,p2)
    b = getDistanceTwoPoints(p1,p3)
    c = getDistanceTwoPoints(p2,p3)
    alpha = np.arccos((b**2+c**2-a**2)/(2*b*c))
    return alpha



'''
Get the projection for the perpendicular orbital plane:
  project the perp orb plane to the (y',z)-plane.
  The z value will stay the same, the y'-value will be equal to the "radius": y'=sqrt(x^2+y^2)
'''
def getProjection(x,y,z,plane):
    r = []
    for i in range(len(x)):
        r.append(getProjectionElement(x[i],y[i],z[i],plane))
    return r


def getProjectionElement(x,y,z,plane):
    # components plane      --> normal = (a,b,c)
    a = plane[0]
    b = plane[1]
    c = plane[2]
    d = plane[3]

    distToPlane = getDistancePlane(plane, [x,y,z])
    r = np.sqrt(x**2+y**2)
    newy = np.sqrt(r**2-distToPlane**2)

    normal = (1/(np.sqrt(a**2+b**2+c**2)))*np.array([a,b,c])
    newOrigin = distToPlane*normal


    if x >= newOrigin[0] and y >= newOrigin[1]:
        newy = -newy
    if x < newOrigin[0] and y >= newOrigin[1]:
        newy = -newy
    #if newOrigin[0] >= 0 and newOrigin[1] >= 0:
        #r = -np.sqrt(x**2+y**2)
    #if newOrigin[1] >= 0 and newOrigin[0] < 0:
        #r = -np.sqrt(x**2+y**2)
    return newy


