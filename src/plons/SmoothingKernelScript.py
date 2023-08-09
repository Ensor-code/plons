import healpy               as hp
import numpy                as np
import numba                as nb

from scipy.spatial import cKDTree


'''
Gives the function f which is used by the function 'smoothingKernelPhantom'
From: Price, D. J., Wurster, J., Tricco, T. S., et al. 2018, PASA, 35, e031; Eq. 17
'''
@nb.njit()
def get_f(d):
    f = np.zeros_like(d)
    
    for i in range(len(d)):
        ind1 = (d[i] >= 0) & (d[i] < 1)
        f[i][ind1] = 1-(3/2)*np.square(d[i][ind1])+(3/4)*np.power(d[i][ind1],3)

        ind2 = (d[i] >= 1) & (d[i] < 2)
        f[i][ind2] = (1/4)*np.power(2-d[i][ind2],3)

    return f


Cnorm = 1/np.pi
hfact = 1.2


'''
Smoothing kernel from phantom from the phantom paper Price 2018 (Sect. 2.1.5, 2.1.6)
     r_a is the position [x,y,z] of the point in the healpix sphere
     r_b_list is an array of the [x,y,z] positions of all closest neighbours of a
     h is the smooting length of the point b
RETURNS:
     np.array with the values of the smoothing kernel W_ab for every closest neighbour b 
'''
@nb.njit()
def smoothingKernelPhantom(dist,h_list):
    d = np.abs(dist/h_list)
    W_ab = (Cnorm/hfact**3)*get_f(d)
    return W_ab

'''
Get a parameter in a certain point (not one that exists), using the PHANTOM smoothing kernel
    par = sum(par_i*W_i)
    'param' is the string of the parameters
'''
def getParamSmoothingKernel(closest_points, W_list, param):
    param_list = param[closest_points]
    par =  param_list*W_list
    par_final = np.sum(par, axis = 1)
    return par_final

'''
Get a 2D array containing the [x,y,z]-components of all pixels 
          shape == 'r' --> spherical splice
          shape == 'z' --> planar slice: orbital plane (xy-plane)
          shape == 'x' --> yz-plane
          shape == 'y' --> xz-plane
'''
def getPixels(shape, n, r, data, bound):

    pixCoord = np.zeros(n)
    # for a spherical slice centred around the AGB star (viewpoint from AGB surface) --> use healpy
    if shape == 'r':
        # define the amount of pixels for the spherical plot: always of the form 12*n**2
        npix = 12*n**2
        # get array with 3D positions of the pixels, centred around 0, radius 1
        pix = np.array(hp.pixelfunc.pix2vec(hp.npix2nside(npix), range(npix))).transpose()
        
        # Modify unit sphere to desired sphere
        if r == 'comp':
            radius = np.linalg.norm(np.array(data['AGBcoord'])-np.array(data['compCoord']))-2*data['rAccrComp']
        else:
            r = float(r)
            radius = r 
        # modify the unit spherical slice to a spherical slice centered around the
        #    AGB star, with radius the orbital separation
        pixCoord = radius*pix + np.array(data['AGBcoord'])
    
    # for a planar slice: orbital plane (xy-plane)
    if shape == 'z':
        # get pixels
        pix = np.linspace(-np.abs(bound),np.abs(bound), n)

        pixCoord = []
        for i in range(len(pix)):
            for j in range(len(pix)):
                pixCoord.append([pix[i],pix[j],0])

        pixCoord = np.array(pixCoord)
    
    # xz-plane
    if shape == 'y':
        # get pixels
        pix = np.linspace(-np.abs(bound),np.abs(bound), n)

        pixCoord = []
        for i in range(len(pix)):
            for j in range(len(pix)):
                pixCoord.append([pix[i], 0, pix[j]])

        pixCoord = np.array(pixCoord)
        #print(pixCoord)
    
    # x-line (y=0=z) yz-plane
    if shape == 'line_x':
        pix = np.linspace(-np.abs(bound),np.abs(bound), n)
        # pix_y = np.linspace(-r,r,n)

        pixCoord = []
        for i in range(len(pix)):
            # for j in range(len(pix_y)):
            pixCoord.append([pix[i],0,0])

        pixCoord = np.array(pixCoord)
        
    # y-line (x=0=z)
    if shape == 'line_y':
        pix = np.linspace(-np.abs(bound),np.abs(bound), n)
        # pix_y = np.linspace(-r,r,n)

        pixCoord = []
        for i in range(len(pix)):
            # for j in range(len(pix_y)):
            pixCoord.append([0,pix[i],0])

        pixCoord = np.array(pixCoord)
    
    # z-line (y=posAGB=x)
    if shape == 'line_z':
        pix = np.linspace(-np.abs(bound),np.abs(bound), n)
        # pix_y = np.linspace(-r,r,n)

        pixCoord = []
        for i in range(len(pix)):
            # for j in range(len(pix_y)):
            pixCoord.append([data['posAGB'][0], data['posAGB'][1], pix[i]])

        pixCoord = np.array(pixCoord)

    return pixCoord

'''
Smooth a plane described in the XYZ mesh
      'X', 'Y', 'Z' are meshgrids in the XYZ plane, on which params need to be smoothed
      'dumpData'    is the dictionary with all the data of the dump
      'observables' is a list of strings with the name of the parameter
      'neighbours' is an integer representing how many neighbours are used in the calculation
'''
def smoothMesh(X, Y, Z, dumpData, observables, neighbours = 100):
    pixCoord = convertFromMesh(X, Y, Z)
    smooth = getSmoothingKernelledPix(neighbours, dumpData, observables, pixCoord)
    X, Y, Z, smooth = convertToMesh(pixCoord, smooth, observables)
    return smooth


'''
Get the smoothed values for the params, using the smoothing kernel with neighbours
      'neighbours' is an integer representing how many neighbours are used in the calculation
      'params'  is a list of strings with the name of the parameter
      'data'    is the dictionary with all the data of the dump
'''
def getSmoothingKernelledPix(neighbours, data, params, pixCoord):

    # define all nearest neighbours
    tree = cKDTree(data['position'])

    # for every pixel in the spherical slice (sphere), get its "neighbor" nearest neighbours
    (distances, closest_points) = tree.query(pixCoord, neighbours)

    '''
    get values of certain parameters in the grid for the nearest neighbours
        gives the position (x,y,z) for all the indices of the closest neighbours (closest_closest_points_60)
        Result: 3D array: - per point in the healsphere
                          - of every closest neighbour of that point
                          - the x, y, z values of th position of the closest neighbour
    '''
    # position_closest_points = data['position'][closest_points]
    h_closest_points        = data['h'       ][closest_points]

    # Get the smoothing kernel W_ab for all nearest neighbours for every pixel in 'sphere'
    
    W_ab = smoothingKernelPhantom(distances, h_closest_points)

    results = dict()
    for param in params:
        # Get the corresponding values for the parameter 'param' for every pixel, by param = sum(param_i*W_i), for i ~ [0,neighbours]
        results[param] = getParamSmoothingKernel(closest_points, W_ab, data[param]) 

    return results

def convertFromMesh(X, Y, Z):
    n = len(X)

    r = np.zeros((n**2, 3))

    i0 = 0
    i1 = n
    for i in range(n):
        r[:,0][i0: i1] = np.transpose(X[:n, i])
        r[:,1][i0: i1] = np.transpose(Y[:n, i])
        r[:,2][i0: i1] = np.transpose(Z[:n, i])
        i0 += n
        i1 += n

    return r
            
'''
Transform one dimensional arrays to a 2D mesh. 
    'n' is the length of the arrays
    'r' is the position vector, x=r[:,0], y=r[:,1] and z=r[:,2] the regular one dimensional arrays
    params is a list of strings for the parameters of interest
    data_params is a dictionary for the data of params
'''
def convertToMesh(r, data_params, params):
    n = int(np.sqrt(len(r[:,0])))

    X = np.zeros(shape=(n, n))
    Y = np.zeros(shape=(n, n))
    Z = np.zeros(shape=(n, n))

    i0 = 0
    i1 = n
    for i in range(n):
        X[:n, i] = np.transpose(r[:,0][i0: i1])
        Y[:n, i] = np.transpose(r[:,1][i0: i1])
        Z[:n, i] = np.transpose(r[:,2][i0: i1])
        i0 += n
        i1 += n

    for param in params:
        temp_mesh = np.zeros(shape=(n, n))
        i0 = 0
        i1 = n
        for i in range(n):
            temp_mesh[:n, i] = np.transpose(data_params[param][i0:i1])
            i0 += n
            i1 += n
        data_params[param] = temp_mesh

    return X, Y, Z, data_params

@nb.njit()
def rotatePixCoordAroundZ(theta, pixCoord):
    n = len(pixCoord)
    x = pixCoord[:, 0]
    y = pixCoord[:, 1]
    z = pixCoord[:, 2]

    x_rot = x * np.cos(theta) - y * np.sin(theta)
    y_rot = x * np.sin(theta) + y * np.cos(theta)
    z_rot = z

    rotatedArray = np.zeros(shape=(n, 3))
    rotatedArray[:, 0] = x_rot
    rotatedArray[:, 1] = y_rot
    rotatedArray[:, 2] = z_rot

    return rotatedArray

@nb.njit()
def rotateMeshAroundZ(theta, X, Y, Z):
    X_rot = X * np.cos(theta) - Y * np.sin(theta)
    Y_rot = X * np.sin(theta) + Y * np.cos(theta)
    Z_rot = Z

    return X_rot, Y_rot, Z_rot

def rotateVelocityAroundZ(theta, velocity):
    VX = velocity['vx']
    VY = velocity['vy']
    VZ = velocity['vz']
    
    velocity['vx'] = VX * np.cos(theta) - VY * np.sin(theta)
    velocity['vy'] = VX * np.sin(theta) + VY * np.cos(theta)
    velocity['vz'] = VZ

    return velocity
