import healpy               as hp
import numpy                as np

from scipy.spatial import cKDTree


'''
Gives the function f which is used by the function 'smoothingKernelPhantom'
From: Price, D. J., Wurster, J., Tricco, T. S., et al. 2018, PASA, 35, e031; Eq. 17
'''
def get_f(d):
    if d >= 0 and d < 1:
        f = 1-(3/2)*d**2+(3/4)*d**3
    if d >= 1 and d < 2:
        f = (1/4)*(2-d)**3
    if d >= 2:
        f = 0
    return f

'''
Smoothing kernel from phantom from the phantom paper Price 2018 (Sect. 2.1.5, 2.1.6)
     r_a is the position [x,y,z] of the point in the healpix sphere
     r_b_list is an array of the [x,y,z] positions of all closest neighbours of a
     h is the smooting length of the point b
RETURNS:
     np.array with the values of the smoothing kernel W_ab for every closest neighbour b 
'''
def smoothingKernelPhantom(r_a,r_b_list,h_list):
    Cnorm = 1/np.pi
    W_ab = []
    for i in range(len(h_list)):
        d = np.abs(np.linalg.norm(r_a-r_b_list[i])/h_list[i])
        W_ab.append((Cnorm/h_list[i]**3)*get_f(d))
    return W_ab

'''
Get a parameter in a certain point (not one that exists), using the PHANTOM smoothing kernel
    par = sum(par_i*W_i)
    'param' is the string of the parameters
'''
def getParamSmoothingKernel(closest_points, W_list,param, massPerRho):
    param_list = param[closest_points]
    par =  param_list*W_list*massPerRho
    par_final = []
    for i in range(len(par)):
        par_final.append(sum(par[i]))
    return np.array(par_final)

'''
Get a 2D array containing the [x,y,z]-components of all pixels 
          shape == 'r' --> spherical splice
          shape == 'z' --> planar slice: orbital plane (xy-plane)
          shape == 'x' --> yz-plane
          shape == 'y' --> xz-plane
'''
def getPixels(shape, n, r, data, bound):
    
    # for a spherical slice --> use healpy
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
        pix = np.linspace(-bound,bound, n)

        pixCoord = []
        for i in range(len(pix)):
            for j in range(len(pix)):
                pixCoord.append([pix[i],pix[j],0])

        pixCoord = np.array(pixCoord)
        #print(pixCoord)
    
    # xz-plane
    if shape == 'y':
        # get pixels
        pix = np.linspace(-bound,bound, n)

        pixCoord = []
        for i in range(len(pix)):
            for j in range(len(pix)):
                pixCoord.append([pix[i],0,pix[j]])

        pixCoord = np.array(pixCoord)
        #print(pixCoord)
    
    # x-line (y=0=z) yz-plane
    if shape == 'line_x':
        pix = np.linspace(-bound,bound, n)
        # pix_y = np.linspace(-r,r,n)

        pixCoord = []
        for i in range(len(pix)):
            # for j in range(len(pix_y)):
            pixCoord.append([pix[i],0,0])

        pixCoord = np.array(pixCoord)
        
    # y-line (x=0=z)
    if shape == 'line_y':
        pix = np.linspace(-bound,bound, n)
        # pix_y = np.linspace(-r,r,n)

        pixCoord = []
        for i in range(len(pix)):
            # for j in range(len(pix_y)):
            pixCoord.append([0,pix[i],0])

        pixCoord = np.array(pixCoord)
    
    # z-line (y=0=x7)
    if shape == 'line_z':
        pix = np.linspace(-bound,bound, n)
        # pix_y = np.linspace(-r,r,n)

        pixCoord = []
        for i in range(len(pix)):
            # for j in range(len(pix_y)):
            pixCoord.append([0,0,pix[i]])

        pixCoord = np.array(pixCoord)

    return pixCoord


'''
Get the values on a healpy sphere for n pixels, 'neighbours' neighbours for a list of parameters 'params'
      'n' and 'neighbours' are integers
      'params'  is a list of strings with the name of the parameter
      'data'    is the dictionary with all the data of the dump
      'r'       is a string with the radius of the wanted sphere
                if 'comp', take the radius up to the companion, else, float('r') --> radius
      'shape' to specify which kind of slice
'''
def getSmoothingKernelledPix(n, neighbours, data, params, r, shape, bound):
    
    rho = data['rho']
    mass = data['mass']
    
    mPrho = mass/rho
    
    pixCoord = getPixels(shape, n, r, data, bound)
    
    # define all nearest neighbours
    tree = cKDTree(data['position'])    
    
    # for every pixel in the spherical slice (sphere), get its 60 nearest neighbours
    (distances, closest_points) = tree.query(pixCoord, neighbours)
    
    '''
    get values of certain parameters in the grid for the nearest neighbours
        gives the position (x,y,z) for all the indices of the closest neighbours (closest_closest_points_60)
        Result: 3D array: - per point in the healsphere
                          - of every closest neighbour of that point
                          - the x, y, z values of th position of the closest neighbour
    '''
    position_closest_points = data['position'][closest_points]
    h_closest_points        = data['h'       ][closest_points]
    mPrho_closest_points    = mPrho[closest_points]
    
    # Get the smoothing kernel W_ab for all nearest neighbours for every pixel in 'sphere'
    W_ab = []
    for i in range(len(pixCoord)):
        W_ab.append(smoothingKernelPhantom(pixCoord[i],position_closest_points[i],h_closest_points[i]))
        
    results = {}
    for i in range(len(params)):
        # Get the corresponding values for the parameter 'param' for every pixel, by param = sum(param_i*W_i), for i ~ [0,neighbours]
        results[str(params[i])] = getParamSmoothingKernel(closest_points, W_ab, data[params[i]], mPrho_closest_points) 

    if shape == 'r':
        return results
    
    if shape == 'z' or shape == 'y' or shape == 'line_x' or shape == 'line_y' or shape == 'line_z':
        pixCoord = np.array(pixCoord).transpose()
        return results, pixCoord[0], pixCoord[1], pixCoord[2]
    
    
    
    
    
