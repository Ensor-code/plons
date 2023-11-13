import numpy                    as np
import plons
import os
import k3d
from ipywidgets import widgets,interact


#import necessary plons scripts
import plons.ConversionFactors_cgs    as cgs
import plons.PhysicalQuantities       as pq
import plons.GeometricalFunctions     as gf


import warnings
warnings.filterwarnings("ignore")


'''
Makes 3D visualisation of all pts, with colorbar for parameter 'par' with colorbar limits 'colMin' and 'colMax'

pts is a 3D array with the x,y,z locations
you can alter the point size 'pt_size' and opacity 'opty'to create a better visualisation according to your needs and the size of the model
(These parameters, together with the colorbar limits can also be adapted interactively in the plot)
'''
def make3DPlot(par,pts,colMin,colMax,pt_size =2, opty = 0.8 ):
    # Rescale to log scale to simplify plotting
    lpar = np.log10(par)

    # Visualise the points
    plot = k3d.plot()#name=title)
    plt_points = k3d.points(positions=pts.astype(np.float32),       # matrix met xyz van alles
                            colors=lpar.astype(np.float32),
                            attribute=lpar.astype(np.float32),        # om opaciteit te berekenen, log om mooi te schalen
    #                         opacity_function=[minld,maxld],
                            color_range=[colMin,colMax],
                            point_size=pt_size,
                            shader='3d',
                            opacity=opty)

    basic_color_maps = [attr for attr in dir(k3d.basic_color_maps) if not attr.startswith('__')]
    @interact(x=widgets.Dropdown(options=basic_color_maps, value=basic_color_maps[0], description='Basic ColorMap:'))
    def g(x):
        plt_points.color_map = getattr(k3d.basic_color_maps, x)

    plot += plt_points
    plot.display()

'''
Makes 3D visualisation of all points for which a 1D velocity v_dir is limited in absolute value by vmin and vmax
This is usefull to understand which data is used to calculate certain velocity channel maps in radiative transfer postprocessing

there is colorbar for 'par' values, with limits 'colMin' and 'colMax'
points is a 3D array with the x,y,z locations
you can alter the point size 'pt_size' and opacity 'opty'to create a better visualisation according to your needs and the size of the model
(These parameters, together with the colorbar limits can also be adapted interactively in the plot)
if you give values for rmax and rmin, only points within this radius range will be plot
'''

def makePlot_vDir(v_dir, r, par, pts, vmin, vmax,rmax = 0,rmin = 0,pt_size = 2,opty=0.8):#, title):

    par = par [np.abs(v_dir) < vmax]
    pts = pts [np.abs(v_dir) < vmax]
    r   = r   [np.abs(v_dir) < vmax]
    v_dir  = v_dir  [np.abs(v_dir) < vmax]
    
    
    par = par [np.abs(v_dir) > vmin]
    pts = pts [np.abs(v_dir) > vmin]
    r   = r   [np.abs(v_dir) > vmin]
    v_dir  = v_dir  [np.abs(v_dir) > vmin]
    
    if rmax > 0:
        par = par   [r < rmax]
        pts = pts   [r < rmax]
        v_dir  = v_dir    [r < rmax]   
        r   = r     [r < rmax]
    if rmin > 0:
        par = par   [r > rmin]
        pts = pts   [r > rmin]
        v_dir  = v_dir    [r > rmin]   
        r   = r     [r > rmin]
   
    lpar   = np.log10(par)

    # Visualise the points
    plot = k3d.plot()#name=title)
    plt_points = k3d.points(positions=pts.astype(np.float32),       # matrix met xyz van alles
                            colors=lpar.astype(np.float32),
                            attribute=lpar.astype(np.float32),        # om opaciteit te berekenen, log om mooi te schalen
    #                         opacity_function=[minld,maxld],
                            # color_range=[colMin,colMax],
                            point_size=pt_size,
                            shader='3d',
                            opacity=opty)

    basic_color_maps = [attr for attr in dir(k3d.basic_color_maps) if not attr.startswith('__')]
    @interact(x=widgets.Dropdown(options=basic_color_maps, value=basic_color_maps[0], description='Basic ColorMap:'))
    def g(x):
        plt_points.color_map = getattr(k3d.basic_color_maps, x)

    plot += plt_points
    plot.display()

