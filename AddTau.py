import pandas as pd
import numpy as np 

def addTau(filename):
    ind = filename.rindex('/')+1
    loc = filename[:ind]
    name = filename[ind:]

    header=pd.read_csv(loc+name+'.ascii',sep='&',nrows=12,header=None)
    data = pd.read_csv(loc+name+'.ascii', delim_whitespace=True,skiprows=12,header=None, names=['x', 'y', 'z', 'particle mass', 'h', 'rho',  'temp', 'v_x', 'v_y', 'v_z', 'u', 'alpha', 'div v', 'itype'])
    points = pd.read_csv(loc+'points_'+name+'.txt', delim_whitespace=True,header=None, names=('x','y','z'))
    tauspoints = pd.read_csv(loc+'taus_'+name+'_inwards.txt',header=None, names=['tau'])
    points=points.join(tauspoints)

    list = np.zeros(len(data))
    j=0
    for i in range(len(data)-2):
        if (data['x'][i] - points['x'][j] < 1e-10) and data['y'][i] - points['y'][j] < 1e-10 and data['z'][i] - points['z'][j] < 1e-10:
            list[i] = points['tau'][j]
            j+=1
    data = data.join(pd.DataFrame({'tau':list}))
    header.to_csv(loc+name+'.ascii', sep='\t',index=False,header=False, mode='w')
    data.to_csv(loc+name+'.ascii', sep='\t', index=False, header=False, mode='a')
