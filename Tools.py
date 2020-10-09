import os
import math                 as math
import numpy                as np
import scipy.ndimage.filters




'''
Returns the index [int] of the element in the array of which a[index]-value is minimal.
'''
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



'''
Create_bins returns an equal-width (distance) partitioning. 
    It returns an ascending list of tuples, representing the intervals.
    A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0 
    and i < quantity, satisfies the following conditions:
        (1) bins[i][0] + width == bins[i][1]
        (2) bins[i-1][0] + width == bins[i][0] and
            bins[i-1][1] + width == bins[i][1]
'''
def create_bins(lower_bound, width, quantity):
    bins = []
    for low in range(lower_bound, lower_bound + quantity*width + 1, width):
        bins.append((low, low+width))
    return bins



#""" bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
    #binning returns the smallest index i of bins so that
    #bin[i][0] <= value < bin[i][1]
#"""
# searches the bin containing the value

def find_bin(value, bins):

    if value < bins[0][0]:
        return 0
    elif value > bins[-1][1]:
        return len(bins) - 1
    
    else:
        for i in range(0, len(bins)):
            if bins[i][0] <= value < bins[i][1]:
                return i
    
        return -1    


'''
Put the data of a chosen parameter 'data' in bins, according to the parameter 'sortData'.
'''
def bin_data(data,sortData):
    # The parameter 'sortData' sets the number of bins by its maximum
    nb = max(sortData)
    # create bins
    number_of_bins = int( round( nb + 1) )
    bins = create_bins( lower_bound=0, width=1, quantity=number_of_bins )

    # create dictionary with keys: the bin [au] and value: the speed in that bin
    binned_data = {}               
    for i in range(number_of_bins):             # make keys
        binned_data[i] = []
        
    for i in range(len(data)):                 # make values
        bin_index = find_bin(sortData[i], bins)
        binned_data[bin_index].append(data[i])
        
    return binned_data, number_of_bins


#smoothens values to plot
def smoothen(values,sigma):
    #Sigma chosen random, since trying to set sigma=standard_deviations did not work
    fit=scipy.ndimage.filters.gaussian_filter(values, sigma)
    return fit

    

