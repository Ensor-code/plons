import numpy                    as np
import scipy.ndimage.filters
import bisect


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
    upper_bound = lower_bound + quantity * width
    return [(low, low + width) for low in range(lower_bound, upper_bound + 1, width)]


'''
bins is a list of tuples, like [(0,20), (20, 40), (40, 60)],
binning returns the smallest index i of bins so that
bin[i][0] <= value < bin[i][1]
searches the bin containing the value
'''
def find_bin(value, bins):
    if value < bins[0][0]:
        return 0
    if value > bins[-1][1]:
        return len(bins) - 1
    i = bisect.bisect_left([x[0] for x in bins], value)
    if i == 0:
        return 0
    elif i == len(bins):
        return len(bins) - 1
    else:
        return i - 1

'''
Put the data of a chosen parameter 'data' in bins, according to the parameter 'sortData'.
'''
def bin_data(data,sortData):
    # The parameter 'sortData' sets the number of bins by its maximum
    nb = np.max(sortData)
    # create bins
    number_of_bins = int(np.round( nb + 1 ) )
    bins = create_bins( lower_bound=0, width=1, quantity=number_of_bins )

    # create dictionary with keys: the bin [au] and value: the speed in that bin
    binned_data = {i: [] for i in range(number_of_bins)}  # use dictionary comprehension for creating keys

    for i in range(len(data)):  # make values
        bin_index = find_bin(sortData[i], bins)
        binned_data[bin_index].append(data[i])

    return binned_data, number_of_bins

'''
smoothens values to plot
'''
def smoothen(values,sigma):
    #Sigma chosen random, since trying to set sigma=standard_deviations did not work
    fit=scipy.ndimage.filters.gaussian_filter(values, sigma)
    return fit


'''
Given any number of dictionaries, shallow copy and merge into a new dict,
precedence goes to key value pairs in latter dictionaries.
'''
def merge_dicts(*dict_args):
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


'''
Function to transform the axis labels to multiples of pi/2,
giving 'pi' as a label
'''
def format_func1(value, tick_number):
    # find number of multiples of pi/2
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)


'''
Function to transform the axis labels to multiples of pi/4,
giving 'pi' as a label
'''
def format_func2(value, tick_number):
    # find number of multiples of pi/2
    N = int(np.round(4 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/4$"
    elif N == 2:
        return r"$\pi/2$"
    elif N == 3:
        return r"$3\pi/4$"
    elif N == 4:
        return r"$\pi$"
    elif N % 4 > 0:
        return r"${0}\pi/4$".format(N)
    else:
        return r"${0}\pi$".format(N // 4)
