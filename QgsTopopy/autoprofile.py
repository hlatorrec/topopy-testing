# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def get_slope_rval(x, y):
    linreg = linregress(x, y)
    return linreg.slope, linreg.rvalue


def fit_channel(dist, height, side_pts):
    if side_pts < 5:
        print('Too few side points for local fits')
        return 0, 0
    elif side_pts > height.size/2 - 1:
        print("Side points' reach excedes channel lenght")
        return 0, 0
    else:
        slopes = np.array([])
        rvalues = np.array([])
        for i in range(side_pts, height.size - side_pts):
            x = dist[i-side_pts : i+side_pts]
            z = height[i-side_pts : i+side_pts]
            slope, rval = get_slope_rval(x, z)
            slopes = np.append(slopes, slope)
            rvalues = np.append(rvalues, rval)
        
        return slopes, rvalues


def meanmed_compare(dist, height, side_pts_range):
    small_side_pts, big_side_pts = side_pts_range
    means = np.array([])
    medians = np.array([])
    for i in range(small_side_pts, big_side_pts):
        slopes, rvalues = fit_channel(dist, height, i)
        mean = np.mean(rvalues)
        median = np.median(rvalues)
        means = np.append(means, mean)
        medians = np.append(medians, median)
    
    return means, medians
        
###############

channel = np.loadtxt('channel.txt')
z = np.trim_zeros(channel, 'b')
x = np.arange(0, z.size, 1)

#linreg = linregress(x, z)

side_pts = 25
fit_x = x[side_pts : x.size - side_pts]
slopes, rvalues = fit_channel(x, z, side_pts)
rmedian = np.ones(rvalues.size)*np.median(rvalues)
rmean = np.ones(rvalues.size)*np.mean(rvalues)

plt.plot(fit_x, rvalues)
plt.plot(fit_x, rmean)
plt.plot(fit_x, rmedian)

#pts_range = (10, 80)
#means, medians = meanmed_compare(x, z, pts_range)

#meanmed_x = np.arange(10, 80)
#plt.plot(meanmed_x, means)
#plt.plot(meanmed_x, medians)

#diff = means - medians

#plt.plot(meanmed_x, diff)
