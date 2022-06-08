# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def find_crossings(array, value):
    crossings = np.where(np.diff(np.signbit(array-value)))[0] #Para encontrar cambios de signo 
    return crossings


def get_slope_rval(x, y):
    linreg = linregress(x, y)
    return linreg.slope, linreg.rvalue


def fit_channel(dist, height, side_pts):
    if side_pts < 10:
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
        

def find_main_slopes(dist, height):
    side_pts = 9
    gap = 0
    while gap < 0.002:
        side_pts += 1
        slope, rvalues = fit_channel(dist, height, side_pts)
        gap = np.mean(rvalues) - np.median(rvalues)

    crossings = find_crossings(rvalues, np.mean(rvalues))
    dist_index = crossings + side_pts
    #Inserto el cero y el final para los tramos de recta
    dist_index = np.append(dist_index, dist.size-1)
    dist_index = np.insert(dist_index, 0, 0)
    #Faltaría considerar el caso en el que rvalues[0] > rmean
    #Eso solo pasa si side_pts es demasiado pequeño
    #o si el canal no tiene un primer tramo lo sificientemente largo
    sections = []
    i = 0
    while i < 4:
        start = dist_index[i]
        end = dist_index[i+1]
        x = dist[start:end]
        z = height[start:end]
        linreg = linregress(x, z)
        #Basta con guardar pendiente, cte y los puntos para saber hasta dónde usar cada recta
        recta = (linreg.slope, linreg.intercept, dist[start], dist[end])
        sections.append(recta)
        i+=2
        
    return sections
        
    
    
    

###############

channel = np.loadtxt('channel.txt')
z = np.trim_zeros(channel, 'b')
x = np.arange(0, z.size, 1)

rectas = find_main_slopes(x, z)
der = np.gradient(z, x)
plt.plot(x,der)
for recta in rectas:
    rx = x[recta[2] : recta[3]]
    ry = recta[0]*np.ones(rx.size)
    plt.plot(rx,ry)

#linreg = linregress(x, z)

'''
side_pts = 10
fit_x = x[side_pts : x.size - side_pts]
slopes, rvalues = fit_channel(x, z, side_pts)
rmedian = np.ones(rvalues.size)*np.median(rvalues)
rmean = np.ones(rvalues.size)*np.mean(rvalues)

gap = np.mean(rvalues) - np.median(rvalues)

plt.plot(fit_x, rvalues)
plt.plot(fit_x, rmean)
plt.plot(fit_x, rmedian)
crossings = find_crossings(rvalues, rmedian) #Para encontrar cambios de signo
'''

#pts_range = (10, 80)
#means, medians = meanmed_compare(x, z, pts_range)

#meanmed_x = np.arange(10, 80)
#plt.plot(meanmed_x, means)
#plt.plot(meanmed_x, medians)

#diff = means - medians

#plt.plot(meanmed_x, diff)
