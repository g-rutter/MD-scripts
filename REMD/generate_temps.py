'''Little script to return a list of temperatures of replicas, based on interpolating a function you give'''
from numpy import linspace, round, zeros, interp
import numpy

#settings
n_replicas=30
samples=100000
x_points = [275, 288, 305, 322, 350]
density_points = [0.5,1.0,3.0,1.0,0.5]

linear_temps = linspace(275,350,samples)
density = interp( linear_temps, x_points, density_points )
density_norm = (density*(n_replicas-1))/density.sum()

nonlinear_temps=zeros(n_replicas)
count=1.0
j=0

for i in range(samples):

   if count >= 1.0:
      nonlinear_temps[j]=linear_temps[i]
      count = 0.0
      j += 1

   count+=density_norm[i]

nonlinear_temps[-1]=linear_temps[-1]

for temp in round(nonlinear_temps, 2):
   print temp
