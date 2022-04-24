import numpy as np
import matplotlib.pyplot as plt
import random as rd
from mpl_toolkits import mplot3d
import timeit

'''
This code is the check that the generation is well made.
This code only implemented a point-like generation, but it could be generalised.
There is also only few photons for this check, in order to see something in the
output figure.
'''

start_time = timeit.default_timer()

'''

Density and cross-section function dependent of the distance "distance".

'''
def density(distance):
    return 1
def cross_section(distance):
    return 1
'''
Function to calculate the length between two points, knowning the expected optical depth tau and the starting point "point" in [x,y,z] format
'''
def calculate_L(point,tau):
    tau_s = 0
    length_point = np.sqrt(point[0]**2+point[1]**2+point[2]**2)
    delta_l =  tau/(density(length_point)*cross_section(length_point))/100
    increment = 0
    while tau_s < tau :
        tau_s+=tau_s +(density(length_point+delta_l)+density(length_point) )*delta_l/2
        increment += 1
    return increment*delta_l

class photon:
    def __init__(self,number_of_scatterings,starting_point):
        self.number_of_scatterings = number_of_scatterings
        self.starting_point        = starting_point
        self.list_of_points        = []
        self.angle_list            = []
        self.tau                   = []
        self.list_of_points.append(starting_point)

        for i in range(self.number_of_scatterings):
            phi   = 2*np.pi*rd.uniform(0,1)
            theta =   np.pi*rd.uniform(0,1)
            self.angle_list.append([phi,theta])
            self.tau.append(-np.log(1-rd.uniform(0,1)))

        for j in range(self.number_of_scatterings):
            self.starting_point = [self.starting_point[0]+calculate_L(self.starting_point,self.tau[j])*np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0]), self.starting_point[1]+calculate_L(self.starting_point,self.tau[j])*np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0]), self.starting_point[2]+calculate_L(self.starting_point,self.tau[j])*np.cos(self.angle_list[j][1])]
            self.list_of_points.append(self.starting_point)

    def plot_3D(self):
        '''
        You'll need to add :
        fig = plt.figure()
        ax  = plt.axes(projection='3d')

        onto the code before calling plot_3D.
        '''
        ax.plot3D([self.list_of_points[w][0] for w in range(len(self.list_of_points))], [self.list_of_points[w][1] for w in range(len(self.list_of_points))], [self.list_of_points[w][2] for w in range(len(self.list_of_points))])

if __name__ == '__main__':
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
    number_test    = 5    #number of allowed scattering before terminating the simulation
    number_photons = 100  #number of photons generated
    for i in range(number_photons):
        print(' '+str(100*i/(number_photons))+'%', end="\r")
        photon1     = photon(number_test,[0,0,0])
        norm_photon = np.sqrt(photon1.list_of_points[-1][0]**2+photon1.list_of_points[-1][1]**2+photon1.list_of_points[-1][2]**2)
        photon1.plot_3D()
    ax.set_xlim(-0.5,0.5)
    ax.set_ylim(-0.5,0.5)
    ax.set_zlim(-0.5,0.5)
    stop_time = timeit.default_timer()
    print(' ')
    print('time of simulation (in s): ', stop_time - start_time)
    print(' ')
    plt.savefig('images/3Dplot.png',dpi=600)
