import numpy as np
import matplotlib.pyplot as plt
import random as rd
from mpl_toolkits import mplot3d


'''

Density and cross-section function dependent of the distance "distance".

'''

def density(distance):
    return 1

def cross_section(distance):
    #return np.exp(-distance)
    return 0.01 #-distance**3

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
        self.test = []

        for i in range(self.number_of_scatterings):
            phi   = 2*np.pi*rd.uniform(0,1)
            theta =   np.pi*rd.uniform(0,1)
            self.angle_list.append([phi,theta])
            self.tau.append(-np.log(1-rd.uniform(0,1)))

        vx0=np.sin(self.angle_list[0][1])*np.cos(self.angle_list[0][0])
        vy0=np.sin(self.angle_list[0][1])*np.sin(self.angle_list[0][0])
        vz0=np.cos(self.angle_list[0][1])
        self.list_of_points.append([self.starting_point, [vx0,vy0,vz0]])

        for j in range(self.number_of_scatterings):
            self.starting_point = [self.starting_point[0]+calculate_L(self.starting_point,self.tau[j])*np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0]),
                                   self.starting_point[1]+calculate_L(self.starting_point,self.tau[j])*np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0]),
                                   self.starting_point[2]+calculate_L(self.starting_point,self.tau[j])*np.cos(self.angle_list[j][1])]
            vx=np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0])
            vy=np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0])
            vz=np.cos(self.angle_list[j][1])
            self.vxvyvz = [vx,vy,vz]
            self.list_of_points.append([self.starting_point, [vx, vy, vz]])
    def plot_3D(self):
        '''
        You'll need to add :
        fig = plt.figure()
        ax  = plt.axes(projection='3d')

        onto the code before calling plot_3D.
        '''
        ax.plot3D([self.list_of_points[w][0][0] for w in range(len(self.list_of_points))], [self.list_of_points[w][0][1] for w in range(len(self.list_of_points))], [self.list_of_points[w][0][2] for w in range(len(self.list_of_points))])


if __name__ == '__main__':
    ax = plt.axes()
    ax.set_facecolor("black")
    distance_from_telescope = 100
    for i in range(10000):
        photon1 = photon(5, [0.1,0.1,0.1])
        print(' '+str(100*i/(10000))+'%', end="\r")
        if photon1.vxvyvz[0]!=0 and photon1.vxvyvz[1]!=0 and photon1.vxvyvz[2]>0 :
            ax.scatter(photon1.starting_point[0],photon1.starting_point[1],s=0.3,color='white',alpha=0.1)
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    plt.savefig('prout.png',dpi=300)
    plt.show()
