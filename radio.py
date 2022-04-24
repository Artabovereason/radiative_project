import numpy as np
import matplotlib.pyplot as plt
import random as rd
import timeit
from mpl_toolkits import mplot3d
start_time = timeit.default_timer()

def wavelength_to_rgb(wavelength):
    '''This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    gamma=0.8
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B)

def planck_law(frequency,temperature):
    return 2*(frequency**3)/(np.exp(frequency/temperature)-1)

'''
    Photon is the class for the photon source,
    number_of_scatterings : is the number of scattering that will be computed.
    starting_point        : takes a [x,y,z] array that define the starting point of the 1st scattering.
    frequency             : defines the frequency of the photon, 'random' gives a random color.
                            'random-visible'    : uniform distribution of random in the visible spectrum.
                            'random-planck'     : distribution of random following Planck's distribution.
                            'random-gaussian'   : distribution of random following a gaussian distribution.
                            'random-lorentzian' : distribution of random following a lorentzian distribution.
'''
class photon:
    def __init__(self,number_of_scatterings,starting_point,frequency):
        self.frequency             = frequency
        if   self.frequency == 'random-visible':
            self.frequency         = rd.uniform(380,750)
        elif self.frequency == 'random-planck':
            self.frequency         = 0
        elif self.frequency == 'random-gaussian':
            self.frequency         = rd.gauss(500, 100)
        elif self.frequency == 'random-lorentzian':
            self.frequency         = 0
        else:
            print('This frequency generation is not implemented yet')
        self.number_of_scatterings = number_of_scatterings
        self.starting_point        = starting_point
        '''
        list_of_points   : list of all scattering points.
        angle_list       : list of all the trajectory angles (theta,phi).
        tau              : optical depth.
        angle_trajectory : defines the angle of the trajectory (vx,vy) in the (x,y) plane.
        '''
        self.list_of_points        = []
        self.angle_list            = []
        self.tau                   = []
        self.angle_trajectory      = 0
        self.vxvyvz=[0,0,1]
        for j in range(self.number_of_scatterings):
            self.tau.append(-np.log(1-rd.uniform(0,1)))
            '''
            This allows us to add an absorption probability to occur.
            We make the [vx,vy,vz] (i.e. the impulsion in each direction) equal to 0
            and the trajectory angle to pi/2 so that the condition that the photon
            hit the telescope is never satisfied.
            '''
            point=self.starting_point
            x=point[0]
            y=point[1]
            z=point[2] #extracts the starting point of the photon
            if (x<0.5 and x>-0.5 and ((y>-5 and y<-0.5) or (y>0.5 and y<5)) and z>1 and z<2):
                self.absorption=0.5 # If the photon goes inside the absorbing medium (finger) it may be absorbed with probability 0.5
            elif(x<1.5 and x>-1.5 and y<0.5 and y>-0.5 and z>1 and z<2):
                self.absorption=1 #If the photon goes inside the ring it is absorbed
            else:
                self.absorption=0 #Otherwise it is not absorbed
            absorption_random = rd.uniform(0,1)
            if absorption_random >1-self.absorption: #Absorption probability
                self.vxvyvz           = [0,0,0]
                self.angle_trajectory = np.pi/2 #We remove the photon in this case
                break
            else:
                if(int(physics.density(point, self.frequency)*physics.cross_section(point, self.frequency))!=0):
                    self.starting_point   = [self.starting_point[0]+physics.calculate_L(self.starting_point,self.tau[j],self.vxvyvz, self.frequency)*np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0]), self.starting_point[1]+physics.calculate_L(self.starting_point,self.tau[j],self.vxvyvz, self.frequency)*np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0]),self.starting_point[2]+physics.calculate_L(self.starting_point,self.tau[j],self.vxvyvz, self.frequency)*np.cos(self.angle_list[j][1])] #New point for the photon

                    self.vxvyvz           = [np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0]), np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0]), np.cos(self.angle_list[j][1])] #New direction

                    self.list_of_points.append([self.starting_point, self.vxvyvz])
                    self.angle_trajectory = np.arctan( float(self.vxvyvz[0])/float(self.vxvyvz[1]) )

                else:
                    self.starting_point=[self.starting_point[0]+self.vxvyvz[0]*self.tau[j-1], self.starting_point[1]+self.vxvyvz[1]*self.tau[j-1], self.starting_point[2]+self.vxvyvz[2]*self.tau[j-1]]
            phi   = 2*np.pi*rd.uniform(0,1)
            theta =   np.pi*rd.uniform(0,1)
            self.angle_list.append([phi,theta])


    def plot_3D(self):
        '''
        You'll need to add :
        fig = plt.figure()
        ax  = plt.axes(projection='3d')

        onto the code before calling plot_3D.
        '''
        ax.plot3D([self.list_of_points[w][0][0] for w in range(len(self.list_of_points))], [self.list_of_points[w][0][1] for w in range(len(self.list_of_points))], [self.list_of_points[w][0][2] for w in range(len(self.list_of_points))])

class physics_geometry:
    def __init__(self, R, d, x ):
        self.self = self
        self.planet_radius=R
        self.planet_distance=d
        self.planet_position=x
    '''
    Density and cross-section function dependent of the distance 'distance' and the wavelength 'wavelength'
    '''
    def density(self, point, c):
        if (x<0.5 and x>-0.5 and ((y>-5 and y<-0.5) or (y>0.5 and y<5)) and z>1 and z<2):
            return 1000 #High density in the finger
        elif(x<1.5 and x>-1.5 and y<0.5 and y>-0.5 and z>1 and z<2):
            return 10000 #Very high density in the ring
        else:
            return 0 #No scattering outside the finger


    def cross_section(self, point, c):
    #return np.exp(-distance)
        if (x<1 and x>-1 and ((y>-5 and y<-0.5) or (y>0.5 and y<5)) and z>1 and z<2):
            return 1 #Non zero cross section in the finger
        elif(x<1.5 and x>-1.5 and y<0.5 and y>-0.5 and z>1 and z<2):
            return 10000 #High cross section in the ring
        else:
            return 0

    '''Function to calculate the length between two points, knowning the expected optical depth tau and the starting point "point" in [x,y,z] format'''

    def calculate_L(self, point,tau, dir, c):
        tau_s = 0
        length_point = np.sqrt(point[0]**2+point[1]**2+point[2]**2)
        delta_l =  tau/(self.density(point, c)*self.cross_section(point, c))/100
        increment = 0
        point1=[0,0,0]
        point1[0]=point[0]+dir[0]*delta_l
        point1[1]=point[1]+dir[1]*delta_l
        point1[2]=point[2]+dir[2]*delta_l
        while tau_s < tau :
            tau_s+=tau_s +(self.density(point1, c)+self.density(point, c) )*delta_l/2
            increment += 1
        return increment*delta_l



if __name__ == '__main__':
    '''
    ax.set_facecolor("black")  : set the background of the plot to black.
    length_telescope_to_object : length from the telescope opening to the source
    number_photons             : number of photons we wish to generate in our simulation
    temperature_of_object      : this is used only in the Planck's distribution for the frequency generation in photon class.
    number_sources             : number of point source that will be considered.
    center_sources             : each source has a defined center given by the user.
    physics                    : define the density (density), cross section (cross_section) and the next scattering point calculation (calculate_L).
    '''
    ax = plt.axes()
    ax.set_facecolor("black")

    length_telescope_to_object = 10
    #number_photons             = 50000
    temperature_of_object      = 1000
    number_sources             = 1
    center_sources             = [[0,0,0]]
    Ns=[]
    R_star=2
    for n in range(1):
        x_p=-10+n*2
        physics                    = physics_geometry(0.5, 10, 0)
        for j in range(number_sources):
            N_photons=0
            for i in range(number_photons):
                x=rd.uniform(-10,10)
                y=rd.uniform(-10,10)
                z=0
                photon1 = photon(50, [x,y,z],'random-gaussian')
                print(' '+str(100*i/(number_photons))+'%', end="\r")                #Shows the % of completion of the code
                point=photon1.starting_point
                if photon1.vxvyvz[2]>0 :
                    direction=photon1.vxvyvz
                    x=direction[0]*((length_telescope_to_object-point[2])/direction[2])+point[0]
                    y=direction[1]*((length_telescope_to_object-point[2])/direction[2])+point[1] #Point of intersection of trajectory and plane of the telescope
                    if(x<10 and x>-10 and y > -10 and y<10) :
                        ax.scatter(x,y,s=0.3,color="white",alpha=0.1)#color='white',alpha=0.1) #If photon goes inside the telescope: detected
                        N_photons+=1
                    else :
                        pass
        Ns.append(N_photons)
        print(Ns)

    plt.xlim(-10,10)
    plt.ylim(-10,10)
    plt.savefig('images/radio_detection.png',dpi=300)
    #plt.show()

    stop_time = timeit.default_timer()
    print(' ')
    print('time of simulation (in s): ', stop_time - start_time)
    print(' ')
