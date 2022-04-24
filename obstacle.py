import numpy as np
import matplotlib.pyplot as plt
import random as rd
import timeit

import scipy.stats
from scipy.stats import rv_continuous


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

class planck_gen(rv_continuous):
    "Planck Law distribution distribution"
    def _pdf(self, x):
        temperature = 6000
        h           = 6.62607004 * 10**(-34)
        c           = 299792458
        kb          = 1.38064852 * 10**(-23)
        if x <0:
            return 0
        else:
            return (1/(5.670400*10**(-8)*temperature**4/np.pi))*(2*h*c**2/(x**5))*(np.exp(h*c/(x*kb*temperature))-1)**(-1)


'''
    Photon is the class for the photon source,
    number_of_scatterings : is the number of scattering that will be computed.
    starting_point        : takes a [x,y,z] array that define the starting point of the 1st scattering.
    frequency             : defines the frequency of the photon, 'random' gives a random color.
                            'random-visible'    : uniform distribution of random in the visible spectrum.
                            'random-planck'     : distribution of random following Planck's distribution at a given temperature.
                            'random-gaussian'   : distribution of random following a gaussian distribution centered on a given frequency.
                            'monochromatic'     : all photons are of the same wavelength
'''
class photon:
    def __init__(self,number_of_scatterings,starting_point,frequency):
        self.frequency     = frequency
        if   self.frequency == 'random-visible':
            self.frequency = rd.uniform(380,750)
        elif self.frequency == 'random-planck':
            planck         = planck_gen(a=380*10**(-9),b=750*10**(-9))
            self.frequency = 10**9*planck.rvs()
        elif self.frequency == 'random-gaussian':
            gauss          = rd.gauss(500, 100)
            while gauss >750 or gauss <380:
                gauss      = rd.gauss(500, 100)
            self.frequency = gauss
        elif self.frequency == 'monochromatic':
            self.frequency = 700
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

        for j in range(self.number_of_scatterings):
            phi   = rd.uniform(0,2*np.pi)
            theta = rd.uniform(0,  np.pi)
            self.angle_list.append([phi,theta])
            self.tau.append(-np.log(1-rd.uniform(0,1)))
            '''
            This allows us to add an absorption probability to occur.
            We make the [vx,vy,vz] (i.e. the impulsion in each direction) equal to 0
            and the trajectory angle to pi/2 so that the condition that the photon
            hit the telescope is never satisfied.
            '''
            absorption_random = rd.uniform(0,1)
            if absorption_random >0.9: #Absorption probability
                self.vxvyvz           = [0,0,0]
                self.angle_trajectory = np.pi/2
                break
            else:
                self.starting_point   = [self.starting_point[0]+physics.calculate_L(self.starting_point,self.tau[j],self.frequency)*np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0]),
                                         self.starting_point[1]+physics.calculate_L(self.starting_point,self.tau[j],self.frequency)*np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0]),
                                         self.starting_point[2]+physics.calculate_L(self.starting_point,self.tau[j],self.frequency)*np.cos(self.angle_list[j][1])]
                self.vxvyvz           = [np.sin(self.angle_list[j][1])*np.cos(self.angle_list[j][0]),
                                         np.sin(self.angle_list[j][1])*np.sin(self.angle_list[j][0]),
                                         np.cos(self.angle_list[j][1])]
                self.list_of_points.append([self.starting_point, self.vxvyvz])
                self.angle_trajectory = np.arctan( float(self.vxvyvz[0])/float(self.vxvyvz[1]) )
                if write_f == 'true':
                    file.write(str(self.list_of_points[j][0][0])+str(' ')+str(self.list_of_points[j][0][1])+str(' ')+str(self.list_of_points[j][0][2])+str(' '))
                    file.write(str(self.list_of_points[j][1][0])+str(' ')+str(self.list_of_points[j][1][1])+str(' ')+str(self.list_of_points[j][1][2])+str(' '))
                    file.write(str(self.angle_list[j][0])+str(' ')+str(self.angle_list[j][1]))
                    file.write('\n')
            '''
            if obstacle_we.name == 'planet':
                abs = obstacle_we.planet(self.starting_point, self.vxvyvz)
                if abs==1:
                    self.vxvyvz = [0,0,-1]
            elif obstacle_we.name == 'absorbing_area':
                abs = obstacle_we.absorbing_area(self.starting_point, self.vxvyvz)
                if abs==1:

                    self.vxvyvz = [0,0,-1]
                    self.frequency = 400
                    break
            '''


    def plot_3D(self):
        '''
        You'll need to add :
        fig = plt.figure()
        ax  = plt.axes(projection='3d')

        onto the code before calling plot_3D.
        '''
        ax.plot3D([self.list_of_points[w][0][0] for w in range(len(self.list_of_points))], [self.list_of_points[w][0][1] for w in range(len(self.list_of_points))], [self.list_of_points[w][0][2] for w in range(len(self.list_of_points))])

class physics_geometry:
    def __init__(self):
        self.self = self
    '''
    Density and cross-section function dependent of the distance 'distance' and the wavelength 'wavelength'
    '''
    def density(self,distance,wavelength):
        return 1
    def cross_section(self,distance,wavelength):
        return 0.1
    '''
    Function to calculate the length between two points, knowning the expected optical depth tau and the starting point 'point' in [x,y,z] format
    '''
    def calculate_L(self,point,tau,wavelength):
        tau_s          = 0
        length_point   = np.sqrt(point[0]**2+point[1]**2+point[2]**2)
        delta_l        = tau/(self.density(length_point,wavelength)*self.cross_section(length_point,wavelength))/100
        increment      = 0
        while tau_s < tau :
            tau_s     += tau_s +(self.density(length_point+delta_l,wavelength)+self.density(length_point,wavelength) )*delta_l/2
            increment += 1
        return increment*delta_l

class obstacle:
    def __init__(self,name,size_x,size_y,distance,position):
        self.name            = name
        self.planet_radius   = 2*size_x
        self.planet_distance = distance
        self.planet_position = position

    def planet(self, point, vec_traj):
        z = self.planet_distance
        if vec_traj[2]>0:
            x = vec_traj[0]*((self.planet_distance-point[2])/vec_traj[2])+point[0]
            y = vec_traj[1]*((self.planet_distance-point[2])/vec_traj[2])+point[1]
            if (np.sqrt((x-self.planet_position)**2+y**2)<=self.planet_radius):
                return 1
            else:
                return 0
    def absorbing_area(self,point,vec_traj):
        if np.abs(point[0])<2 and np.abs(point[1])<2 and np.abs(point[2])>4:
            return 1

if __name__ == '__main__':
    '''
    ax.set_facecolor("black")  : set the background of the plot to black.
    length_telescope_to_object : length from the telescope opening to the source
    number_photons             : number of photons we wish to generate in our simulation
    number_sources             : number of point source that will be considered.
    center_sources             : each source has a defined center given by the user.
    physics                    : define the density (density), cross section (cross_section) and the next scattering point calculation (calculate_L).
    user_mode                  : 'code' or 'terminal'
                                 'code mean' that all the preliminary information are given to the code on the code,
                                 'terminal' mean that all the preliminary information are asked to the user upon running the python script.
    write_f                    : 'true' or 'false' if true, save data on a .txt
    size_source                : this parameter occur only if the source is not 'perfect-point', and thus the generation of the photons occur in.
    dimension_telescope        : the captation zone is a square of (x,y) dimension (dimension_telescope_x,dimension_telescope_y)
    number_capted_photons      : allow us to then compute the number of photons in capted in the telescope.
    '''
    ax = plt.axes()
    ax.set_facecolor("black")
    user_mode = 'code'
    physics   = physics_geometry()
    write_f   = 'false'
    if write_f == 'true':
        file      = open('structure.txt', 'w')
    '''
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    ax.set_zlim(-2,2)
    '''
    number_scatterings    = 5
    number_capted_photons = 0
    obstacle_we = obstacle('absorbing_area',1,0,10,0)
    if   user_mode == 'code'    :
        length_telescope_to_object = 100
        number_photons             = 100000
        number_sources             = 1
        center_sources             = [[0,0,0]]
        generation_type            = 'perfect-point'
        if generation_type != 'perfect-point':
            size_source_x          = 1
            size_source_y          = 5
            size_source_z          = 0.1
        generation_frequency       = 'monochromatic'
        dimension_telescope_x      = 20
        dimension_telescope_y      = 20
    elif user_mode == 'terminal':
        length_telescope_to_object = float(input('What\'s the length between the telescope and the source ? (in meters)'))
        number_photons             = float(input('Number of photons in source ?'))
        number_sources             = float(input('Number of sources'))
        center_sources             = []
        for i in range(number_sources):
            cache                  = input('Center of the source number '+str(i)+' in the format [x,y,z]')
            center_sources.append(cache)
        generation_type            = input('Type of geometry of the source generation ? (perfect-point, perfect-spherical, perfect-cubic)')
        if generation_type != 'perfect-point':
            size_source_x = input('What is the size of your source in x ?')
            size_source_x = input('What is the size of your source in y ?')
            size_source_z = input('What is the size of your source in z ?')
        generation_frequency       = input('Type of frequency generation ? (random-visible, random-planck, random-gaussian, monochromatic)')
        dimension_telescope_x      = input('Size of the input telescope in x? (base is 20)')
        dimension_telescope_y      = input('Size of the input telescope in y? (base is 20)')
    else:
        print('The user mode you requested is not implemented yet.')
    if write_f == 'true':
        file.write('START Parameters\n')
        file.write(str('length_telescope_to_object = ')+str(length_telescope_to_object)+str('\n'))
        file.write(str('number_photons = ')      +str(number_photons)          +str('\n'))
        file.write(str('number_sources = ')      +str(number_sources)          +str('\n'))
        file.write(str('generation_type = ')     +str(generation_type)         +str('\n'))
        file.write(str('generation_frequency = ')+str(generation_frequency)    +str('\n'))
        if generation_type != 'perfect-point':
            file.write(str('size_source_x = ')     +str(size_source)           +str('\n'))
            file.write(str('size_source_y = ')     +str(size_source)           +str('\n'))
            file.write(str('size_source_z = ')     +str(size_source)           +str('\n'))
        file.write(str('dimension_telescope_x = ') +str(dimension_telescope_x) +str('\n'))
        file.write(str('dimension_telescope_y = ') +str(dimension_telescope_y) +str('\n'))
        file.write(str('END Parameters\n'))
        file.write(str('x y z vx vy vz phi theta\n')+str('\n'))

    for j in range(number_sources):
        for i in range(number_photons):
            if   generation_type == 'perfect-point':
                photon1 = photon(number_scatterings, center_sources[j],generation_frequency)
                print(' '+str(100*i/(number_photons*number_sources))+'%', end="\r")                #Shows the % of completion of the code

                if photon1.vxvyvz[2]>0 :
                    point     = photon1.starting_point
                    direction = photon1.vxvyvz
                    x_dir     = direction[0]*((length_telescope_to_object-point[2])/direction[2])+point[0]
                    y_dir     = direction[1]*((length_telescope_to_object-point[2])/direction[2])+point[1]
                    if( np.abs(x_dir)<dimension_telescope_x and np.abs(y_dir)<dimension_telescope_y) :
                        ax.scatter(x_dir,y_dir,s=0.3,color=wavelength_to_rgb(photon1.frequency),alpha=0.1)
                        number_capted_photons += 1
                    else :
                        pass



            elif generation_type == 'perfect-spherical':
                angle1  = rd.uniform(0,2*np.pi)
                angle2  = rd.uniform(0,  np.pi)
                if size_source_x != size_source_y and size_source_y != size_source_z:
                    print('The size of the source is not spherical thus the generation is false.')
                photon1 = photon(number_scatterings, [center_sources[j][0]+size_source_x*np.cos(angle1)*np.sin(angle2),center_sources[j][1]+size_source_x*np.sin(angle1)*np.sin(angle2),center_sources[j][2]+size_source_x*np.cos(angle2)],generation_frequency)
                print(' '+str(100*i/(number_photons*number_sources))+'%', end="\r")                #Shows the % of completion of the code
                if photon1.vxvyvz[2]>0 :
                    point     = photon1.starting_point
                    direction = photon1.vxvyvz
                    x_dir     = direction[0]*((length_telescope_to_object-point[2])/direction[2])+point[0]
                    y_dir     = direction[1]*((length_telescope_to_object-point[2])/direction[2])+point[1]
                    if( np.abs(x_dir)<dimension_telescope_x and np.abs(y_dir)<dimension_telescope_y) :
                        ax.scatter(x_dir,y_dir,s=0.3,color=wavelength_to_rgb(photon1.frequency),alpha=0.1)
                        number_capted_photons += 1
                    else :
                        pass

            elif generation_type == 'perfect-cubic':
                photon1 = photon(number_scatterings, [center_sources[j][0]+rd.uniform(-size_source_x,size_source_x),center_sources[j][1]+rd.uniform(-size_source_y,size_source_y),center_sources[j][2]+rd.uniform(-size_source_z,size_source_z)],generation_frequency)
                print(' '+str(100*i/(number_photons*number_sources))+'%', end="\r")                #Shows the % of completion of the code
                if photon1.vxvyvz[2]>0 :
                    point     = photon1.starting_point
                    direction = photon1.vxvyvz
                    x_dir     = direction[0]*((length_telescope_to_object-point[2])/direction[2])+point[0]
                    y_dir     = direction[1]*((length_telescope_to_object-point[2])/direction[2])+point[1]
                    if( np.abs(x_dir)<dimension_telescope_x and np.abs(y_dir)<dimension_telescope_y) :
                        ax.scatter(x_dir,y_dir,s=0.3,color=wavelength_to_rgb(photon1.frequency),alpha=0.1)
                        number_capted_photons += 1
                    else :
                        pass
            elif generation_type == 'perfect-rod':
                photon1 = photon(number_scatterings, [center_sources[j][0]+rd.uniform(-size_source_x,size_source_x),center_sources[j][1]+rd.uniform(-size_source_y,size_source_y),center_sources[j][2]+rd.uniform(-size_source_z,size_source_z)],generation_frequency)
                print(' '+str(100*i/(number_photons*number_sources))+'%', end="\r")                #Shows the % of completion of the code
                if photon1.vxvyvz[2]>0 :
                    point     = photon1.starting_point
                    direction = photon1.vxvyvz
                    x_dir     = direction[0]*((length_telescope_to_object-point[2])/direction[2])+point[0]
                    y_dir     = direction[1]*((length_telescope_to_object-point[2])/direction[2])+point[1]
                    if( np.abs(x_dir)<dimension_telescope_x and np.abs(y_dir)<dimension_telescope_y) :
                        ax.scatter(x_dir,y_dir,s=0.3,color=wavelength_to_rgb(photon1.frequency),alpha=0.1)
                        number_capted_photons += 1
                    else :
                        pass



    stop_time = timeit.default_timer()
    print(' ')
    print('time of simulation (in s): ', stop_time - start_time)
    print('There was '+str(number_photons)+' photons generated and '+str(number_capted_photons)+' were capted meaning '+str(100*number_capted_photons/number_photons)+'%')
    print(' ')
    plt.xlim(-1.1*dimension_telescope_x,1.1*dimension_telescope_x)
    plt.ylim(-1.1*dimension_telescope_y,1.1*dimension_telescope_y)
    #plt.show()
    plt.savefig('images/obstacle.png',dpi=300)
    #plt.show()
    plt.clf()

    '''
    #number of photons 100000 the telescope captation is 40*40
    #distance between planet and source is 10 and between planet and telescope is also 10.
    value_radius_planet    = [0          , 0.1       , 0.2       , 0.3       , 0.4       , 0.5       , 0.6       , 0.7       , 0.8       , 0.9       , 1         , 1.5      , 2        , 2.5      , 5       , 10 ]
    value_number_photons_1 = [20685/5000 , 4030/1000 , 3911/1000 , 3745/1000 , 3567/1000 , 3361/1000 , 3045/1000 , 2749/1000 , 2568/1000 , 2183/1000 , 1890/1000 , 876/1000 , 404/1000 , 186/1000 , 15/1000 , 0  ]
    error_bar = []
    for i in range(len(value_number_photons_1)):
        error_bar.append(value_number_photons_1[i]*5/100)
    plt.errorbar(value_radius_planet,value_number_photons_1,yerr=error_bar)
    plt.xlabel('radius of the planet')
    plt.ylabel('% of captation')
    plt.savefig('test.png',dpi=300)
    '''
