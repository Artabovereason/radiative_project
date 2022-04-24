import numpy as np
import matplotlib.pyplot as plt
import random as rd
import timeit
import scipy.stats
from scipy.stats import rv_continuous
from mpl_toolkits import mplot3d
import seaborn as sns
sns.set(style='ticks')

start_time = timeit.default_timer()
'''
    This code allow to check that the way we are generating the random frequency is indeed
    generated following Planck's law density.

    There is a normalization factor in the _pdf(self,x) : '(1/(5.670400*10**(-8)*temperature**4/np.pi))' (were 5.670400*10**(-8) is Stefan's constant)
    to ensure that the Planck's Law will be used as a normalized probability density.

    The class planck_gen is based on the rv_continuous subclass of the scipy.stats class.
    We can adapt the whole code at a different temperature using the temperature variable.
'''

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

def pdf(x):
    temperature = 6000
    h           = 6.62607004 * 10**(-34)
    c           = 299792458
    kb          = 1.38064852 * 10**(-23)
    if x.any() < 0:
        return 0
    else:
        return (2*h*c**2/(x**5))*(np.exp(h*c/(x*kb*temperature))-1)**(-1)

gaussian    = planck_gen(a=0.,b=3.5*10**(-6))
N           = 50000
save_random = gaussian.rvs(size=N)
to_plot     = np.linspace(0,3.5*10**(-6),1000)
fig,ax      = plt.subplots()
ax2         = ax.twinx()

ax.hist(save_random, 150, histtype='barstacked', density=True, alpha=0.4, edgecolor='none',label='Histogram')
ax2.plot(to_plot, pdf(to_plot),color='orange',label='Real Planck Distribution')
ax.set( xlabel='Wavelength [m]', ylabel='Number of randomly generated number in bin')
ax2.set(                  ylabel='Spectral radiance [W.sr$^{-1}$.m$^{-2}$.$\mu$m$^{-1}$]')
fig.suptitle('Representation of the Planck\'s law at the temperature $T=6000$ K')
plt.fill_between([0.38*10**(-6),0.75*10**(-6)],[10**(15),10**(15)],color='red',alpha=0.1,label='Visible light')


plt.xlim(0,3*10**(-6))
plt.ylim(0,4*10**(13))
ax.legend(loc='upper left')
ax.set_ylim(0,1.8*10**6)
plt.legend()
plt.savefig('images/Planck_generation.png',dpi=600)

stop_time  = timeit.default_timer()
print(' ')
print('time of simulation (in s): ', stop_time - start_time)
print(' ')
