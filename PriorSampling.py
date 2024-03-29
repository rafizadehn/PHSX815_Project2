import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.interpolate import interp1d as interp
from scipy.stats import invgamma
from scipy.special import gamma, gammainc

# import our Random class from Random.py file
sys.path.append(".")
from Random import Random

def boltz(v,m,T): # defines a function that returns a PROBABILITY for a specific set of parameters v, m, and T
    kB = 1.38e-23
    return (m/(2*np.pi*kB*T))**1.5 * 4*np.pi * v**2 * np.exp(-m*v**2/(2*kB*T))

def inv_boltz(v,m,T): # inverse function inverts it for ease of distribution calculation
    kB = 1.38e-23
    a = np.sqrt(kB*T/m)
    return erf(v/(np.sqrt(2)*a)) - np.sqrt(2/np.pi)* v* np.exp(-v**2/(2*a**2))/a

def velocities(n, seed):
    # this uses a random number generator from 
    # Dr. Christopher Rogan's GitHub repository.
    # If you cannot access this, numpy can be used with:
    # np.random.rand(n)
    random = Random(seed)
    vals = []
    for i in range(0, n):
        vals.append(random.rand())
    vel = inv_cdf(vals)
    return vel

def inv_gamma(x, alpha, beta, seed):
    np.random.seed(seed)
    return (beta ** alpha) / (gamma(alpha)) * (1/(x**alpha)) * np.exp(-beta/x)

if __name__ == "__main__":
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s -temp1 [temperature (k)] -temp2 [temperature (k)] -seed [seed] -Nmeas [number of particles sampled] -mass [mass in amu]" % sys.argv[0])
        print
        sys.exit(1)

    # default alpha
    alpha = 1

    # default beta
    beta = 2
    
    # default seed
    seed = 4836

    # default sample size
    Nsamp = 10000
 
    # read the user-provided inputs from the command line (if there)
    if '-seed' in sys.argv:
        p = sys.argv.index('-seed')
        seed = int(sys.argv[p+1])
    if '-alpha' in sys.argv:
        p = sys.argv.index('-alpha')
        alpha = float(sys.argv[p+1])
    if '-beta' in sys.argv:
        p = sys.argv.index('-beta')
        beta = float(sys.argv[p+1])

    # samples a temperature value from the inverse gamma 
    # distribution from a given seed value
    x_shift = 275
    x_values = np.linspace(x_shift+0.1, x_shift + 75, 1000)

    y_values = inv_gamma(x_values - x_shift, alpha, beta, seed)
    

    rng = np.random.default_rng(seed)

    # plot the inverse gamma distribution
    plt.plot(x_values, y_values)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Probability Density Function')
    plt.title('Inverse Gamma Function (alpha={}, beta={})'.format(float(alpha),float(beta)), fontweight = 'bold')
    plt.show()

    # puts the values used in the script into a single array 
    # so it can then be exported to a text file
    param = [alpha, beta, seed, rng.uniform(x_shift, x_shift+75)]

    # exports paramter array to text file
    with open(r'PriorParameters.txt', 'w') as fp:
        for item in param:
            fp.write("%s\n" % item)
    
