import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.interpolate import interp1d as interp

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

def invgamma(alpha, beta, seed):
    np.random.seed(seed)
    return 1 / np.random.gamma(alpha, 1 / beta) + 275

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
    seed = 5555
 
    
    # read the user-provided inputs from the command line (if there)
    if '-seed' in sys.argv:
        p = sys.argv.index('-seed')
        seed = int(sys.argv[p+1])
    if '-alpha' in sys.argv:
        p = sys.argv.index('-alpha')
        alpha = int(sys.argv[p+1])
    if '-beta' in sys.argv:
        p = sys.argv.index('-beta')
        beta = int(sys.argv[p+1])

    x = invgamma(alpha, beta, seed)
        
    param = [alpha, beta, seed, x]

    with open(r'PriorParameters.txt', 'w') as fp:
        for item in param:
            fp.write("%s\n" % item)
