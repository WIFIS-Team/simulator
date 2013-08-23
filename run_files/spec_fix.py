#functions used to convolve in spectral direction (blur as well??)
from numpy import *
import module as mod
import scipy
from scipy import ndimage
from scipy.special import erf
from scipy import ndimage


def convolveSpec(spec, nx, dwidth, sigma,axisn=-1):
	#axis shows which axis to convolve on
	
        l = int(ceil(2.0*4.0*sigma/dwidth)) #(ie. 4 sigma, *2 for total length)
        if l%2 == 0:
                l+=1

        kernel = normal_weight(l, 0.0, sigma, dwidth)
        spec = scipy.ndimage.filters.convolve1d(spec, kernel, axis=axisn, mode = 'mirror')
        return spec

def normal_weight(l, mu, sigma, dz):
	
        #dz is the step size
        #l must be odd pos-integer.
        r = zeros(l)
        num_new =  int(l)/2
        denom = sqrt(2.0)*sigma
        numer = dz/2.0 - mu
        r[num_new] = 2.0*erf(numer/denom)

        for i in range(l/2):
                j = i+1
                numer1 = (i + 1.5)*dz - mu
                numer2 = (i + 0.5)*dz - mu
                a = erf(numer1/denom) - erf(numer2/denom)
                r[num_new+j] = a
                r[num_new-j] = a

        r = array(r)
        r = r/sum(r)

        return r	
