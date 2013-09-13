#functions used to convolve in spectral direction used by object.py and sky.py, closely based off functions writen by Stephen Ro
from numpy import *
import module as mod
import scipy
from scipy import ndimage
from scipy.special import erf
from scipy import ndimage

#function ConvolveSpec - convolve spectra to right resolution
#inputs:	spec - spectrum to be convolved, can be any size array, axis to be convolved on (spectral defined later)
#		nx - number of pixels
#		dwidth - lenght of pixel in nm
#		sigma - standard deviation in nm
#		asisn - axis to convolve over defualts to the last axis
#returns convolved spectrum

def convolveSpec(spec, nx, dwidth, sigma,axisn=-1):
	#axis shows which axis to convolve on
	
        l = int(ceil(2.0*4.0*sigma/dwidth)) #(ie. 4 sigma, *2 for total length)
        if l%2 == 0:
                l+=1

        kernel = normal_weight(l, 0.0, sigma, dwidth)
        spec = scipy.ndimage.filters.convolve1d(spec, kernel, axis=axisn, mode = 'mirror')
        return spec
#function: normal_weight , creates guassian psf kernel (uses error functions)
#inputs:	l - length of kernel found in convolveSpec
#		mu - mean
#		sigma - standard deviation
#		dz - size of pixel
#returns normalized guassian kernel
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
