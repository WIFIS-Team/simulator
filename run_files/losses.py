#takes convolved and rebinned data cubes and applies losses from sky transmission, throughput, QE and takes units from pix/s/m^2/arcsec^2/nm to just phot/sec/pix.

from pylab import *
import pdb
import scipy.ndimage
from scipy.interpolate import interp1d

def getfluxvalues(cube,object,through=ones(2048)): #if using doing sky don't have any effect from sky throughput so simply multiply by one for objecg through is the sky throughput array


#####NOTE: current throughput values aren't actually quite right, include QE don't include telescope loses!

	#get optics throughput array
	wl,thop=loadtxt(object.THROUGHPUT_FILE,unpack=1)
	#throughput array isn't the right size, interpolate to make it the righs size. 
	x=arange(object.BANDMIN,object.BANDMAX,object.DBANDWIDTH)
	wl=wl*1000
	thop=interp1d(wl,thop)
	thop=thop(x)
	savetxt(object.HOME+object.PARAM+'supplements/throughput.txt',thop)

	#different sintax for emission or real cube
	if isinstance(cube,ndarray)==1:
		cube*=thop
		#cube*=object.QE
		cube*=object.COL_AREA

	else:
		cube.cube*=through[newaxis,:]
		cube.cube*=thop[newaxis,:]
		#cube.cube*=object.QE
		cube.cube*=object.COL_AREA 
		

	return cube
