#takes convolved and rebinned data cubes and applies losses from sky transmission, throughput, QE and takes units from pix/s/m^2/arcsec^2/nm to just phot/sec/pix.

from pylab import *
import pdb
import scipy.ndimage
from scipy.interpolate import interp1d

###function: getfluxvalues
###inputs: 	cube --> the sky or science images to apply scinece to sky is a 1d array of flux values, science is a data cube
###		object --> the dictionary of all the useful simulator info
###		through --> sky transmission array, defaults to ones since doesn't get applied to sky (this might not be necissary anymore?)
###outputs:	cube with appropriate losses applies same dimentions as input. 
def getfluxvalues(cube,object,through=ones(2048)): #if using doing sky don't have any effect from sky throughput so simply multiply by one for objecg through is the sky throughput array


#####NOTE: current throughput values aren't actually quite right, include QE don't include telescope loses!

	#get optics throughput array
	wl,thop=loadtxt(object.THROUGHPUT_FILE,unpack=1)
	#throughput array isn't the right size, interpolate to make it the righs size. 
	x=arange(object.BANDMIN,object.BANDMAX,object.DBANDWIDTH)
	wl=wl*1000
	thop=interp1d(wl,thop)
	thop=thop(x)
	#save array for future use
	savetxt(object.HOME+object.PARAM+'supplements/throughput.txt',thop)

	# need different sintax for sky emission or science cube so look at each seperatly
	#if cube is an array it's the sky
	if isinstance(cube,ndarray)==1:
		#multiply by instrument throughput
		cube*=thop
		
		#might need to apply QE seperatly, right now throughput includes QE
		#cube*=object.QE
		#multiply by collecting area
		cube*=object.COL_AREA
	#if not its the science cube
	else:
		#multiply along spectral axis by throughput of optics and sky
		cube.cube*=through[newaxis,:]
		cube.cube*=thop[newaxis,:]
		
		#might need to apply QE seperatly, right now throughput includes QE
		#cube.cube*=object.QE
		#multiply by collecting area
		cube.cube*=object.COL_AREA 
		

	return cube
