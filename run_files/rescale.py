#take image cube and convolve in spatial and spectral directions rebin to WIFIS pixel scale

from pylab import *
import datacube
import module as mod
from time import time
import scipy
from scipy import ndimage
from scipy import signal
import downsample as ds
from spec_fix import *
import pdb

def scaling(cube,object,psfname,filename,savenm):

	#convolve in spatial direction
	psfcube=constructPSF(cube,psfname)
	del cube
	#get data
	pcube=psfcube.cube



###rebin in both directions
	
	print 'changing pixel scale'
	NX=object.NX/2
	NY=object.NY
	NZ=object.NZ
	
	new_cube=pcube
	#also each slice in x is two thick
	new_cube=repeat(new_cube,2,axis=0)
	new_cube/=2.

	#update arcsx,y,z
	NX*=2
	
	arcx=object.SLICE/2*arange(NX)
	arcy=object.YSCALE*arange(NY)
	
	spec=arange(object.BANDMIN,object.BANDMAX,object.DBANDWIDTH)
	
	psfcube.arcsy=arcy
	psfcube.arcsx=arcx
	psfcube.z=spec
	psfcube.cube=new_cube

	
	if savenm == True & isinstance(filename, str):
		try:	
			psfcube.SaveFITS(filename)
		except:
			print 'PSF: failed filename. Cannot save.'
	

	return psfcube

def constructPSF(cube,psfname):
	start = time()
	arcsx = cube.arcsx
	arcsy = cube.arcsy
	z = cube.z

	psf_dict = mod.getDictionary(psfname)
#	PSF_NSAMPLE = psf_dict['PSF_NSAMPLE']
#	PSF_XSIZE = psf_dict['PSF_XSIZE']
#	PSF_YSIZE = psf_dict['PSF_YSIZE']
	PSF_FWHM = psf_dict['PSF_FWHM']

	psfccube = applyCubePSF(cube,PSF_FWHM)
	psfcube = datacube.Cube(arcsx, arcsy, z, psfccube, psf_dict)


	print 'PSF Convolution complete. dt='+str(int(abs(time()-start)))+'s  '
	
	return psfcube

def applyCubePSF(cube, fwhm):

	arcsx = cube.arcsx
	arcsy = cube.arcsy	
	z = cube.z
	dcube = cube.cube
	
	dx = abs(arcsx[1]-arcsx[0])
	dy = abs(arcsy[1]-arcsy[0])

	dcube = scipy.swapaxes(dcube, 0, 2)  #x <-> z: (z,y,x)
	dcube = scipy.swapaxes(dcube, 1, 2)  #x <-> y: (z,x,y)

#	print 'Convolving '+str(len(dcube[0]))+'x'+str(len(dcube[0][0]))+' plane with '+str(int(xsize/dx))+'x'+str(int(ysize/dx))+' PSF. DOWNSAMPLING = '+str(nsample)


	ccube = []
	sigma=fwhm/2.3548201
	#print 'spatial (Y,X)',dy/sigma,dx/sigma
	sigmax=sqrt(sigma**2-(dx/(2.3548200*2))**2)
	sigmay=sqrt(sigma**2-(dy/(2.3548200*2))**2)


        lx = int(ceil(2.0*4.0*sigmax/dx)) #(ie. 4 sigma, *2 for total length)
        if lx%2 == 0:
                lx+=1

        ly = int(ceil(2.0*4.0*sigmay/dy)) #(ie. 4 sigma, *2 for total length)
        if ly%2 == 0:
                ly+=1

	psf=array(outer(normal_weight(lx,0.0,sigmax,dx), normal_weight(ly,0.0,sigmay,dy)))

	
	for i in range(len(z)):
		plane = scipy.signal.fftconvolve(dcube[i], psf,mode='same')
		ccube.append(plane)


	ccube = asarray(ccube)
	ccube = scipy.swapaxes(ccube, 0, 2)  #z <-> x: (y,x,z)
	ccube = scipy.swapaxes(ccube, 0, 1)  #x <-> y: (x,y,z)

	return ccube




