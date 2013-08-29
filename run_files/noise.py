from numpy import *

import datacube
import sky
import module as mod
import sys
from time import time
import pyfits as py
import pdb


def ApplyNoise(cube, skySpectrum, object,  fileprefix, pid):
	start = time()

	savetxt(object.HOME+object.PARAM+'supplements/skySpectrum.txt',skySpectrum)
	nx = len(cube.arcsx) #cube x
	ny = len(cube.arcsy) #cube y
	nz = len(cube.z)     #cube z
	
	#The cube's spectral length should match the extracted sky spectrum 
	if nz != len(skySpectrum):
		print 'skySpectrum length does not match nz given'
		sys.exit()

	#Warning: SKY = sky, SCI = object
	#Integration time per frame
	INT_SKY = object.INT_SKY
	INT_SCI = object.INT_SCI
	#Number of sky frames
	NFRAMES_SKY = object.NFRAMES_SKY
	NFRAMES_SCI = object.NFRAMES_SCI

	THERM=object.THERM*object.QE
	DARK_CURRENT = object.DARK_CURRENT
	READ_NOISE = object.READ_NOISE	
	
		#dcube is the science image
	dcube = cube.cube  #[electrons/s]

	
	#make array of photon noise sources at each sky wavelength for the right integration time
	#componants are sky flux, thermal flux from instrument and dark current
	varSpectrum=(skySpectrum+THERM+DARK_CURRENT)*INT_SKY*NFRAMES_SKY

	skyCube = []

	for mu in varSpectrum: #Stepping through wavelength space. Generate nx*ny spatial plane.
		plane = random.poisson(mu, nx*ny)+random.normal(0,sqrt(READ_NOISE**2*NFRAMES_SKY),nx*ny) #readnnoise is gaussian, add seperatly
		plane = plane.reshape(nx,ny)
		skyCube.append(plane)

	#Rearrange axes to standard config (x,y,z)
	skyCube = array(skyCube)
	skyCube = swapaxes(skyCube, 0, 2)  #y, x, z
	skyCube = swapaxes(skyCube, 0, 1)  #x, y, z

	#add sky to dcube (use broadcasting so that sky spectrum gets added to each spatial pixel
	dcube+=skySpectrum[newaxis,:]
	#add other flux sourses
	varcube=(dcube+THERM+DARK_CURRENT)*INT_SCI*NFRAMES_SCI
	#apply 'noise'
	sciCube=random.poisson(varcube)+reshape(random.normal(0,sqrt(READ_NOISE**2*NFRAMES_SCI),nx*ny*nz),(nx,ny,nz))

	#check for saturation:
	#blank thing to hold saturation map
	satmap=zeros((object.NX/2,object.NY,object.NZ))
	#add up along x axis
	for i in arange(object.NX/2):
		row=sum(sciCube[i*2:(i+1)*2,:,:],axis=0)/NFRAMES_SCI
		#where is it saturated
		sat=where(row>=object.WELL_DEPTH)
		#set those points ==1
		row*=0
		row[sat]+=1
		satmap[i,:,:]=row
	#check if saturated anywhere at all:
	if amax(satmap>=1):
		#split satmap into same shape as data
		satmap=repeat(satmap,2,axis=0)
		im=py.PrimaryHDU(swapaxes(satmap,0,-1))
		im.header.add_comment('satruated pixels =1')
		im=py.HDUList([im])
		im.writeto(fileprefix+'saturation_map.fits',clobber=1)
		
	#Final Subtracted Image:

	sciCuber = sciCube/(INT_SCI*NFRAMES_SCI) - skyCube/(INT_SKY*NFRAMES_SKY) #Now electron/second


	###Convert to Flux again??? righ no e-/s /gain for ph./s
	

        #SAVE Cube(s)

        if object.SAVE_IND == True:
                datacube.Cube(cube.arcsx, cube.arcsy, cube.z, sciCube/(INT_SCI*NFRAMES_SCI)).SaveFITS(fileprefix+'sci+sky_psf_'+pid+'.fits')
		datacube.Cube(cube.arcsx, cube.arcsy, cube.z, skyCube/(INT_SKY*NFRAMES_SKY)).SaveFITS(fileprefix+'sky_psf_'+pid+'.fits')
	
	hdu_dict = dict( object.dic.items() )
	try:
		hdu_dict['PSFFWHM']=cube.hdu_dict['PSF_FWHM']
	except:
		hdu_dict['PSFFWHM']=cube.hdu_dict['PSFFWHM']
	sciCuber = datacube.Cube(cube.arcsx, cube.arcsy, cube.z, sciCuber, hdu_dict)
	sciCuber.SaveFITS(fileprefix+'noise_psf_'+pid+'.fits')
