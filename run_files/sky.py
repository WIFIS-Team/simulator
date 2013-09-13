##gets sky transmission and throughput files and puts them in a useful form for WIFIS

from numpy import *
import scipy
import scipy.ndimage
from scipy.special import erf
from spec_fix import *
from scipy.interpolate import interp1d
import pdb
import gc

##Function: getSkySpectrum
##input: object - holds basic parameter values
##Returns spectrum of sky emission and transmission as arrays
def getSkySpectrum(object):
	log=open(object.HOME+object.PARAM+'supplements/progress','a')
	#defind wifis spectral range
	wl1=object.BANDMIN
	wl2=object.BANDMAX

#import spectras
	#log.write( '\tget sky transmission and absortipn spectra. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	#trackprog[0]+=1
	
	#spectra file locations
	emfile = object.SKY_EM_FILE
	tranfile= object.SKY_TRAN_FILE
	#how much to keep to either side
	buffernm=(object.BANDMAX-object.BANDMIN)/16.
	
	#load emission file
	wav,emflx=loadtxt(emfile,unpack=1)
	#find areas within buffer region of WIFIS band width
	good=where(((wav>=wl1) | (abs(wav-wl1)<buffernm)) & ((wav<=wl2) | (abs(wav-wl2)<buffernm)))[0]
	#only look at those regions
	wav=wav[good]
	emflx=emflx[good]
	
	#load transmission spectrum elements within desired spectral range
	tranflx=loadtxt(tranfile)[good,1]
	
	#wavelength step of WIFIS
	newstep=object.DBANDWIDTH
	#convolve to right line width
	R=object.RES
	
	sigmaw = (1100./R)/2.3548200 
	sigmas=0.045/2.3548200
	sigma=sqrt(sigmaw**2-sigmas**2)
	
	#extrapolate data half a wavelength step from begining of spectral range to improve convolutions and rebinning.
	wav=append((object.BANDMIN-newstep/2),wav)
	emflx=append(mean(emflx[0:10]),emflx)
	#convolveSpec function in specfix.py
	emis=convolveSpec(emflx,len(wav),wav[2]-wav[1],sigma)
	tran=convolveSpec(tranflx,len(wav[1:]),wav[2]-wav[1],sigma)

	#print sky convolved
	
	del tranflx #to save memory space
	del emflx
	gc.collect()
	
	#only look at WIFIS range
	good=where((wav>=(wl1-newstep/2)) & (wav<=(wl2)))[0]
	wav=wav[good]
	emis=emis[good]
	tran=tran[good]

	#transmission just use interpolation:
	
	wx=arange(wl1,wl2,newstep)

	tran=interp1d(wav,tran)
	tran=tran(wx)

	
	#integrate over each WIFIS pixel in emision
	
	#get boundery points values through interpolation
	intr=interp1d(wav,emis)
	nnew_emis=intr(wx)

	#create master list with interpolated values and origional
	nnew_emis=append(nnew_emis,emis)
	tx=append(wx,wav)
	#sort and get rid of repeats
	u=unique(tx,return_index=1)
	nnew_emis=nnew_emis[u[1]]
	tx=tx[u[1]]
	
	#make wavelength and flux values into complex number for ease of manipulation
	spec=tx+1j*nnew_emis
	#make into a list
	spec=list(spec)
	check=array(spec[:])
	#blank list to be filled later
	wemis=[]
	#have end point
	wx=append(wx,1350)


	#integrate on boxes of length newstep centred on each WIFIS pixel, hense adding point one half step before first wifis pixel
	for i in arange(size(wx)-1):
		intg=array([])
		k=spec[0].real
		while k<(wx[i]+newstep/2):
			a=spec.pop(0)
			k=spec[0].real
			intg=append(intg,(a.imag*(k-a.real)))
		wemis.append(sum(abs(intg)))
	wemis=array(wemis)

	#multiply emission flux by pixel size
	wemis*=object.SLICE/2.*object.YSCALE	
	#save sky transmission file for later use
	savetxt(object.HOME+object.PARAM+'supplements/sky_trans.txt',tran)
	return array(wemis),array(tran)







	
	












