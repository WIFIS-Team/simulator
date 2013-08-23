from numpy import *
import scipy
import scipy.ndimage
from scipy.special import erf
from spec_fix import *
from scipy.interpolate import interp1d
import pdb
import gc


#Returns spectrum of sky emission and transmission
def getSkySpectrum(object):
	log=open(object.HOME+object.PARAM+'supplements/progress','a')
	wl1=object.BANDMIN
	wl2=object.BANDMAX

#import spectras
	#log.write( '\tget sky transmission and absortipn spectra. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	#trackprog[0]+=1
	emfile = object.SKY_EM_FILE
	tranfile= object.SKY_TRAN_FILE
	
	buffernm=(object.BANDMAX-object.BANDMIN)/16.

	wav,emflx=loadtxt(emfile,unpack=1)

	good=where(((wav>=wl1) | (abs(wav-wl1)<buffernm)) & ((wav<=wl2) | (abs(wav-wl2)<buffernm)))[0]
	wav=wav[good]
	emflx=emflx[good]
	tranflx=loadtxt(tranfile)[good,1]

	newstep=object.DBANDWIDTH
	#convolve to right line width
	R=object.RES
	
	sigmaw = (1100./R)/2.3548200 
	sigmas=0.045/2.3548200
	sigma=sqrt(sigmaw**2-sigmas**2)
	
	wav=append((object.BANDMIN-newstep/2),wav)
	emflx=append(mean(emflx[0:10]),emflx)
	emis=convolveSpec(emflx,len(wav),wav[2]-wav[1],sigma)
	tran=convolveSpec(tranflx,len(wav[1:]),wav[2]-wav[1],sigma)

	#print sky convolved
	del tranflx
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

	
	#integrate over each WIFIS pixel
	intr=interp1d(wav,emis)
	nnew_emis=intr(wx)

	nnew_emis=append(nnew_emis,emis)
	tx=append(wx,wav)

	u=unique(tx,return_index=1)
	nnew_emis=nnew_emis[u[1]]
	tx=tx[u[1]]

	spec=tx+1j*nnew_emis
	spec=list(spec)
	check=array(spec[:])
	wemis=[]
	wx=append(wx,1350)



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

	savetxt(object.HOME+object.PARAM+'supplements/sky_trans.txt',tran)
	return array(wemis),array(tran)







	
	












