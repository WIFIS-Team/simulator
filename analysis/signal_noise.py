#calculate signal to noise of image twice once from images (as with real observations) once from theory
#Miranda J, Summer 2013

import sys, os
lib_path = os.path.abspath('/home/miranda/files/WIFIS_sim/run_files/')
sys.path.append(lib_path)

from pylab import *
import pyfits as py
import object as obj
import sky

#objet fils that want to annalyze (full path)
rundir='/home/miranda/files/WIFIS_sim/object/'
#run number want to annalyze (string with /)
runnum='0/'

#make class object to hold all the infos
basisfilename = rundir+runnum+'basis.param'
object=obj.importObject(basisfilename)

#load 'observed' images and sky subtracted one
#sci=getdata('/home/miranda/files/WIFIS_sim/object/0/sci+sky_psf_0.fits')
sky=swapaxes(py.getdata(rundir+runnum+'sky_psf_0.fits'),0,-1)
sub=swapaxes(py.getdata(rundir+runnum+'noise_psf_0.fits'),0,-1)

#Integration time per frame
INT_SKY = object.INT_SKY
INT_SCI = object.INT_SCI
#Number of sky frames
NFRAMES_SKY = object.NFRAMES_SKY
NFRAMES_SCI = object.NFRAMES_SCI


#r=radius of apetture to use
r=3

scube=shape(sub)

#cen=centre of the spatial array
cen=(scube[0]/2,scube[1]/2)

goodx=array([])
goody=array([])
goody.dtype='int'
goodx.dtype='int'
for x in arange(scube[0]):
	for y in arange(scube[1]):
		if sqrt((cen[0]-x)**2+(cen[1]-y)**2)<=r: 
			goodx=append(goodx,int(x))
			goody=append(goody,int(y))

good=(goodx,goody)

SN=[]
for z in arange(scube[2]):
	subz=sub[:,:,z]
	skyz=sky[:,:,z]
	
	signal=sum(subz[good])
	noiseS=sqrt(mean(skyz**2)-(mean(skyz))**2)
	noise=sqrt((noiseS*sqrt(len(goodx)))**2+(sqrt(signal*INT_SCI*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI))**2+(sqrt(len(goodx))*noiseS*sqrt((INT_SKY*NFRAMES_SKY)/(INT_SCI*NFRAMES_SCI)))**2) 
	SN.append(signal/noise)

x=arange(.9,1.35,(1.35-.9)/2048)
plot(x,SN,label='from data')

###calculate signal to noise by theory
#get object data
sci=swapaxes(py.getdata(rundir+runnum+'psf_0.fits'),0,-1)

#get sky and throughput
emis=loadtxt(rundir+'supplements/skySpectrum.txt')
trans=loadtxt(rundir+'supplements/sky_trans.txt')
thop=loadtxt(rundir+'supplements/throughput.txt')

sci*=(trans*thop)[newaxis,:]
sci*=object.COL_AREA 

#load other noise
therm=object.THERM*object.QE
dark = object.DARK_CURRENT
read= object.READ_NOISE	

SN_thry=[]
for z in arange(scube[2]):
	sciz=sci[:,:,z]

	signalt=sum(sciz[good])
	noise1=sqrt(((signalt+(dark+therm+emis[z])*len(goodx))*INT_SCI+(read**2*len(goodx)))*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI)
	noise2=sqrt((((dark+therm+emis[z])*len(goodx))*INT_SKY+read**2*len(goodx))*NFRAMES_SKY)/(INT_SKY*NFRAMES_SKY)
	noiset=sqrt(noise1**2+noise2**2)
	SN_thry.append(signalt/noiset)

plot(x,SN_thry,label='from theory')
legend(loc=2)
title('Signal to Noise over a '+str(r)+' radius appetrure centered at pixel ('+str(cen[0])+','+str(cen[1])+')')
xlabel('Wavelength ($\mu$m)')
ylabel('S/N')
#savefig('signal_noise.png')
show()


