#calculate signal to noise of image twice once from images (as with real observations) once from theory and take the data cube produced by simulator and manipulate it to geet the 'initial' light curve back using input throughput and transmission spectra (and by using observations of known stars)
#Miranda J, Summer 2013

import sys, os
lib_path = os.path.abspath('/home/miranda/files/WIFIS_sim/run_files/')
sys.path.append(lib_path)

from pylab import *
import pyfits as py
import object as obj
import sky

#rebin images for better signal to noise?
rebinx=2
rebiny=2


#objet fils that want to annalyze (full path)
rundir='/home/miranda/files/WIFIS_sim/object_empty/'
#run number want to annalyze (string with /)
runnum='0/'

#make class object to hold all the infos
basisfilename = rundir+runnum+'basis.param'
object=obj.importObject(basisfilename)

#load 'observed' images and sky subtracted one
#sci=getdata('/home/miranda/files/WIFIS_sim/object/0/sci+sky_psf_0.fits')

sub=swapaxes(py.getdata(rundir+runnum+'noise_psf_0.fits'),0,-1)

scube=shape(sub)

#r=radius of apetture to use in arcsec
r=0.75
ry=ceil(r/object.YSCALE)
rx=ceil(r/(object.SLICE/2.))


#get object data
sci=swapaxes(py.getdata(rundir+runnum+'psf_0.fits'),0,-1)

scin=zeros((floor(scube[0]/rebinx),floor(scube[1]/rebiny),floor(scube[2])))

#cen=centre of the spatial array
cen=(scube[0]/(2),scube[1]/(2))

goodx=array([])
goody=array([])
goody.dtype='int'
goodx.dtype='int'
edgex=array([],dtype='int')
edgey=array([],dtype='int')
skyn=zeros((floor(scube[0]/rebinx),floor(scube[1]/rebiny),floor(scube[2])))
subn=zeros((floor(scube[0]/rebinx),floor(scube[1]/rebiny),floor(scube[2])))
for x in arange((floor(scube[0]/rebinx))):
	for y in arange(floor(scube[1]/rebiny)):

		subn[x,y,:]=apply_over_axes(sum,sub[x*rebinx:x*rebinx+rebinx,y*rebiny:y*rebiny+rebiny,:],(0,1))


for x in arange(scube[0]):
	for y in arange(scube[1]):
		#if (cen[0]-x)**2/rx**2+(cen[1]-y)**2/ry**2<=1: #ellipse
		if ((cen[0]-x)**2<=rx**2) & ((cen[1]-y)**2<=ry**2): #rectangle
			goodx=append(goodx,int(x))
			goody=append(goody,int(y))
		#if (cen[0]-x)**2/rx**2+(cen[1]-y)**2/ry**2==1: #ellipse edge:
		if (((cen[0]-x)**2==rx**2) & ((cen[1]-y)**2<=ry**2)) | (((cen[0]-x)**2<=rx**2) & ((cen[1]-y)**2==ry**2)) : #rectangle edge
			edgex=append(edgex,int(x))
			edgey=append(edgey,int(y))


sncube=shape(sky)
#Integration time per frame
INT_SKY = object.INT_SKY
INT_SCI = object.INT_SCI
#Number of sky frames
NFRAMES_SKY = object.NFRAMES_SKY
NFRAMES_SCI = object.NFRAMES_SCI

#get sky and throughput
emis=loadtxt(rundir+'supplements/skySpectrum.txt')
trans=loadtxt(rundir+'supplements/sky_trans.txt')
thop=loadtxt(rundir+'supplements/throughput.txt')

good=(goodx,goody)


#load other noise
therm=object.THERM*object.QE
dark = object.DARK_CURRENT
read= object.READ_NOISE	



wl=arange(.9,1.35,(1.35-.9)/(2048))

c= 2.998e14 #micromters/s
h=6.626e-34 #Js

	
figure()
N=sqrt(apply_over_axes(mean,(sub**2),(0,1))-apply_over_axes(mean,(sub),(0,1))**2)*sqrt(20)*10
N=reshape(N,2048) #photons/s
N/=object.COL_AREA  #ph/s/m^2
N*=(h*c)/(wl) #J/s/M^2
N*=1000            #erg/s/cm^2
N=N/trans/thop
semilogy(wl,N)

noise1=sqrt((((dark+therm+emis))*INT_SCI+(read**2))*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI)

noise=sqrt(2*noise1**2)

noise=noise*sqrt(20)
noise/=object.COL_AREA  #ph/s/m^2
noise*=(h*c)/(wl) #J/s/M^2
noise*=1000*10
noise/=(thop*trans)         #erg/s/sm^2
#semilogy(wl,noise)


figure()
plot(N-noise)
show()

