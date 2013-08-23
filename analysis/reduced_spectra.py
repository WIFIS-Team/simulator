##take the data cube produced by simulator and manipulate it to geet the 'initial' light curve back
##do this using input throughput and transmission spectra (and by using observations of known stars)

from pylab import *
import pyfits as py

#objet fils that want to annalyze (full path)
rundir='/home/miranda/files/WIFIS_sim/object/'
#run number want to annalyze (string with /)
runnum='0/'

#get final data cube
hdu=py.open(rundir+runnum+'noise_psf_0.fits')
img=swapaxes(hdu[0].data,0,-1)

#get throughput and transmission
trans=loadtxt(rundir+'supplements/sky_trans.txt')
thop=loadtxt(rundir+'supplements/throughput.txt')

#collecting area of telescope
col=3.534


img/=(trans*thop)[newaxis,:]
img/=col 

hdu[0].data=swapaxes(img,0,-1)
hdu.writeto(rundir+runnum+'spectrum_theory.fits')



