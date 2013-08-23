###scale galaxy spectrum to match with observations (particularly the 1" slit used in Mould et al's paper)

from pylab import *
import pyfits as py
from scipy.special import gammainc
#make object file to get all the important stuff
import sys, os
lib_path = os.path.abspath('/home/miranda/files/WIFIS_sim/run_files/')
sys.path.append(lib_path)

import object as obj

#objet fils that want to annalyze (full path)
rundir='/home/miranda/files/WIFIS_sim/object/'
#run number want to annalyze (string with /)
runnum='0/'

#make class object to hold all the infos
basisfilename = rundir+runnum+'basis.param'
Object=obj.importObject(basisfilename)

wl,spec=loadtxt(rundir+'supplements/object_spec_EnergyFluxDensity.txt',unpack=1)
model=swapaxes(py.getdata(rundir+runnum+'galfit_model.fits'),0,-1)
model=model/sum(model) #normalize model

allflux=sum(spec)
	
#normalize spectrum
spec=spec/sum(spec)


	
ellipse=swapaxes(py.getdata(rundir+'supplements/ellipse_mask.fits'),0,-1)

fel=sum(model*ellipse)
	
ce=allflux/fel	
model*=ce

#take slice
sh=shape(model)
flux=sum(model[sh[0]/2-5:sh[0]/2+5,:])

spec*=flux*0.1-flux*0.01 #multiply back into spectrum turn into erg/cm/A

#wl in micro meters
wl*=0.001
wrange=where((wl>=(0.935))&(wl<=2.45))[0]
spec=spec[wrange]
wav=wl[wrange]

plot(wav,spec)
title('NGC 7562, spectrum used')
ylabel('Flux(ergs/cm$^2$/s/$\.A$)')
xlabel('$\lambda$ ($\mu$m)')
savefig('flux_for_slit_comp.eps')
show()






