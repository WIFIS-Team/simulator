###make python figures of images from simulation fits files 

from pylab import *
import pyfits as py
from matplotlib.colors import LogNorm 


#objet fils that want to annalyze (full path)
rundir='/home/miranda/files/WIFIS_sim/object/'
#run number want to annalyze (string with /)
runnum='0/'

###show galfit model:

mod=swapaxes(py.getdata(rundir+runnum+'galfit_model.fits'),0,-1)
#mod=log10(mod)
im=imshow(mod,interpolation='none',cmap='gray',norm=LogNorm(vmin=amin(mod), vmax=amax(mod)+10))

colorbar(im)

savefig('galfit_model.png')


show()
