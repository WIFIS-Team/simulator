#find average flux in simulated images
from pylab import *

flux=array([927619.80932767701,927675.4692680802,927498.21996878227,927657.6611504145,927620.4692008479,927598.63132049784,927540.77233989444,927566.0677043664,927558.22211155551,927679.47349538829])

f=mean(flux)

stdv=sqrt(mean(flux**2)-mean(flux)**2)
error=stdv/sqrt(len(flux))

print f,error
t=927591.040473
print max(abs((flux-t)/t*100))
print (abs((f-t)/t*100))

