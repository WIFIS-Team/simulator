#create the spectrum text file to be used in the simulator in order to model spiral galaxy's emission lines in arms
#Miranda Jarvis
#summer 2013

from pylab import *

#flux in units of 10e15 erg/cm^2/s

#flux of lines. line order is: 
#[[SIII] [SIII] Pa(delta) HeI Pa(gamma) [FeII] Pa(beta)]
line_wav=array([0.907, 0.953, 1.005, 1.083, 1.093, 1.257, 1.282 ]) #wavelength microns
line_wav*=1000 #nm

#flux from paper
###NGC 2903
line_flux=array([7.48, 9.52, 1.81, 7.95,  6.02, 3.53,  9.49]) #N2 aperture from NGC2903

#wavelength range of wifis
bandmin=0.9
bandmax=1.4

#steps size
stp_sz=0.0001

#flux without x103-15
line_flux*=10e-15
#to J/m^2
line_flux*=0.001 #J/m^2/s
#line_flux/=stp_sz #J/m^2/s/micron

#per arcsecond?
line_flux/=(0.8*2)

#wavelength array
wl=arange(bandmin*10000,bandmax*10000,stp_sz*10000)/10000
#wl in nm
wl*=1000 

#flux array (empty, put in lines later)
flux=zeros(len(wl))

#load elliptical spectrum for continuum
wav,ellps=loadtxt('elliptical_3.txt',unpack=1)



#find scale factor to go from spec to right magnitude
m=9.9

Fv=3.01e-9
Jmin=1160
Jmax=1360

#flux based of magnitude
Fs=Fv/(10**(m/2.5))

#flux from imported spec over J band
jrange=where((wav>=Jmin)&(wav<=Jmax))[0] #region in J band
Fspec=mean(ellps[jrange]) #flux

#find scale factor
c=Fs/Fspec

#multiply whole spectrum by scale factor
ellps*=c

ellps/=(pi*(18)**2)

#look at each wavelength of emission line and put that flux there in flux
for i in arange(len(line_wav)):
 
    #where is line on wl array?
    spot=where(wav<=line_wav[i])[0]
    spot=spot[-1]
    print ellps[spot],line_flux[i]
    ellps[spot]+=line_flux[i]


#convert to right units . . .


#flux without x103-15
#flux*=10e-15
#to J/m^2
#flux*=0.001 #J/m^2/s
#flux/=stp_sz #J/m^2/s/micron




#plot(wl,flux)

savetxt('spiral_arms+continuum.txt',matrix([wav,ellps]).transpose())

show()
