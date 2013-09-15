#calculate signal to noise of image twice once from images (as with real observations) once from theory and take the data cube produced by simulator and manipulate it to geet the 'initial' light curve back using input throughput and transmission spectra (and by using observations of known stars)
#Miranda J, Summer 2013

import sys, os


##put all into a function, save all images as png files

def doAnalysis(path,directory,numruns):

	lib_path = os.path.abspath(path+'run_files/')
	sys.path.append(lib_path)

	from pylab import *
	import pyfits as py
	import object as obj
	import sky

#rebin images for better signal to noise?
#rebinx=2
#rebiny=2


#objet fils that want to annalyze (full path)
	rundir=path+directory+'/'
#loop over all run numbers want to annalyze (string with /) 
	for runnum in arange(numruns):
		runnum=str(runnum)+'/'

#make class object to hold all the infos
		basisfilename = rundir+runnum+'basis.param'
		object=obj.importObject(basisfilename)

#load 'observed' images and sky subtracted one
#sci=getdata('/home/miranda/files/WIFIS_sim/object/0/sci+sky_psf_0.fits')
		sky=swapaxes(py.getdata(rundir+runnum+'sky_psf_0.fits'),0,-1)
		sub=swapaxes(py.getdata(rundir+runnum+'noise_psf_0.fits'),0,-1)

		scube=shape(sub)

#r=radius of apetture to use in arcsec
		d=2
		dy=int(ceil(d/object.YSCALE))
		dx=int(ceil(d/(object.SLICE/2.)))
		rx=dx/2
		ry=dy/2

#get object data
		sci=swapaxes(py.getdata(rundir+runnum+'psf_0.fits'),0,-1)

		#scin=zeros((floor(scube[0]/rebinx),floor(scube[1]/rebiny),floor(scube[2])))

#cen=centre of the spatial array
		cen=(scube[0]/(2),scube[1]/(2))
		cen=array(cen,dtype='float')
		goodx=array([])
		goody=array([])
		goody.dtype='int'
		goodx.dtype='int'
		edgex=array([],dtype='int')
		edgey=array([],dtype='int')
		#skyn=zeros((floor(scube[0]/rebinx),floor(scube[1]/rebiny),floor(scube[2])))
		#subn=zeros((floor(scube[0]/rebinx),floor(scube[1]/rebiny),floor(scube[2])))
		#for x in arange((floor(scube[0]/rebinx))):
		#	for y in arange(floor(scube[1]/rebiny)):
		#		skyn[x,y,:]=apply_over_axes(sum,sky[x*rebinx:x*rebinx+rebinx,y*rebiny:y*rebiny+rebiny,:],(0,1))
		#		subn[x,y,:]=apply_over_axes(sum,sub[x*rebinx:x*rebinx+rebinx,y*rebiny:y*rebiny+rebiny,:],(0,1))
		#		scin[x,y,:]=apply_over_axes(sum,sci[x*rebinx:x*rebinx+rebinx,y*rebiny:y*rebiny+rebiny,:],(0,1))
		if dx%2==0: cen[0]+=0.5
		if dy%2==0: cen[1]+=0.5
		for x in arange(scube[0]):
			for y in arange(scube[1]):

		#if (cen[0]-x)**2/rx**2+(cen[1]-y)**2/ry**2<=1: #ellipse doesnt' work
				if ((cen[0]-x)**2<=rx**2) & ((cen[1]-y)**2<=ry**2): #rectangle
					goodx=append(goodx,int(x))
					goody=append(goody,int(y))
		#if ceil(abs(cen[0]-x))**2/rx**2+ceil(abs(cen[1]-y))**2/ry**2==1: #ellipse edge doesn't work
				if ((ceil(abs(cen[0]-x))**2==rx**2) & ((cen[1]-y)**2<=ry**2)) | (((cen[0]-x)**2<=rx**2) & (ceil(abs(cen[1]-y))**2==ry**2)) : #rectangle edge
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

		good1=(goodx,goody)
		shift=scube[1]-4-ry-cen[1]
		good2=(goodx,goody+int(shift))


#load other noise
		therm=object.THERM*object.QE
		dark = object.DARK_CURRENT
		read= object.READ_NOISE	



		tn=trans*thop
		scit=sci*(array(tn))[newaxis,:]
		scit*=object.COL_AREA 

		wl=arange(.9,1.35,(1.35-.9)/(2048))

#get final data cube
		hdu=py.open(rundir+runnum+'noise_psf_0.fits')
		img=swapaxes(hdu[0].data,0,-1)


#collecting area of telescope
		col=object.COL_AREA

		c= 2.998e17 #nanometers/x
		h=6.626e-34 #Js

		img/=(array(tn))[newaxis,:]
		img/=col 
	        hdu[0].data=swapaxes(img,0,-1)
		hdu.writeto(rundir+runnum+'spectrum_theory.fits')

		i=0
		thing=('Centre', 'Edge')
		comp=zeros((2,2048))
		for good in (good1,good2):
	
			SN_thry=[]
			SN=[]
			for z in arange(sncube[2]):
				subz=sub[:,:,z]
				skyz=sky[:,:,z]
		
				signal=sum(subz[good])
				noiseS=sqrt(mean(skyz**2)-(mean(skyz))**2)
				noise=sqrt(2*(noiseS*sqrt(len(goodx)))**2+(sqrt(signal*INT_SCI*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI))**2) 
				SN.append(signal/noise)


		#spect=sum(img[good],axis=0) #ph/s/m^2
				spe=signal
				spe/=col
				spe*=(h*c)/(wl[z]*1000) #J/s/m^2
				spe*=1000 #erg/s/cm^2
				cthing=(spe)/(signal/noise)*10
				comp[i,z]=cthing

				sciz=scit[:,:,z]

				signalt=sum(sciz[good])
				noise1=sqrt(((signalt+(dark+therm+emis[z])*len(goodx))*INT_SCI+(read**2*len(goodx)))*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI)
				noise2=sqrt((((dark+therm+emis[z])*len(goodx))*INT_SKY+read**2*len(goodx))*NFRAMES_SKY)/(INT_SKY*NFRAMES_SKY)
				noiset=sqrt(noise1**2+noise2**2)
				SN_thry.append(signalt/noiset)



			figure()

			subplot2grid((3,1),(0,0),rowspan=2)
			plot(wl,SN,label='Calculated from Simulation')
			plot(wl,SN_thry,label='From Theoretical Values')
			legend(loc=2)
			title('Signal to Noise over a '+str(d)+'" arcsec aperture at the '+thing[i])
			xlabel('Wavelength ($\mu$m)')
			ylabel('S/N')
			tick_params(axis='x',labelbottom='off')
			locator_params(axis='y',prune='lower')
			
			subplot2grid((3,1),(2,0))

			subplots_adjust(hspace=0)
			plot(wl,(array(SN)-array(SN_thry))/array(SN_thry)*100,'r-',label='Percent Difference')
			legend(loc=2)
	
			xlabel('Wavelength ($\mu$m)')
			savefig(rundir+runnum+'signal_noise_'+thing[i]+'.png')



########################################
####################light curve




	#find spectrum around source use cen position and apeture from signal to noise
			spec=mean(img[good],axis=0)

	#convert to flux density (ph/s/micron)
			spec/=(object.DBANDWIDTH/1000)



	#convert to energy from photons (J/s/micron/m^2)
			spec*=(h*c)/(wl*1000)

#convert to erg/s/A
#spec*=0.1

#load model spec in J/s/micron
#mspec=loadtxt(rundir+'supplements/object_spec_EnergyFluxDensity.txt')
#mx=arange(0.9,1.35,(1.35-.9)/len(mspec))
	
	#get theory spectrum from PSF image
			mspec=mean(sci[good],axis=0)
	#convert to flux density (ph/s/micron)
			mspec/=(object.DBANDWIDTH/1000)
	#convert to energy from photons (J/s/micron)
			mspec*=(h*c)/(wl*1000)
 
			figure()
			subplot2grid((3,1),(0,0),rowspan=2)
			plot(wl,spec,label='Reduced from Simulation')
			plot(wl,mspec,label='Original fed into Simulator')
			legend()
			tick_params(axis='x',labelbottom='off')
			locator_params(axis='y',prune='lower')
			title('Mean spectrum over a '+str(d)+'" arcsec appetrure on the '+thing[i])
			ylabel('Flux (10$^{-15}$J/s/m$^2$/$\mu$m)')
	
			subplot2grid((3,1),(2,0))
			plot(wl,(spec-mspec)/mspec*100,'r-',label='Percent Difference')
			legend()
			subplots_adjust(hspace=0)
			xlabel('Wavelength ($\mu$m)')

	
			savefig(rundir+runnum+'reduced_spect_'+thing[i]+'.png')

	
	#######################################
	####calculate flux for one signal to noise
	#spec:(J/s/micron/m^2)




	###come up with a plot showing dominant noise source
			dat=mean(scit[good],axis=0)
			scp=dat/(INT_SCI*NFRAMES_SCI)
			skp=emis/(INT_SCI*NFRAMES_SCI)+emis/(INT_SKY*NFRAMES_SKY)
			dp=dark/(INT_SCI*NFRAMES_SCI)+dark/(INT_SKY*NFRAMES_SKY)
			tp=therm/(INT_SCI*NFRAMES_SCI)+therm/(INT_SKY*NFRAMES_SKY)
			rd=read**2/(INT_SKY**2*NFRAMES_SKY)+read**2/(INT_SCI**2*NFRAMES_SCI)
	
			noise1=sqrt(((dat+(dark+therm+emis))*INT_SCI+(read**2))*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI)
			noise2=sqrt((((dark+therm+emis))*INT_SKY+read**2)*NFRAMES_SKY)/(INT_SKY*NFRAMES_SKY)
			noise=(noise1**2+noise2**2)

			figure()
			plot(scp/noise,label='source noise')
			plot(skp/noise,label='sky noise')
			plot(dp/noise,label='dark noise')
			plot(tp/noise,'m-',label='thermal noise')
			plot(rd/noise,'c-',label='read noise')

	#sum at each wavelelnght should = one:
	#print min(rd/noise+tp/noise+dp/noise+skp/noise+scp/noise),max(rd/noise+tp/noise+dp/noise+skp/noise+scp/noise)

			legend()
			title('Contribution of Various Noise Sources at the '+thing[i])
			ylabel('Percentage of Total Noise')
			xlabel('Pixel (wavelength)')
			savefig(rundir+runnum+'noise_'+thing[i]+'.png')
			i+=1
#####

#####

#plot how flux at 10sigma varies
#figure()

#semilogy(wl,comp[1,:])

#ylabel('Flux at 10 $\sigma$ (erg/s/cm$^2$)')
#xlabel('Wavelength ($\mu$m)')



#################################
###signal to noise on single pixel vs rebinned pixel

#sky and sub
#skyn and subn
#sn_rebin=zeros(2048)
#sn_multi=zeros((rebinx,rebiny,2048))
#for z in arange(2048):
#	noiseS=sqrt(mean(sky[:,:,z]**2)-(mean(sky[:,:,z]))**2)
#	signalr=subn[cen[0]/rebinx,cen[1]/rebiny,z]
#	noiseSr=sqrt(mean(skyn[:,:,z]**2)-(mean(skyn[:,:,z]))**2)
#	noiser=sqrt(2*(noiseSr)**2+(sqrt(signalr*INT_SCI*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI))**2) 
#	sn_rebin[z]=signalr/noiser
#	for x in arange(rebinx)+cen[0]:
#		for y in arange(rebiny)+cen[1]:
#		
#			signalm=sub[x,y,z]
#			noisem=sqrt(2*(noiseS)**2+(sqrt(signalm*INT_SCI*NFRAMES_SCI)/(INT_SCI*NFRAMES_SCI))**2) 
#			sn_multi[x-cen[0],y-cen[1],z]=signalm/noisem
#sn_multi=apply_over_axes(mean,sn_multi,(0,1))
#sn_multi=sn_multi.reshape((2048))



#figure()
#subplot2grid((3,1),(0,0),rowspan=2)
#plot(wl,sn_multi,label='Original')
#plot(wl,sn_rebin,label='Re-binned')
#legend(loc=2)
#tick_params(axis='x',labelbottom='off')
#locator_params(axis='y',prune='lower')

#ylabel('S/N')
	
#subplot2grid((3,1),(2,0))
#plot(wl,sn_rebin/sn_multi,'r-',label='Times Improved')
#legend(loc=2)
#subplots_adjust(hspace=0)
#ylabel('Re-binned/Original')
#xlabel('Wavelength ($\mu$m)')

######################
###show final image summed over all wavelelngths with apeture overlayed
#		figure(figsize=(8, 3))
#		sub=swapaxes(py.getdata(rundir+runnum+'noise_psf_0.fits'),0,-1)
#		sub=sum(sub,axis=2)

#		edge=(edgex,edgey)
#		sub[edge]=amax(sub)+5
#		edge2=(edgex,edgey+int(shift))
#		sub[edge2]=amax(sub)+5
#		imshow(sub,interpolation='none',cmap='gray')

#		savefig(rundir+runnum+'final_img.png')

#show()
