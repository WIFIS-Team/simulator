#deals with creating dictionary items and initial object data cube(s)

from pylab import *
import module as mod
import datacube
import subprocess
import pyfits as py
from scipy.interpolate import interp1d
from scipy.special import gammainc
from scipy.signal import fftconvolve
from spec_fix import *
import pdb

###creates the object class which holds dictionary items to be called later.
###this class is called ty the impotObject and importThing functions defined bellow
class Object:
	def __init__(self,dic):
		self.dic = dic
		
		for i in range(len(dic)):
			vars(self)[dic.keys()[i]] = dic.values()[i]






###function: importObject
###input: basisfilename - filename of a parameter file
###takes a parameter file and creats on object with each entry in that file callable as dictionary entries.
###return 'object' class containing information from parameter file
def importObject(basisfilename):
	#uses module.py to create dictionary
	basis_dict = mod.getBasisDictionary(basisfilename)
	#creates object using that dictionary
	object = Object(basis_dict)

	return object
###function: importThing (similar to importObject)
###input: 	oid - number identifying which of the simulated objects this 'thing' represents
###		object - the object class item created by importObject
###return 'object' class containing info for that 'thing'
def importThing(oid,object):
	#string of the filename of the parameter file of the object number oid
	obfile=object.RUN_DIR+object.OBJ_PREFIX+str(oid)+object.TAIL
	#use module.py to create dictionary
	obdict=mod.getDictionary(obfile)
	#make object class out of dictionary
	thing=Object(obdict)
	return thing

###function: makeCube
###input:	object - thing created by importObject
###		rundir - directory with simulation in it
###		psfname - name of psf file being used for this image (can run multiple psfs in paralel)
###		pid - psf number being used for this image
###		thing - specific object in simulation being modeled in this image (all are later summed)
###		oid - number of the specific object being modeled
###		trackprog - step number used for the progress log file
###output: object data cube and updated trackprog value

def makeCube(Object,rundir,psfname,pid,thing,oid,trackprog): #make galaxy model in galfit, how to fit spectum to it????????????????
	log=open(Object.HOME+Object.PARAM+'supplements/progress','a')
	log.write('================================================================================\n')
	log.write( 'making data cube for object '+str(oid)+' with psf '+str(pid)+':\n')
	##make galfit template file
	log.write( 'Use Galfit:\n')
	log.write( '\tmake template. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	FOV,ellipse,s=makeTemplate(Object,thing,rundir)
	log.write( '\trun galfit. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	##run galfit

	with open(Object.HOME+Object.PARAM+'supplements/galfit_output', "a") as outfile:
    		subprocess.call([Object.GALFIT,rundir+'galfit_template'], stdout=outfile)
	
	log.write('\tload galfit image. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	##load galfit model file
	model=py.getdata(rundir+'galfit_model.fits') #units are electrons/s (totalled over all wavelength range??)
	model=swapaxes(model,0,1)

	NX=Object.NX/2	#make one x pixel per slice will be spread into two later
	NY=Object.NY
	NZ=Object.NZ

	model=model/sum(model)
	arcx=Object.SLICE*arange(NX)
	arcy=Object.YSCALE*arange(NY)
	
	wx=arange(Object.BANDMIN,Object.BANDMAX,Object.DBANDWIDTH)

	log.write( '================================================================================\n')	
	log.write('Load sample spectrum and manipulate until useful:\n')
	log.write('\t Load spectrum. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	#load sample spectrum:
	wl,fl=loadtxt(Object.OBJECT_FILE+thing.SPEC,unpack=1)
	
	log.write( '\tscale spectrum. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	#find scale factor to go from spec to right magnitude
	m=thing.MJ

	Fv=Object.LVEGA
	Jmin=Object.JMIN
	Jmax=Object.JMAX

	#flux based of magnitude
	Fs=Fv/(10**(m/2.5))

	#flux from imported spec over J band
	jrange=where((wl>=Jmin)&(wl<=Jmax))[0] #region in J band
	Fspec=mean(fl[jrange]) #flux

	#find scale factor
	c=Fs/Fspec

	#multiply whole spectrum by scale factor
	
	fl*=c

	###convolve in spectral direction

	R=Object.RES

	sigma =(1100./R)/(2.3548200)

	sigma=sqrt(sigma**2-((wl[1]-wl[0])/(2.3548200*2))**2)
	log.write( '\tconvolving in spectral axis. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1

	spec=convolveSpec(fl,len(wl),wl[1]-wl[0],sigma)

	savetxt(Object.HOME+Object.PARAM+'supplements/object_spec_EnergyFluxDensity.txt',matrix([wl,fl]).transpose())
	#find part of spectrum in wifis range
	wrange=where((wl>=(Object.BANDMIN-Object.DBANDWIDTH/2.))&(wl<=Object.BANDMAX))[0]
	spec=fl[wrange]
	wav=wl[wrange]
	del wl
	del fl
	

	#convert to photons from energy
	c= 2.998e17 #nanometers/x
	h=6.626e-34 #Js
	spec*=wav/(h*c)
	
	#change pixel scale of spec to match wifis . . . 


	
	#integrate over each WIFIS pixel
	intr=interp1d(wav,spec)
	new_spec=intr(wx)

	nnew_emis=append(new_spec,spec)
	tx=append(wx,wav)

	u=unique(tx,return_index=1)
	nnew_emis=nnew_emis[u[1]]
	tx=tx[u[1]]

	spect=tx+1j*nnew_emis
	spect=list(spect)
	wemis=[]
	wx=append(wx,1350)
	log.write('\tintegrate to get to right pixel scale. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1

	
	for i in arange(size(wx)-1)+1:
		intg=array([])
		k=spect[0].real
		while k<(wx[i]):
			a=spect.pop(0)
			k=spect[0].real
			intg=append(intg,(a.imag*0.001*(k-a.real)))
		wemis.append(sum(abs(intg)))
	spec=array(wemis)
	savetxt(Object.HOME+Object.PARAM+'supplements/object_spec_final.txt',spec)
	z=wx[:-1]

	log.write( '================================================================================\n')
	log.write('Combine into usable data cube:\n')
	#put flux into the normalized model (to save on compuatation time)
	log.write( '\tscale spatial model to right flux. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	allflux=sum(spec) #LRe=flux within ellipse
	#print allflux
	#normalize spectrum
	spec=spec/sum(spec)

	#find flux in ellipse for model
	fel=sum(model*ellipse)
	#imshow(model*ellipse)
	#show()
	ce=allflux/fel
	
	#save ellipse mask

	im=py.PrimaryHDU(ellipse) 
	im=py.HDUList([im])
	im.writeto(Object.HOME+Object.PARAM+'supplements/ellipse_mask.fits',clobber = True)
	
	#put in flux
	model*=ce
	#print 'flux in ellipse in model after multiplying in right amount',sum(model*ellipse)
	log.write( '\tconvolve spatial image. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
	#convolve model in spatial direction (seeing)
	model,psf_dict=constructPSF(model,psfname,Object)
	#print 'after convolving', sum(model*ellipse)

	log.write('\tcrop to right FOV. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
        ##get rid of extra pixels in model
        model=model*FOV
	model=model[where(model!=0)]

	model=reshape(model,(NX*Object.SLICE/.11,NY*Object.YSCALE/.11))
	#print 'size', shape(model)
	#print 'crop model', sum(model)
        #rebin galfit image

	
	dx=Object.SLICE/s
	dy=Object.YSCALE/s
	log.write( '\trebin image. Step '+str(trackprog[0])+' of '+str(trackprog[1])+' \n')
	trackprog[0]+=1
        gmodel=zeros((NX,NY))
        for x in arange(NX):
                for y in arange(NY):
			gmodel[x,y]=sum(model[x*dx:(x+1)*dx,y*dy:(y+1)*dy])
	model=gmodel
	#print 'rebin model', sum(model)

	#make right sized cube
	cube=ones((ceil(NX),ceil(NY),ceil(NZ)))

	cube*=spec
	
	cube=cube*model[:,:,newaxis]
	#print 'add spec', sum(cube)
	#also each slice in x is two thick
	cube=repeat(cube,2,axis=0)
	cube/=2.
	#print 'split x',sum(cube)
	#update arcsx,y,z
	NX*=2
	
	arcx=Object.SLICE/2*arange(NX)
	arcy=Object.YSCALE*arange(NY)

	hdu_dict = dict(thing.dic.items()+psf_dict.items())
       	cube=datacube.Cube(arcx,arcy,z,asarray(cube), hdu_dict)
	
	if Object.SAVE_OBJCUBE == True:

		try:	
			cube.SaveFITS(rundir+Object.OBJ_PREFIX+str(oid)+'_'+Object.PSF_PREFIX+pid+'.fits')
		except:
			print 'object: failed filename. Cannot save.'

	return cube,trackprog
	
###caled by makeCube to create galfit template file
###returns masks of the FOV on the galfit image, the ellips making the object out to defined radius and the pixel scale of the galfit image
def makeTemplate(Object,thing,rundir):

	#size of model pixels in arcsecs
	s=0.11

	v9=thing.v9 #axis ratio
	v10=(-1)*thing.v10
	t=v10*pi/180	

	gsx=int(thing.v1_1) #arcseconds left or right galaxy shifted in WIFIS FOV
	gsy=int(thing.v1_2) #arcseconds up or down galaxy shifted in WIFIS FOV

	#make each into an integer number of model pixels
	gsx=round(gsx/s) #galaxy shift in pixels
	gsy=round(gsy/s) #galaxy shift in pixels
	
	
	NX=Object.NX/2*Object.SLICE/s
	NY=Object.NY*Object.YSCALE/s	#size of WIFIS FOV in pixels (both axes)
	
	#if ((gsx>NX/2) | (gsy>NY/2)):
	#	print 'error galaxy shifted outside FOV'
	#	sys.exit()		

	modsz=Object.MODSZ
	
	a=thing.v4/s
	
	b=v9*a

	elx=[]
	ely=[]

	#check if a is big enough to include buffer
	
	ar=max(NY/2+abs(gsy)+NY*modsz,NX/2+abs(gsx)+NX*modsz)

	#make array the right size of semimajor axis
	area=zeros((ar*2,ar*2))

	#middle
	m=(ar)

	for x in arange(ar*2):
	    for y in arange(ar*2):
	        xc=x-m
	        yc=y-m
	        xp=((xc)*cos(t)-(yc)*sin(t))
	        yp=((xc)*sin(t)+(yc)*cos(t)) 
	        if (xp)**2/a**2+(yp)**2/b**2<=1:
	            elx.append(x)
	            ely.append(y)


	#make ellipse and FOV masks

	#make mask
	mFOV=zeros((ar*2,ar*2))
	mFOV[m-(gsx+NX/2):m+(NX/2-gsx),m-(gsy+NY/2):m+(NY/2-gsy)]=1
	ellipse=(elx,ely)
	mellipse=zeros((ar*2,ar*2))
	mellipse[ellipse]=1


	template=open(rundir+'galfit_template','w')
	template.write('B) '+rundir+'galfit_model.fits \n')
	template.write('H) 1    '+str(int(ar*2))+'   1    '+str(int(ar*2))+'\n')
	template.write('J) '+str(20)+'\n')
	template.write('K) '+str(s)+'  '+str(s)+'\n')
	template.write('O) regular \n')
	template.write('P) 1 \n')
	template.write('0) '+thing.TYPE+'\n')
	template.write('1) '+str(m)+' '+str(m)+'\n') #how far object is from edge of WIFIS FOV + how far edge of model is from edge of WIFIS
	template.write('3) '+str(4)+'\n')
	template.write('4) '+str(thing.v4/s)+'\n')
	template.write('5) '+str(thing.v5)+'\n')
	template.write('6) '+str(thing.v6)+'\n')
	template.write('7) '+str(thing.v7)+'\n')
	template.write('9) '+str(thing.v9)+'\n')
	template.write('10) '+str(thing.v10)+'\n')
#	if thing.SPIRAL == True:
#		template.write('R0) '+str(thing.r0)+'\n')
#		template.write('R1) '+str(thing.r1/s)+'\n')
#		template.write('R2) '+str(thing.r2/s)+'\n')
#		template.write('R3) '+str(thing.r3)+'\n')
#		template.write('R4) '+str(thing.r4)+'\n')
#		template.write('R9) '+str(thing.r9)+'\n')
#		template.write('R10) '+str(thing.r10)+'\n')
	
	template.close	
	
	return mFOV,mellipse,s

###function: constructPSF (called by makeCube)
###imputs:	img - image to apply psf over
###		psfname - file name of psf parameter file
###		object -object with useful dictionary terms
###returns image with psf applied and psf dictionary
def constructPSF(img,psfname,Object):
	#start = time()

	#get psf dictionary
	psf_dict = mod.getDictionary(psfname)
#	PSF_NSAMPLE = psf_dict['PSF_NSAMPLE']
#	PSF_XSIZE = psf_dict['PSF_XSIZE']
#	PSF_YSIZE = psf_dict['PSF_YSIZE']
	PSF_FWHM = psf_dict['PSF_FWHM']

	psfcube = applyCubePSF(img,PSF_FWHM,Object)


	#print 'PSF Convolution complete. dt='+str(int(abs(time()-start)))+'s  '
	
	return psfcube,psf_dict

###function: applyCubePSF (called by constructPSF)
###inputs:	img - to be convolved with psf
###		fwhm - of psf
###		object
###returns convolved image
def applyCubePSF(img, fwhm,Object):

	dx = Object.SLICE
	dy = Object.YSCALE


	sigma=fwhm/2.3548201
	#print 'spatial (Y,X)',dy/sigma,dx/sigma
	sigmax=sqrt(sigma**2-(dx/(2.3548200*2))**2)
	sigmay=sqrt(sigma**2-(dy/(2.3548200*2))**2)


        lx = int(ceil(2.0*4.0*sigmax/dx)) #(ie. 4 sigma, *2 for total length)
        if lx%2 == 0:
                lx+=1

        ly = int(ceil(2.0*4.0*sigmay/dy)) #(ie. 4 sigma, *2 for total length)
        if ly%2 == 0:
                ly+=1

	psf=array(outer(normal_weight(lx,0.0,sigmax,dx), normal_weight(ly,0.0,sigmay,dy)))
	
	img = fftconvolve(img, psf,mode='same')
	
	return img

	
