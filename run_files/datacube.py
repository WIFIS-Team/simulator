###writen by stephen Ro for Galino
###makes image arrays into 'data cubes' to be used by pyfits
###some edits by Miranda Jarvis

from numpy import *
import module as mod
import pyfits
import sys
#strictly holds a cube



import copy_reg
import types

###datacube class useful for storing data and saving fits files in particular
class Cube:
	#x,y are spatial coordinates -- in arcsecond units. 
	#z is the OBSERVED wavelength
	def __init__(self, arcsx, arcsy, z, cube, hdu_dict = {}):
		if (len(arcsx) != len(cube)):
			print 'X-axis and cube size do not match'
		elif (len(arcsy) != len(cube[0])):
			print 'Y-axis and cube size do not match'
		elif (len(z) != len(cube[0][0])):
			print 'Z-axis and cube size do not match'
		else:	
			
			self.arcsx = arcsx
			self.arcsy = arcsy
			self.z = z
			self.cube = cube
			self.hdu_dict = {}
			self.hdu_dict = dict(self.hdu_dict.items() + hdu_dict.items())


	def addCube(self, addcube):
		self.cube = self.cube + addcube

		
#	def getPhysX(self, redshift):
#		return mod.binwidthkpc(self.arcsx, self.redshift) #in kpc
		
#	def getPhysY(self, redshift):
#		return mod.binwidthkpc(self.arcsy, self.redshift) #in kpc
		
	def addDict(self, new_dict):
		self.hdu_dict = dict(self.hdu_dict.items() + new_dict.items())

	def getDict(self):
		return self.hdu_dict
		
	##used to save fits files with headers defined by hdu_dict items
	def SaveFITS(self, filename, hdu_dict = {}):
		self.hdu_dict =dict(self.hdu_dict.items() + hdu_dict.items())
		dum_hdu = self.hdu_dict
		keys = dum_hdu.keys()
		values = dum_hdu.values()
		for i in range(len(keys)):
			if keys[i].find('_') > -1:
				dum_hdu.pop(keys[i])
				keys[i] = keys[i].translate(None, "_")
				dum_hdu[keys[i]] = values[i]
		hdu = self.GenerateFITS()
		self.WriteFITS(filename, hdu)

	def GenerateFITS(self):
		lenslet = abs(self.arcsx[1]-self.arcsx[0])
		dz = abs(self.z[1]-self.z[0])
		dat=self.cube
		dat=swapaxes(dat,0,2)
		hdu = pyfits.PrimaryHDU(dat)
		header = hdu.header
		#header.update('CRPIX1',  0, comment = 'reference pixel location')
		#header.update('CRPIX2',  0, comment = 'reference pixel location')
		#header.update('CRPIX3',  0, comment = 'reference pixel location')

		#header.update('CUNIT1','nm      ', comment = 'Vacuum wavelength unit is nanometers')
		#header.update('CUNIT2','deg     ', comment = 'Degrees')
		#header.update('CUNIT3','deg     ', comment = 'Degrees')

		#header.update('CDELT1', dz, comment = 'nm/channel')
		#header.update('CDELT2', lenslet/360./60./60., comment = 'deg/pixel')
		#header.update('CDELT3', lenslet/360./60./60., comment = 'deg/pixel')

		#header.update('CRVAL1', self.z[0], comment = 'reference pixel location')
		#header.update('CRVAL2', 0., comment = 'reference pixel location')
		#header.update('CRVAL3', 0., comment = 'reference pixel location')

		#header.update('LENSLET', lenslet			,comment = 'Pixel angular width [arcsec]')		

		try:
			for i in range(len(self.hdu_dict)):
				if len(self.hdu_dict.keys()[i]) < 9:
					header.update(self.hdu_dict.keys()[i], self.hdu_dict.values()[i])
				else:
					header.update((self.hdu_dict.keys()[i])[0:7], self.hdu_dict.values()[i])
		except:
			print 'GenerateFITS(): Missed variable values? Will continue to write FITS file, may miss some header items.'

		return hdu


	def WriteFITS(self, filename, hdu):
		hdulist = pyfits.HDUList([hdu])
		#hdulist = pyfits.CompImageHDU(hdu.data, hdu.header)
		hdulist.writeto(filename, clobber = True)
		print 'FITS file saved: ' + filename + '\n\n'


#	def UpdateHeader(self, varname, value, comment = ''):
		


