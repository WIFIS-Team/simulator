from scipy import integrate
from scipy.special import erf
from scipy import ndimage
from numpy import *
import imp
import sys
import scipy
import datacube
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
def getBasisDictionary(basisfilename):
	basis = imp.load_source('basis', basisfilename)
	keys = basis.__dict__.keys()
	values = basis.__dict__.values()
	basis_dict = {}

	for i in range(len(keys)):
		basis_dict[keys[i]] = values[i]
		if '__' in keys[i]: #Two Underscores!
			basis_dict.pop(keys[i])
	return basis_dict

#Assumes only two columns in file: key and value
def getDictionary(filename):
	infile = open(filename, 'r')

	key = []
	value = []
	for line in infile:
		if not line.startswith("#") and line.strip():
			data = line.split()
			if len(data) != 2:
				print 'filename= '+filename+' does not have a strict key,value structure (ie. two columned). Goodbye.'
				sys.exit()	
			else:
				key.append(data[0])
				try:
					a = float(data[1])
					value.append(a)
				except ValueError:
					value.append(data[1])
	infile.close()
	d = dict(zip(key,value))
	return d






