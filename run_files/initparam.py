from numpy import *
import sys
import shutil
import imp
import os
import stat
import glob  #A very useful pythonic version of grep. 
#f,i,fl, g, v, b = filter, instr, flux, gal, vel, basic
import module as mod
import pdb

NUMCPUS=2
psf = 'psf'


#Get basis_top.param file -- Admin file holding file names/tags
def constructParam(bparamname,directory,nruns,psffwhm,objectstuff,basicstuff,objectname,basicname,savestuff,savename):
	bparam = imp.load_source('bparam',bparamname)
	print bparam.PSPACE_FILE
	pspacefile = bparam.PSPACE_FILE #Get user's parameter space file paramspace.param

	#pspace = []
	#pspacename = []
	#pspacename, pspace = getPSpace(pspacefile)	

	#add '/' to directory
	directory=''.join((directory,'/'))	

#	pspace_dict = {}
#	for i in range(len(pspace)):
#		pspace_dict[pspace[i][0]] = pspace[i][1]
#		print pspace[i][0],  pspace[i][1]

	#Extract the PSF parameters, and delete them from pspace.
	#pspace_psf = getPSpaceElements(pspacename, pspace, psf)
	#delPSpaceElements(pspacename, pspace, psf)

	#num_fullcube, dum = generateParam(bparam, pspace, pspacename, createRun, pid = 0)

	#create three lists from all combinations of basic parameters
	basisstuff=[]
	for x in basicstuff[0]:
		for y in basicstuff[1]:
			for z in basicstuff[2]:
				for w in basicstuff[3]:
					basisstuff.append([x,y,z,w])

	for i in range(nruns):
		rundir = bparam.HOME+directory 
		#make directories
		
  		os.makedirs(rundir+ str(i) + '/')
		
		#make .param files
		npsf=mkpsffile(rundir,psffwhm,bparam,i)
		nobj=makeobjfile(rundir,objectstuff,objectname,bparam,i)
		#steps per run:
		numsteps=((11*nobj)*npsf+4)*nruns
		makebasisfile(rundir,basisstuff[i],basicname,bparam,directory,savestuff,savename,i)
		#func_createPSF = lambda w, x, y, z: createPSF(w, x, y, z, rundir)
		#num_psf, dum = generateParam(bparam, pspace_psf, 0, func_createPSF, pid = 0)
		writeExecutionScript(bparam.BASIS_PREFIX+'.param', i, rundir,numsteps,nruns)
	


	outname = bparam.HOME+directory +'run.script'
	outfile = open(outname, 'wb')
	outfile.write('for ((i = 0; i < ' + str(nruns) + '; i++))\n')
	outfile.write('do\n')
#	outfile.write('\techo "HI EVERYONE $i"\n')
	outfile.write('\tipython '+bparam.HOME+directory+'$i/run.py\n')
	outfile.write('done\n')
	outfile.write('echo "All Runs of Simulation Completed."\n')
	outfile.close()
	os.chmod(outname, stat.S_IRWXU)



#write psf.param files
def mkpsffile(rundir,psffwhm,bparam,nim):
	npsf=len(psffwhm)
	for i in arange(len(psffwhm)):
		pfile=open(rundir+ str(nim) + '/'+bparam.PSF_PREFIX+str(i)+bparam.TAIL,'w')
		pfile.write('PSF_FWHM'+'\t\t' + psffwhm[i] + '\n')
		pfile.close()
	return npsf

def makeobjfile(rundir,objectstuff,objectname,bparam,num):
	nobj=len(objectstuff[0])
	for i in arange(len(objectstuff[0])):
		ofile=open(rundir+ str(num) + '/'+bparam.OBJ_PREFIX+str(i)+bparam.TAIL,'w')
		for k in arange(len(objectstuff)):
			ofile.write(objectname[k]+'\t\t'+objectstuff[k][i]+'\n')
	return nobj

def makebasisfile(rundir,basicstuff,basicname,bparam,directory,savestuff,savename,num):
	outname=rundir+ str(num) + '/'+bparam.BASIS_PREFIX+bparam.TAIL
	copyFile(bparam.BASIS_TOP_TEMPLATEFILE, outname)
	basisoutfile = open(outname, 'a') #append
	#Additional lines:

	for i in arange(len(basicstuff)):
		basisoutfile.write(basicname[i]+'\t=\t'+basicstuff[i]+'\n')

	basisoutfile.write('SKY_TRAN_FILE = HOME + \'sky/mktrans_zm_\' + SKY_TRAN + \'.dat\'\n')
	basisoutfile.write('SKY_EM_FILE = HOME + \'sky/mk_skybg_zm_\' + SKY_EM + \'_ph.dat\'\n')

	basisoutfile.write('PARAM\t=\t\''+directory+'\'\n')
	basisoutfile.write('RUN_DIR\t=\t\''+rundir+ str(num) + '/'+'\'\n')
	for i in arange(len(savename)):
		basisoutfile.write(savename[i]+'\t=\t'+savestuff[i]+'\n')

	basisoutfile.write('\n\n\n')

	infile = open(bparam.BASIS_END_TEMPLATEFILE, 'r')
	for line in infile:
		basisoutfile.write(line)

	basisoutfile.write('\n\n')
	basisoutfile.close()




#writeExecutionScript:
#Description: The automated execution of GALINO is by shell script command. The command simply executes each generator run.py in each run folder. 
#Here is a template of of each run.py file. 
def writeExecutionScript(basisname, pid, rundir,numsteps,nruns):
	basisfilename = rundir+str(pid)+'/'+basisname
	basis_dict = mod.getBasisDictionary(basisfilename)

	
	outname = rundir+ str(pid) + '/'+'run.py'
	outfile = open(outname, 'w')
	outfile.write('import sys, os\n')
	outfile.write('lib_path = os.path.abspath(\'' + basis_dict['HOME'] + 'run_files/\')\n')
	outfile.write('sys.path.append(lib_path)\n\n')
	
	

	outfile.write('import re\n')
	outfile.write('import imp\n')
	outfile.write('from numpy import *\n')
	outfile.write('import glob\n')
	outfile.write('from time import time\n')
	outfile.write('import pp\n')
	outfile.write('import pyfits, rescale, losses,pdb\n')
	outfile.write('import module as mod\n')
	outfile.write('import sky\n')
	outfile.write('import noise\n')
	outfile.write('import pickle\n')
	outfile.write('import object as obj\n\n')

	outfile.write('try:\n')

	outfile.write('\ttimes = []\n')
	outfile.write('\tstart = time()\n')


	outfile.write('\tbasisfilename = \'' + basisfilename + '\'\n')
	outfile.write('\tobject=obj.importObject(basisfilename)\n')


	outfile.write('\tlog=open(object.HOME+object.PARAM+\'supplements/progress\',\'a\')\n')
	outfile.write('\t#what step out of how many\n')
	outfile.write('\tnumsteps='+str(numsteps)+'\n')
	outfile.write('\tstepnum=1+'+str((pid)*(numsteps/nruns))+'\n')
	outfile.write('\ttrackprog=[stepnum,numsteps]\n')

	outfile.write('\n\n')
	outfile.write('\ttimes.append(time()-start)\n')


	outfile.write('\toidl=[]\n')
	outfile.write('\tfor obj_name in glob.glob(\''+rundir+ str(pid) + '/'+basis_dict['OBJ_PREFIX']+'*.param\'):  \n') 
	outfile.write('\t\t oidl.append(re.findall(r\'[0-9]+\',obj_name)[-1])\n') 


	outfile.write('\trcubes=[]\n')

	outfile.write('\tfor psfname in glob.glob(\'' + rundir + str(pid) + '/'+ basis_dict['PSF_PREFIX'] +'*.param\'): \n')
	outfile.write('\t\tpid = re.findall(r\'[0-9]+\', psfname)\n')
	outfile.write('\t\tobcube=[]\n')
	
	outfile.write('\t\tfor oid in oidl:\n')
	outfile.write('\t\t\tthing= obj.importThing(oid,object)\n')
	outfile.write('\t\t\ttempcube,trackprog=(	obj.makeCube(object,\''+rundir+ str(pid) + '/'+'\',psfname,pid[-1],thing,oid,trackprog))\n')
	outfile.write('\t\tobcube.append(tempcube)\n')
	outfile.write('\t\tcube=[]\n')
	#outfile.write('\tlog.write(\'combining individual objects for psf \'+str(pid[-1])+\'. Step \'+str(trackprog[0])+\' of \'+str(trackprog[1])+\'\\n\')\n')
	#outfile.write('\ttrackprog[0]+=1\n')
	outfile.write('\t\tfor i in arange(len(obcube)): cube.append(obcube[i].cube)\n')
	outfile.write('\t\tcube=array(cube)\n')
	outfile.write('\t\tcube=sum(cube,axis=0)\n')
	outfile.write('\t\tobcube=obcube[0]\n')
	outfile.write('\t\tobcube.cube=cube\n')
	outfile.write('\t\tif object.SAVE_PSFCUBE==1: obcube.SaveFITS(\''+rundir+ str(pid) + '/'+'\'+object.OBJ_PREFIX+\'_\'+object.PSF_PREFIX+pid[-1]+\'.fits\')\n')
	outfile.write('\t\trcubes.append(obcube)\n')
	outfile.write('\tlog.write( \'================================================================================\\n \')\n')
	outfile.write('\tlog.write( \'have input image, apply additional observational effects:\\n\')\n')
	outfile.write('\tlog.write(\'\\tget sky transmission and absortipn spectra. Step \'+str(trackprog[0])+\' of \'+str(trackprog[1])+\'\\n\')\n')
	outfile.write('\ttrackprog[0]+=1\n')
	outfile.write('\tems,thr=sky.getSkySpectrum(object)\n\n')
	outfile.write('\tlog.write(\'\\tfind number of electrons read by detector from sky emission. Step \'+str(trackprog[0])+\' of \'+str(trackprog[1])+\'\\n\')\n')
	outfile.write('\ttrackprog[0]+=1\n')
	outfile.write('\tems=losses.getfluxvalues(ems,object)\n\n')	

	outfile.write('\tfor psfname in glob.glob(\'' + rundir+ str(pid) + '/'+ 'psf_*.param\'): \n')
	outfile.write('\t\tpid = re.findall(r\'[0-9]+\', psfname)\n')
	outfile.write('\t\tcube = rcubes[int(pid[-1])]\n')
	outfile.write('\t\tlog.write( \'\\tfind number of electrons read by detector from source (psf=\'+str(pid[-1])+\'). Step \'+str(trackprog[0])+\' of \'+str(trackprog[1])+\'\\n\')\n')
	outfile.write('\t\ttrackprog[0]+=1\n')
	outfile.write('\t\tcube=losses.getfluxvalues(cube,object,through=thr)\n')
	outfile.write('\t\tlog.write( \'\\tapply noise and find final observed image (psf=\'+str(pid[-1])+\'). Step \'+str(trackprog[0])+\' of \'+str(trackprog[1])+\'\\n\')\n')
	outfile.write('\t\ttrackprog[0]+=1\n')
	outfile.write('\t\tnoise.ApplyNoise(cube,ems, object,\''+rundir+ str(pid) + '/'+'\', pid[-1])\n')


	outfile.write('\ttimes.append(time()-times[-1])\n')
	outfile.write('\tif os.path.exists(object.RUN_DIR+\'saturation_map.fits\'): log.write(\'ERROR: some pixels saturated, see saturation_map.fits\\n\')\n')
	outfile.write('\tlog.write( \'Simulation completed.\\n\')\n')
	outfile.write('except Exception, error:\n')
	outfile.write('\tlog=open(\''+rundir+'supplements/progress\',\'a\')\n')
	outfile.write('\tlog.write(\'WARNING! An Error Occured(\'+str(error)+\'\\n\')\n') 

def copyFile(inname, outname):
	shutil.copyfile(inname,outname)






