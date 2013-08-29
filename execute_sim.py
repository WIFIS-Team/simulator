#run entire simulation in one go!!!
import sys 
import imp
from numpy import *
import os
import subprocess
import pdb


pathname='/Users/suresh/Source/wifis_simulator/' #string

path=os.path.abspath(pathname+'run_files/')
sys.path.append(path)

#read in system variables . . . 

#path=sys.argv[0] path to execute . . . 
savename=['SAVE_OBJCUBE','SAVE_PSFCUBE','SAVE_IND']
savestuff=[(sys.argv[1]),(sys.argv[2]),(sys.argv[3])]

directory=sys.argv[4]

basicstuff=[(sys.argv[5]).split(','),(sys.argv[6]).split(','),(sys.argv[7]).split(','),(sys.argv[8]).split(',')]

psffwhm=(sys.argv[9]).split(',')



objectstuff=[(sys.argv[10]).split(','),(sys.argv[11]).split(','),(sys.argv[12]).split(','),(sys.argv[13]).split(','),(sys.argv[14]).split(','),(sys.argv[15]).split(','),(sys.argv[16]).split(','),(sys.argv[17]).split(','),(sys.argv[18]).split(','),(sys.argv[19]).split(','),(sys.argv[20]).split(','),(sys.argv[21]).split(',')]


objectname=['MJ','SPEC','LABEL','TYPE','v1_1','v1_2','v4','v5','v6','v7','v9','v10']

basicname=['INT_SKY','INT_SCI','NFRAMES_SCI','NFRAMES_SKY']

#*************all values still strings***************

nruns=len(basicstuff[0])*len(basicstuff[1])*len(basicstuff[2])*len(basicstuff[3])

import initparam
initparam.constructParam(pathname+'basis_top.param',directory,nruns,psffwhm,objectstuff,basicstuff,objectname,basicname,savestuff,savename) ##will need to take in the directory name as well will also need to work out how many runs from user input, how many psf's too? maybe later don't want to write psf.param files so how to figure out how many times to do it with what values

###run run.script
log=open(pathname+directory+'/supplements/progress', "a") 
subprocess.call([pathname+directory+'/run.script'],shell=1,executable='/bin/bash')

log.write('All Simulations Complete, files now available for download')

subprocess.Popen(['find '+pathname+directory+'/ -type f \( -name "*.param" -o -name "*.py" -o -name "*.script" -o -name "*.paramc" \) -delete'],shell=1,executable='/bin/bash')


path=os.path.abspath(pathname+'/')
sys.path.append(path)

subprocess.call(['tar','-cvf',pathname+directory+'/'+directory+'.tar', '../'+directory+'/'])



