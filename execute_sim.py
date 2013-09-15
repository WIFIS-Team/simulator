#run entire simulation in one go!!!
#this code is called by app.py, primarilly it calls initparam.py which sets directories up and then runs each run.py file created
import sys 
import imp
from numpy import *
import os
import subprocess
import pdb

###primary running directory
#pathname='/Users/suresh/Source/wifis_simulator/' #string
pathname='/home/miranda/files/GitHub_files/wifis_simulator/' #for miranda's computer ;)

##make change current path to run_files directory (needed to run initparam.py
path=os.path.abspath(pathname+'run_files/')
sys.path.append(path)

#read in system variables provided from app.py, stored for now as strings in sys.argv

###sort into sets for ease of pushing to initparam.py, each set has one variable which has the name of each parameter and one which has the calue of each parameter as a string
###if multiple values given for a parameter it is given as a string with commas seperating different entries, split this string into list (in this case have lists as entries in list)

#denotes which files to save
savename=['SAVE_OBJCUBE','SAVE_PSFCUBE','SAVE_IND']
savestuff=[(sys.argv[1]),(sys.argv[2]),(sys.argv[3])]

#denotes basic properties of observation
basicname=['INT_SKY','INT_SCI','NFRAMES_SCI','NFRAMES_SKY']
basicstuff=[(sys.argv[5]).split(','),(sys.argv[6]).split(','),(sys.argv[7]).split(','),(sys.argv[8]).split(',')]

#denotes properties of objects to observe
objectstuff=[(sys.argv[10]).split(','),(sys.argv[11]).split(','),(sys.argv[12]).split(','),(sys.argv[13]).split(','),(sys.argv[14]).split(','),(sys.argv[15]).split(','),(sys.argv[16]).split(','),(sys.argv[17]).split(','),(sys.argv[18]).split(','),(sys.argv[19]).split(','),(sys.argv[20]).split(','),(sys.argv[21]).split(',')]
objectname=['MJ','SPEC','LABEL','TYPE','v1_1','v1_2','v4','v5','v6','v7','v9','v10']



directory=sys.argv[4] #this one isn't in a grouping
psffwhm=(sys.argv[9]).split(',') #neither is this

#*************all values still strings***************

#calculate how many runs will result from parameters given
nruns=len(basicstuff[0])*len(basicstuff[1])*len(basicstuff[2])*len(basicstuff[3])

#import initparam.py
import initparam
#run initparam, piping the necissary variables
initparam.constructParam(pathname+'basis_top.param',directory,nruns,psffwhm,objectstuff,basicstuff,objectname,basicname,savestuff,savename) ##will need to take in the directory name as well will also need to work out how many runs from user input, how many psf's too? maybe later don't want to write psf.param files so how to figure out how many times to do it with what values

##open the progress log
log=open(pathname+directory+'/supplements/progress', "a")
###run run.script, does all the runs set up by initparam)
subprocess.call([pathname+directory+'/run.script'],shell=1,executable='/bin/bash')


#delete files that we don't want to pass to the user (python codes, parameter files, run script . . .)
#subprocess.Popen(['find '+pathname+directory+'/ -type f \( -name "*.param" -o -name "*.py" -o -name "*.script" -o -name "*.paramc" \) -delete'],shell=1,executable='/bin/bash')

#change active path out of run_files
path=os.path.abspath(pathname+'/')
sys.path.append(path)

#tar the directory that the simulations were in so the user can download it
subprocess.call(['tar','-cvf',pathname+directory+'/'+directory+'.tar', '../'+directory+'/'])

#do analysis and save images
#go to right dir
path=os.path.abspath(pathname+'analysis/')
sys.path.append(path)
#import code
from analysis import *

#do it
doAnalysis(pathname,directory,nruns)

#write done message to progress log
log.write('All Simulations Complete, files now available for download')



