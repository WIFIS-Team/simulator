#set up the web application for running the WIFIS simulator code, uses cherrypy and calls execute_sim.py to run simulation with user defined parameters
#Miranda J
#summer 2013

import cherrypy
import glob
import os.path
import os
import re
import glob
import sys
from numpy import *
import linecache
from cherrypy.lib.static import serve_file
import subprocess
import pdb

class InputExample:

	def index ( self ):
		##home page, lets user set up conditions for simulation
		# Read in HTML form code from file
		htmlstring = open('inputform.html', 'r').read()

		return(htmlstring)

	def RunCode(self,button,specf=None,save=None,directory=None,intsky=None,intsci=None,nframesci=None,nframesky=None,psffwhm=None,magnitude=None,spec=None,label=None,typef=None,v1_1=None,v1_2=None,v4=None,v5=None,v6=None,v7=None,v9=None,v10=None):
		#clicking submit buttosn on home page instigates this code, depends on which button pressed
		
		#if start simulation button
		if button=='start simulation':	
			##check that all parameters are given if not revert to defaults
			
			#if spectrum of previously uploaded file given, add 'temp' to path name
			if spec:
				spec=spec.split(',')
				for i in arange(len(spec)):
					spec[i]='temp/'+spec[i]
				spec=','.join(spec)
					
			if not intsky: intsky='360'
			if not intsci: intsci='360'
			if not nframesci: nframesci='10'
			if not nframesky: nframesky='10'
			if not psffwhm: psffwhm='2'
			if not magnitude: magnitude='9.7'
			if not spec: spec='elliptical_3.txt'
			if not label: label='____'
			if not typef: typef='sersic'
			if not v1_1: v1_1='0'
			if not v1_2: v1_2='0'
			if not v4: v4='22.95'
			if not v5: v5='4'
			if not v6: v6='0'
			if not v7: v7='0'
			if not v9: v9='0.96'
			if not v10: v10='15'
		
		
			#check if a directory was added
			if not directory:
				return 'error: no directory added'
			#if directory given already in use return error
			if os.path.exists('../'+directory):
				return 'Error: directory specified already exists pick another name'
			else: os.makedirs('../'+directory+'/supplements')

			#make sure have right number of parameters in each line of object stuff, if not return errors
			num=len(label.split(','))
			err=[]
			if len(magnitude.split(','))!=num: err.append('magnitude')
			if len(spec.split(','))!=num: err.append('spectrum')
			if len(v1_1.split(','))!=num: err.append('x position')
			if len(v1_2.split(','))!=num: err.append('y position')
			if len(v4.split(','))!=num: err.append('scale length')
			if len(v5.split(','))!=num: err.append('variable parameter 1')
			if len(v6.split(','))!=num: err.append('variable parameter 2')
			if len(v7.split(','))!=num: err.append('variable parameter 3')
			if len(v9.split(','))!=num: err.append('axis ratio')
			if len(v10.split(','))!=num: err.append('position angle')
			if len(typef.split(','))!=num: err.append('type')
			if len(err)!=0:
				err=', '.join(err)
				
				return 'not all object parameters have the same number of entries as there are labels (effected parameters: '+err+')'
			
			#start progress log
			log=open('../'+directory+'/supplements/progress','a')
			log.write('Setting up Simulation \n ')
	
			#figure out which files were requested and put into form usable by code
			initCubes=1
			sumCube=1
			skysci=1
			#if save:
			#	if size(array(save))>=2:
			#		save=array(save)
			#		for s in save:
			#			if s=='initCubes': initCubes=1 
			#			if s=='sumCube': sumCube=1
			#			if s=='skysci': skysci=1
					
			#	else: 
			#		if save=='initCubes': initCubes=1 
			#		if save=='sumCube': sumCube=1
			#		if save=='skysci':skysci=1
			#star the simulation by running execute_sim with all of the nescissary files piped to it
			subprocess.Popen(['ipython ../execute_sim.py '+str(initCubes)+' '+str(sumCube)+' '+str(skysci)+' '+directory+' '+intsky+' '+intsci+' '+nframesci+' '+nframesky+' '+psffwhm+' '+magnitude+' '+spec+' '+label+' '+typef+' '+v1_1+' '+v1_2+' '+v4+' '+v5+' '+v6+' '+v7+' '+v9+' '+v10], shell=True,stdin=None, stdout=None, stderr=None, close_fds=True) #opens in background doesn't wait to finish
			
			##present the check progress and download buttons, each links to the class of that name
			return """<html><body><h2>Simulation started, check progress at:</h2><a href="progress?directory=%s">PROGESS</a><br />"""%directory +\
				"""<h2>When Completed, view results at : </h2><a href="showResults?directory=%s">Results</a><br />"""%directory
#"""<h2>When Completed download files at:</h2><a href="download?directory=%s">DOWNLOAD</a><br />"""%directory
		##if the upload files button pressed
		if button=='upload files':
			#make sure can deal with simultanious uploads by making into list
			if type(specf)!=list:
				temp=[]
				temp.append(specf)
				specf=temp
			#look at each element in list
			for i in arange(len(specf)):
				
				#check if there is already a file with that name:
				if os.path.exists('../object_spec/'+specf[i].filename):
					return 'error, a spectrum already exists with that name, try again'
				#save spectra to object_spec/temp
				else:
					temp=''
					while True:
						data=specf[i].file.read()
						temp+=data
						if not data:
							break
					savefile=open('../object_spec/temp/'+specf[i].filename,'wb')
					savefile.write(temp)
					savefile.close()
			return 'spetra saved, use the names of your file(s) in the spectrum section of the main form.'
					

	index.exposed = RunCode.exposed = True



class progress:
	#page went to if progress button pressed on submit page
	def index(self, directory):
		#look at porgress file
		filename = '../'+directory+'/supplements/progress'
		#bits used to display last few lines of progress log, instead show all
		#p=subprocess.Popen(['tail','-n',str(5),filename], stdout=subprocess.PIPE)
		#s,sinput=p.communicate()
		#s=s.split('\n')
    	
		##read all of progress file and format each line such that it displays well 
		with open(filename, "r") as f:
			page = '<ul>%s</ul>' % "\n".join("%s<br/>" % line for line in f)
		##display contents of progress file
		return "<html><body><h2>Progress:</h2><br />"+page
		###stuff from displaying last few lines
		#"<html><body><h2>Progress:</h2><br />\n" + \
		#"%(s0)s<br />\n"%{'s0':s[0]} +\
		#"%(s1)s<br />\n"%{'s1':s[1]} +\
		#"%(s2)s<br />\n"%{'s2':s[2]} +\
		#"%(s3)s<br />\n"%{'s3':s[3]} +\
		#"%(s4)s<br />\n"%{'s4':s[4]}


	index.exposed = True
#for downloading ones favorite spectrum to use in the simulator
class download:
	def index(self,directory):
		#look for right file based on user input
		filename='../'+directory+'/'+directory+'.tar'
		#check if the tar file has been made yet, if so download it to user
		if os.path.exists(filename):
			filename=os.path.abspath(filename)
			return serve_file(filename, "application/x-download", "attachment")
		#if not return error	
		else: return 'Simulation Incomplete, try again later!'	
	index.exposed = True

class showResults:
	def index(self,directory):
		#need to look in each run folder
		files=os.listdir('../'+directory)
		#make into single string
		files=','.join(files)
		#find numbers since these are runs
		files=re.findall(r'[0-9]',files)
		#blank html to return
		html='<html><body>'
		#find png's in each run file make into list

		for x in files:
			x=str(x)
			html=html+'<h2> Results from run %s </h2>'%x

			imgs=glob.glob('../'+directory+'/'+x+'/*.png')
			for img in imgs:
				img=str(img)
				
				subprocess.call(['cp',img,'./images'])
				img=img.split('/')
				
				html=html+'<img src="images/%s" >'%img[-1]
			
			

		html=html+'<br/>'
		print html
		return html #"""<img src='images/final_img.png'>"""
	index.exposed= True


#for miranda
conf = {'/showResults/images': {'tools.staticdir.on': True,
        'tools.staticdir.dir': '/home/miranda/files/GitHub_files/wifis_simulator/web_app/images'}}

if __name__ == '__main__':
	root = InputExample()
	root.progress=progress()
	root.download=download()
	root.showResults=showResults()

	cherrypy.quickstart(root, '/', config=conf)



