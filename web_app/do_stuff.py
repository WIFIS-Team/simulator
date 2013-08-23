#try doing stuff as part of code that takes a long time
from pylab import *
thing=open('thing.txt','w')
s=[]
for i in arange(50000):
	s.append(sin(arange(i)))
	s=s[i:-1]
	thing.write(str(i/50000.*100)+'\n')
thing.close()
	

