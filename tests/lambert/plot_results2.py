#!/usr/bin/python

#plot results from lambert test2
#
#  JW : 10/31/2014

import matplotlib.pyplot as p

day2sec = 24.0*3600.0
sec2day = 1.0/day2sec
mu      = 398600.0  #earth k3/s2

#fig = p.figure(figsize=(20,5))
fig = p.figure(figsize=(20,10))

def read_traj(file):
	
	dt=[]
	alpha=[]
	#a=[]
	energy = []
	
	#  n, tof,  alpha, alpha, ...,
	#  2,     1.,   -2564910495.549713 ,  -1699999955.950112 ,
	
	#long n=0, long n=1,s=1, long n=1,s=2, short n=0, short n=1,s=1, short n=1,s=2
	
	ins = open( file, 'r')
	for line in ins:
		s = line.split(',')
		for i in range(2,len(s)):
			if (s[i].strip()!=''):
				dt.append(sec2day*float(s[1].strip()))   #days
				#dt.append( float(s[1].strip()) )   #sec
				alpha.append(  float(s[i].strip())        )
				energy.append( float(s[i].strip())/(-2.0) )  #energy
				#a.append(mu/float(s[i].strip()))
	ins.close()
	
	return dt,alpha,energy
	
dt,alpha,energy = read_traj('test2.csv')

dt_short,alpha_short,energy_short = read_traj('test2_short.csv')
dt_long,alpha_long,energy_long = read_traj('test2_long.csv')

#initial and final points:
r1x=20000
r1y=20000
r2x=-20000
r2y=10000

#p.plot(dt,alpha,'.',markersize=1)
#p.ylim(0, 18)
#p.ylabel('alpha ($\mathrm{km}^2 \mathrm{s}^{-2}$)')

p.plot(dt_short,energy_short,'b.',markersize=1,label='Short Way Transfer')
p.plot(dt_long,energy_long,'r.',markersize=1,label='Long Way Transfer')
p.ylim(-10,-2)
p.ylabel('Energy ($\mathrm{km}^2 \mathrm{s}^{-2}$)')

p.xlabel('$\Delta t$ (days)')
p.grid('on')
												
p.title("Lambert's Problem : $r_1=["+str(r1x)+","+str(r1y)+"] \/ \mathrm{km}$, $r_2=["+str(r2x)+","+str(r2y)+"] \/ \mathrm{km}$, $\mu=398600 \/ \mathrm{km}^3 \mathrm{s}^{-2}$")
p.legend(loc=4)

p.savefig('lambert2.png')

#test:
#p.ylim(-9,-8)
#p.xlim(0.5,1.5)


