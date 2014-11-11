#!/usr/bin/python

#plot results from lambert test
#
#  JW : 10/29/2014

import matplotlib.pyplot as p

fig = p.figure(figsize=(11,15))

def read_traj(file):
	
	x =[]
	y =[]
	
	ins = open( file, 'r')
	for line in ins:
		s = line.split(',')
		x.append(float(s[0].strip()))
		y.append(float(s[1].strip()))
	ins.close()
	
	return x,y

def plot_traj(file,linestyle='b'):
		
	x,y = read_traj(file)
	
	p.plot(x,y,linestyle, label=file.replace('.txt','').replace('_',' '), linewidth=2.0)
	
	p.legend(loc=4)
	p.axis('equal')
	p.grid('on')
	p.xlabel('x (km)')
	p.ylabel('y (km)')
	
plot_traj('shortway_n=0_s=1.txt','k')
plot_traj('longway_n=0_s=1.txt','k--')
plot_traj('shortway_n=1_s=1.txt','b')
plot_traj('longway_n=1_s=1.txt','b--')
plot_traj('shortway_n=1_s=2.txt','r')
plot_traj('longway_n=1_s=2.txt','r--')

#draw earth:											
circle1=p.Circle((0,0),6378,color='b')
fig = p.gcf()
fig.gca().add_artist(circle1)	

#draw initial and final points:

r1x=20000
r1y=20000
p.plot(r1x, r1y,'k.', markersize=20.0)
arr=p.arrow(0,0,r1x, r1y,  lw=3, length_includes_head=True, head_width=2000)
p.gca().add_patch(arr) 
p.text(r1x/2-5000, r1y/2+1000, '$r_1$', fontsize=20)

r2x=-20000
r2y=10000
p.plot(r2x, r2y,'k.', markersize=20.0)
arr=p.arrow(0,0,r2x,r2y,  lw=3, length_includes_head=True, head_width=2000)
p.gca().add_patch(arr) 
p.text(r2x/2, r2y/2+2000, '$r_2$', fontsize=20)											
												
p.title("Lambert's Problem : $r_1=["+str(r1x)+","+str(r1y)+"] \/ \mathrm{km}$, $r_2=["+str(r2x)+","+str(r2y)+"] \/ \mathrm{km}$, $\mu=398600 \/ \mathrm{km}^3 \mathrm{s}^{-2}$, $\Delta t =1 \/ \mathrm{day}$")

p.savefig('lambert.png')

