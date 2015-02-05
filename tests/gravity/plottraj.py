#
# Read a trajectory file and plot it
#

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

def loadTrajectory(file):
	"""load the trajectory data from the file"""
	import re
	
	# format is et, x, y, z, vx, vy, vz
	
	re1='(\\s+)'                              # White Space
	re2='([+-]?\\d*\\.\\d+)(?![-+0-9\\.])'    # Float
	rg = re.compile(re1+re2+re1+re2+re1+re2+re1+re2+re1+re2+re1+re2+re1+re2,re.IGNORECASE|re.DOTALL)
	
	x = []
	y = []
	z = []
	
	f = open(file, 'r')
	for line in f:
		m = rg.search(line)
		if m:
			if (m.group(4)!=''):
				x = append(x, float(m.group(4)))
				y = append(y, float(m.group(6)))
				z = append(z, float(m.group(8)))
			
	return x,y,z		
	
def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

#load the trajectory [this is the output from the gravity test case]
x,y,z = loadTrajectory('../../bin/log2.txt')

#open the figure:
fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111, projection='3d')

#plot orbit:
ax.plot(x,y,z, label='trajectory')
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')

#plot sphere:
# Note: matplotlib doesn't do the line-hiding correctly for the orbit...
if (False):
	u = np.linspace(0, 2 * np.pi, 100)
	v = np.linspace(0, np.pi, 100)
	x = 6378.0 * np.outer(np.cos(u), np.sin(v))
	y = 6378.0 * np.outer(np.sin(u), np.sin(v))
	z = 6378.0 * np.outer(np.ones(np.size(u)), np.cos(v))
	ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='White', cmap=cm.coolwarm)
	axisEqual3D(ax)
	
plt.savefig('trajectory.png')
plt.show()

