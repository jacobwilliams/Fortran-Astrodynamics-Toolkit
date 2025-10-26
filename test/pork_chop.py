#!/usr/bin/env python
####################################################################################################
"""
Generate an Earth-Mars Pork Chop Plot.

AUTHOR: Jacob Williams, 11/11/2014

"""

####################################################################################################

import os
import datetime
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import dates
import math

####################################################################################################
if __name__ == "__main__":
	"""main program"""	

	#########################
	def read_file(file):

		#........data file.......
		#n_t0,    74
		#n_dt,    21
		#2456658.500    ,      100.0000000    ,      21.76904085    ,      32.92329343    
	     	#...
     	
		jd_j2000 = 2451545.0   #julian date of J2000 epoch
		j2000_epoch = datetime.datetime(2000, 1, 1, 12, 0, 0)
		
		n_t0 = None
		n_tf = None
		it0 = -1
		itf = -1
		t0_prev = -99999
		
		f = open( file, 'r')
		
		for line in f:
			s = line.split(',')
			
			if (s[0]=='n_t0'):
			
				n_t0 = int(s[1].strip())
				
			elif (s[0]=='n_dt'):
			
				n_tf = int(s[1].strip())
				
				#size the arrays:
				t0  = np.zeros(shape=(n_t0, n_tf))
				tf  = np.zeros(shape=(n_t0, n_tf))
				c30_lt_pi = np.zeros(shape=(n_t0, n_tf))
				c30_gt_pi = np.zeros(shape=(n_t0, n_tf))
				c3f_lt_pi = np.zeros(shape=(n_t0, n_tf))
				c3f_gt_pi = np.zeros(shape=(n_t0, n_tf))
				
			else:
									
				d = float(s[0].strip())-jd_j2000  #days since j2000
				dt = float(s[1].strip())
							
				_t0 = j2000_epoch + datetime.timedelta(d) 	#departure epoch
				_tf = j2000_epoch + datetime.timedelta(d) + datetime.timedelta(dt) 	#arrival epoch			
				_c30_lt_pi = float(s[2].strip())
				_c30_gt_pi = float(s[3].strip())
				_c3f_lt_pi = float(s[4].strip())
				_c3f_gt_pi = float(s[5].strip())
				
				if ( _t0 != t0_prev ):
					it0 = it0+1
					itf = 0
					t0_prev=_t0
				else:
					itf = itf + 1
				
				t0[it0,itf]  = matplotlib.dates.date2num(_t0)
				tf[it0,itf]  = matplotlib.dates.date2num(_tf)
				c30_lt_pi[it0,itf] = _c30_lt_pi
				c30_gt_pi[it0,itf] = _c30_gt_pi
				c3f_lt_pi[it0,itf] = _c3f_lt_pi
				c3f_gt_pi[it0,itf] = _c3f_gt_pi
			
		f.close()
	
		return t0,tf,c30_lt_pi,c30_gt_pi,c3f_lt_pi,c3f_gt_pi

	#########################
	# from: http://stackoverflow.com/questions/19876882/print-string-over-plotted-line-mimic-contour-plot-labels
	def label_line(line, label_text, near_i=None, near_x=None, near_y=None, rotation_offset=0, offset=(0,0)):
	    """call 
	        l, = plt.loglog(x, y)
	        label_line(l, "text", near_x=0.32)
	    """

	    def put_label(i):
	        """put label at given index"""
	        i = min(i, len(x)-2)
	        dx = sx[i+1] - sx[i]
	        dy = sy[i+1] - sy[i]
	        rotation = np.rad2deg(math.atan2(dy, dx)) + rotation_offset
	        pos = [(x[i] + x[i+1])/2. + offset[0], (y[i] + y[i+1])/2 + offset[1]]
	        plt.text(pos[0], pos[1], label_text, size=9, rotation=rotation, color = line.get_color(),
	        ha="center", va="center", bbox = dict(ec='1',fc='1'))
	
	    x = line.get_xdata()
	    y = line.get_ydata()
	    ax = line.get_axes()
	    if ax.get_xscale() == 'log':
	        sx = np.log10(x)    # screen space
	    else:
	        sx = x
	    if ax.get_yscale() == 'log':
	        sy = np.log10(y)
	    else:
	        sy = y
	
	    # find index
	    if near_i is not None:
	        i = near_i
	        if i < 0: # sanitize negative i
	            i = len(x) + i
	        put_label(i)
	    elif near_x is not None:
	        for i in range(len(x)-2):
	            if (x[i] < near_x and x[i+1] >= near_x) or (x[i+1] < near_x and x[i] >= near_x):
	                put_label(i)
	    elif near_y is not None:
	        for i in range(len(y)-2):
	            if (y[i] < near_y and y[i+1] >= near_y) or (y[i+1] < near_y and y[i] >= near_y):
	                put_label(i)
	    else:
	        raise ValueError("Need one of near_i, near_x, near_y")
									
	#########################
	def draw_flight_time_line(x,y,c):
		"""draw a transfer time line"""
		
		l, = plt.plot(x,y,c)
		label_line(l, str(int(y[0]-x[0]))+' days', near_i=np.size(x)-1, rotation_offset=-25, offset=[6,5])

	#########################
	def draw_all_flight_time_lines():
		"""draw all the transfer time lines"""
		
		draw_flight_time_line(t0[:,0],tf[:,0],'r')	
		draw_flight_time_line(t0[:,np.size(t0,1)/4],tf[:,np.size(t0,1)/4],'r')
		draw_flight_time_line(t0[:,np.size(t0,1)/2],tf[:,np.size(t0,1)/2],'r')
		draw_flight_time_line(t0[:,3*np.size(t0,1)/4],tf[:,3*np.size(t0,1)/4],'r')
		draw_flight_time_line(t0[:,np.size(t0,1)-1],tf[:,np.size(t0,1)-1],'r')

	#########################
	def generate_pork_chop_data():
		"""
			Call the f2py fortran routine and get the pork chop plot data.
			See porkchop.f90 for details. 
		"""
		
		import porkchop
		
		#time ranges:
		initial_t0		= 0			# departure epoch
		delta_t0   		= 1			#
		final_t0		= 100		#	
		initial_dt		= 100		# time of flight 
		delta_dt		= 1			#
		final_dt		= 400		#
		
		#initial epoch:
		y = 2016
		m = 1
		d = 1
		
		#this is because I don't know if f2py will work for allocatable output arrays:
		n_t0 = len(range(initial_t0,final_t0+1,delta_t0))
		n_dt = len(range(initial_dt,final_dt+1,delta_dt))
		
		#get the data:
		t0,tf,c30_lt_pi,c30_gt_pi,c3f_lt_pi,c3f_gt_pi = porkchop.porkchop( n_t0,n_dt, 
													y,m,d,
													initial_t0,
													delta_t0,
													final_t0,
													initial_dt,
													delta_dt,
													final_dt)
	
		return t0,tf,c30_lt_pi,c30_gt_pi,c3f_lt_pi,c3f_gt_pi

####################################################################################################

	#read file generated by executable:
	#t0,tf,c30_lt_pi,c30_gt_pi,c3f_lt_pi,c3f_gt_pi = read_file('pork.csv')
	
	#call the wrapper to the fortran routine:
	t0,tf,c30_lt_pi,c30_gt_pi,c3f_lt_pi,c3f_gt_pi = generate_pork_chop_data()
	
	#plot options:
	contour_range = list(range(1,30,2))  #contour lines for c3 plots
	contour_range_sum = list(range(1,50,2))  #contour lines for c30+c3f plot
	#contour_range_sum.append(100)

	label_fontsize = 20
	title_fontsize = 20
	x_tick_step = 5.09
	y_tick_step = 20.0
	hfmt = dates.DateFormatter('%m/%d/%y')  #axes are dates
	fmt = '%i'  #format for contour labels

	### rcParams are the default parameters for matplotlib
	matplotlib.rcParams['font.size'] = 10.
	matplotlib.rcParams['font.family'] = 'Serif'
	matplotlib.rcParams['axes.labelsize'] = 10.
	matplotlib.rcParams['xtick.labelsize'] = 10.
	matplotlib.rcParams['ytick.labelsize'] = 10.
	
	#make a plot:
	fig = plt.figure(figsize=(30,10))
	
	ax = fig.add_subplot(131)
	draw_all_flight_time_lines()
	CS_1 = plt.contour(t0, tf, c30_lt_pi , contour_range, colors='k')
	plt.clabel(CS_1, inline=1, fontsize=10, fmt=fmt)
	CS_2 = plt.contour(t0, tf, c30_gt_pi , contour_range, colors='b')
	plt.clabel(CS_2, inline=1, fontsize=10, fmt=fmt)	
	ax.xaxis.set_major_formatter(hfmt)
	ax.yaxis.set_major_formatter(hfmt)
	plt.xlabel('Earth Departure Date',fontsize=label_fontsize)	
	plt.xticks(np.arange(t0.min(), t0.max()+1, x_tick_step))	#set axis ticks
	locs, labels = plt.xticks()	
	plt.setp(labels, rotation=90)
	plt.ylabel('Mars Arrival Date',fontsize=label_fontsize)
	plt.yticks(np.arange(tf.min(), tf.max()+1, y_tick_step))
	plt.title('Earth Departure C3 [$km^2 s^{-2}$]',fontsize=title_fontsize)
	plt.grid('on')
	lines = [ CS_1.collections[0], CS_2.collections[0] ]  #add a legend
	labels = ['< $\pi$ Transfer','> $\pi$ Transfer']
	plt.legend(lines, labels, loc=4)
	
	ax = fig.add_subplot(132)	
	draw_all_flight_time_lines()
	CS_1 = plt.contour(t0, tf, c3f_lt_pi , contour_range, colors='k', label='<pi')
	plt.clabel(CS_1, inline=1, fontsize=10, fmt=fmt)
	CS_2 = plt.contour(t0, tf, c3f_gt_pi , contour_range, colors='b', label='>pi')
	plt.clabel(CS_2, inline=1, fontsize=10, fmt=fmt)	
	ax.xaxis.set_major_formatter(hfmt)
	ax.yaxis.set_major_formatter(hfmt)
	plt.xlabel('Earth Departure Date',fontsize=label_fontsize)
	plt.xticks(np.arange(t0.min(), t0.max()+1, x_tick_step))	#set axis ticks
	locs, labels = plt.xticks()	
	plt.setp(labels, rotation=90)
	plt.yticks(np.arange(tf.min(), tf.max()+1, y_tick_step))
	plt.title('Mars Arrival C3 [$km^2 s^{-2}$]',fontsize=title_fontsize)
	plt.grid('on')
	lines = [ CS_1.collections[0], CS_2.collections[0] ]  #add a legend
	labels = ['< $\pi$ Transfer','> $\pi$ Transfer']
	plt.legend(lines, labels, loc=4)

	ax = fig.add_subplot(133)	
	draw_all_flight_time_lines()
	CS_1 = plt.contour(t0, tf, c30_lt_pi+c3f_lt_pi , contour_range_sum, colors='k', label='<pi')
	plt.clabel(CS_1, inline=1, fontsize=10, fmt=fmt)
	CS_2 = plt.contour(t0, tf, c30_gt_pi+c3f_gt_pi , contour_range_sum, colors='b', label='>pi')
	plt.clabel(CS_2, inline=1, fontsize=10, fmt=fmt)	
	ax.xaxis.set_major_formatter(hfmt)
	ax.yaxis.set_major_formatter(hfmt)
	plt.xlabel('Earth Departure Date',fontsize=label_fontsize)
	plt.xticks(np.arange(t0.min(), t0.max()+1, x_tick_step))	#set axis ticks
	locs, labels = plt.xticks()	
	plt.setp(labels, rotation=90)
	plt.yticks(np.arange(tf.min(), tf.max()+1, y_tick_step))
	plt.title('Earth Departure + Mars Arrival C3 [$km^2 s^{-2}$]',fontsize=title_fontsize)
	plt.grid('on')
	lines = [ CS_1.collections[0], CS_2.collections[0] ]  #add a legend
	labels = ['< $\pi$ Transfer','> $\pi$ Transfer']
	plt.legend(lines, labels, loc=4)
	
	#save it:
	plt.savefig('pork_chop.png')
	
	