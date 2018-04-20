#!/usr/bin/python
#
import os, sys
import matplotlib
matplotlib.use('Agg')
import numpy as np 
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#import math 
#import os

def l2_hgrid(top=399, bottom=0):
	'''
    return a height grid for caliop level-2 aerosol profile data
    '''
	#
	hgrid = np.zeros(top-bottom+1)
	for i in np.arange(bottom,top+1):
	   if ( i < 145  ):
	      hgrid[i] = -0.5 + 0.06 * i
	   elif( i < 345 ):
	      hgrid[i] = 8.2 + 0.06 * (i-145)
	   elif( i <=399 ):
	      hgrid[i] = 20.2 + 0.18 * (i-345)

	return hgrid

def l2_hcenter(top=399, bottom=0):

	hgrid = l2_hgrid(top=top, bottom=bottom)
	hcenter = ( hgrid[1:] + hgrid[:-1] ) / 2.0

	return hcenter


def get_extinction(hdf_file_name, domain):
	
	'''
    :param hdf_file_name: name of HDF file containing calipso sds
    :type hdf_file_name: string
    :returns, f_ext, f_lon, f_lat: 
    WHERE
    	:f_ext: 2D array with the filtered extinction coefficient for the interest domain 
    	:f_lon: 2D array with longitude values for the interest domain 
    	:f_lat: 2D array with latitude values for the interest domain
    
    '''
	#-- Open hdf file and read interest sds
	sd	 = Dataset(hdf_f)
	ext  = sd.variables['Extinction_Coefficient_532']
	cld  = sd.variables['Cloud_Layer_Fraction']
	tlon = sd.variables['Longitude']
	tlat = sd.variables['Latitude']
	
	#-- filter data using the interest domain
	lon = tlon[:,0]
	lat = tlat[:,0]
	#ll_sub_l = np.where((lat >= minlat) & (lat <= maxlat) & (lon >= minlon) & (lon <= maxlon))
	ll_sub_l = np.where((lat >= domain[1]) & (lat <= domain[3]) & \
				(lon >= domain[0]) & (lon <= domain[2]))
	
	f_ext	= ext[ll_sub_l]
	f_cld 	= cld[ll_sub_l]
	f_lon 	= lon[ll_sub_l]
	f_lat 	= lat[ll_sub_l]
	
	# -- take care of values greater than 1 than 0
	f_ext[ f_ext < 0 ] = np.nan 
	f_ext[ f_cld > 0 ] = -999
	
	print( np.min(lon), np.min(lat), np.max(lon), np.max(lat) )
	
	return f_ext, f_lon, f_lat


def plot_calipso(lat, lon, ext_data, **kwargs):
	
	''' 
    :param lat: 2D array containing latitude data
    :type lat: array
    :param lon: 2D array containing longitude data
    :type lon: array
    :param ext_data: 2D array containing the data to plot Extinction Coefficient
    :type ext_data: array
    :param domain: list defining the desired map domain (llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlon)
    :type domain: list, optional
    :param cmap: name color map
    :type cmap: matplotlib color, optional
    :param vmin: minimum value for data in  map
    :type vmin: int, optional
    :param vmax: maximum value for data in map
    :type vmax: int, optional
    :param ymin: minimum value for altitude (y axes)
    :type ymin: int, optional
    :param ymax: maximum value for altitude
    :type ymax: int, optinal
    :param y_label: label for Y axis
    :type y_label: string, optional
    :param title: title of the plot
    :type title: string, optional
    :return: None
    
    '''
	
	domain 	= kwargs.get('domain', [np.min(lon), np.min(lat), np.max(lon), np.max(lat)])
	cmap 	= kwargs.get('cmap', plt.cm.jet)
	vmin 	= kwargs.get('vmin', 0.0)
	vmax 	= kwargs.get('vmax', 0.5)
	ymin 	= kwargs.get('ymin', 0)
	ymax 	= kwargs.get('ymax', 30)
	y_label = kwargs.get('y_label', "")
	title  	= kwargs.get('title', "")
	
	level 	= ext_data.shape[1]
	points 	= ext_data.shape[0]

	xdata = np.zeros((points,level))
	ydata = np.zeros((points,level))
	hcent = l2_hcenter()	

	print( domain )

	for idx in range(points):
		xdata[idx,:] = idx

	for idx, idx_r in enumerate(reversed(range(level))):
		#ydata[:,idx_r] = -0.5 + idx * 0.06
	        ydata[:,idx_r] = hcent[idx]
	
	#fig = plt.figure(figsize=(10,4))
	ax1 = plt.subplot()
	cmap.set_under("gray")
	imdata = np.ma.masked_invalid( ext_data )
	img1  = ax1.pcolormesh(xdata, ydata, imdata, vmin=vmin, vmax=vmax, cmap=cmap)

	# add x, y labels 
	ax1.set_title(title,fontsize=14)
	ax1.set_xlim(0, np.max(xdata)+.1)
	ax1.set_ylim(ymin, ymax)
	ax1.set_ylabel(y_label, fontsize=13)
	
	ax1.xaxis.set_ticks(np.arange(0, np.max(xdata)+5, np.max(xdata)/5.0 ))
	xaxis_idx = [int(idx) for idx in np.arange(0, np.max(xdata)+5, np.max(xdata)/5.0 )]
	xticklabels = [r"%.1f"%lat[i]+"\n" + r"%.1f"%lon[i] for i in xaxis_idx]
	ax1.xaxis.set_ticklabels( xticklabels )
	
	plt.colorbar(img1, extend='min')
	
	
	axins = fig.add_subplot(231)
	axins.set_position([0.1, 0.630, 0.25, 0.23])
	#axins = fig.add_axes([0.05, 0.630, 0.25, 0.23]) 
	map1 = Basemap(llcrnrlon=domain[0],llcrnrlat=domain[1],urcrnrlon=domain[2],
		 urcrnrlat=domain[3],resolution='l', ax=axins)
	map1.drawcoastlines(linewidth=0.5)
	map1.drawstates(linewidth=0.5, linestyle='solid', color='k', 
		antialiased=1, ax=None, zorder=None)
	map1.scatter(lon,lat,1,marker='o',color='b', zorder=2)
	
	ax1.annotate('Cloud', xy=(1, 0), xytext=(50, -10), va='bottom',  
			xycoords='axes fraction', textcoords='offset points')
	ax1.annotate('Lat',   xy=(1, 0), xytext=(30, -15), va='bottom',  
			xycoords='axes fraction', textcoords='offset points')
	ax1.annotate('Lon',   xy=(1, 0), xytext=(30, -28), va='bottom',  
			xycoords='axes fraction', textcoords='offset points')

	

hg = l2_hgrid()
hc = l2_hcenter()
#for i in range( np.size(hg) ):
#   print(i, hg[i], hc[i] )
#sys.exit()

#CAL_LID_L2_05kmAPro-Standard-V4-10.2017-08-23T15-26-33ZD.hdf
#CAL_LID_L2_05kmAPro-Standard-V4-10.2017-08-23T17-05-23ZD.hdf
#CAL_LID_L2_05kmAPro-Standard-V4-10.2017-08-23T18-44-18ZD.hdf
#CAL_LID_L2_05kmAPro-Standard-V4-10.2017-08-23T20-23-09ZD.hdf

dir_data = '/Dedicated/jwang-data/xxu69/data/CALIOP/L2_05kmAPro/'
hdf_file = 'CAL_LID_L2_05kmAPro-Standard-V4-10.2017-08-23T15-26-33ZD.hdf' 
hdf_f    = dir_data + hdf_file
str_time = hdf_file[35:-9]
img_file = 'img/' + str_time + '.png'
print( img_file )
#sys.exit()

#-- define interest domain
minlon = -60 
maxlon = -40 
minlat =  40
maxlat =  66
domain = [minlon, minlat, maxlon, maxlat]

#-- get sds from HDF file
ext, lon, lat = get_extinction(hdf_f, domain)

#-- create an instance of figure and defining its dimensions
fig = plt.figure(figsize=(10,4))

#-- define some plot attributes
cmap 	= plt.cm.rainbow 
y_label = 'Altitude km'
title 	= 'Extinction coefficient ('+str_time+ ')'

#-- call the function that creates calipso plots
plot_calipso(lat, lon, ext, cmap=cmap, y_label=y_label, title=title, domain=domain)

#-- save plot
#fig.suptitle('2006-10-12T07-59-00ZN', fontsize=14)
plt.savefig(img_file, bbox_inches='tight')


plt.cla()
plt.clf()
plt.close(fig)
