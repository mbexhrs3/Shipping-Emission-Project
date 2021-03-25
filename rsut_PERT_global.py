# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 16:09:38 2016

@author: acrnrrs
"""

from netCDF4 import*
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import matplotlib
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors

from scipy.stats import t
import scipy.stats
from scipy import stats
import csv 
import seaborn as sns
sns.set(color_codes=True)
sns.set(font_scale=1.5)  # controls the font size in the plots
 
import cdo 
from cdo import *
cdo = Cdo()
cdo.debug = True

matplotlib.rcParams.update({'font.size': 15,'font.family':'Times New Roman'})


################ Calculations using CDO ########################################################################################################
time_step1 =int(841)
time_step2 = int(1200)
#time_step1 =int(481)
#time_step2 = int(840)
TS1 = str(time_step1)
TS2 =  str(time_step2)

YEAR_moving_avg = np.linspace(6005,6095,91)
YEAR = np.linspace(6001,6100,100)


# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
#m_golbe = Basemap(projection='robin',lon_0=0,resolution='c')
m_golbe = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')


var_name='rsut'

infile = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc"

infile_seasonmean = infile
outfile_seasonmean = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_30years_"+TS1+"_"+TS2+"_seasonmean_global.nc"
# calculates seasonal mean (e.g. calculates 4 values) for given number of years (e.g. 30 years)
#cdo.yseasmean(input=infile_seasonmean, output=outfile_seasonmean)
cdo.yseasmean(input="-seltimestep,"+TS1+"/"+TS2+" /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_seasonmean)


infile_subdomain = infile
outfile_subdomain = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_yearmean_subdomain_global.nc" 
cdo.yearmean(input="/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_subdomain)
#to calculate time series for 100 years for the selected sub-domain

outfile_subdomain1 = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_30years_"+TS1+"_"+TS2+"_yearmean_subdomain_global.nc" 
cdo.yearmean(input="-seltimestep,"+TS1+"/"+TS2+" /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_subdomain1)
#calculates 30 years annual mean values for the given sub-domain and selected time-period TS1 and TS2

infile_subdomain_fldmean = infile
outfile_subdomain_fldmean = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_yearmean_subdomain_fldmean_global.nc" 
cdo.fldmean(input="-yearmean /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_subdomain_fldmean)
#to calculate grid-weighted average over the selected sub-domain 

############################
outfile_60_90N_fldmean = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_yearmean_subdomain_fldmean_60_90N.nc" 
cdo.fldmean(input="-yearmean -selindexbox,1,128,54,64 /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_60_90N_fldmean)
# to calculate mean weighted temperature for 60-90N 

############################# 
outfile_seasmean_DJF = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_30years_"+TS1+"_"+TS2+"_seasmean_DJF_global.nc" 
cdo.timselmean(3,11,9, input=" -seltimestep,"+TS1+"/"+TS2+" /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_seasmean_DJF)


outfile_seasmean_MAM = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_30years_"+TS1+"_"+TS2+"_seasmean_MAM_global.nc" 
cdo.timselmean(3,2,9, input=" -seltimestep,"+TS1+"/"+TS2+" /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_seasmean_MAM)

outfile_seasmean_JJA = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_30years_"+TS1+"_"+TS2+"_seasmean_JJA_global.nc" 
cdo.timselmean(3,5,9, input=" -seltimestep,"+TS1+"/"+TS2+" /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_seasmean_JJA)

outfile_seasmean_SON = "/HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/rsut_pert_30years_"+TS1+"_"+TS2+"_seasmean_SON_global.nc" 
cdo.timselmean(3,8,9, input=" -seltimestep,"+TS1+"/"+TS2+" /HOME/rrs/VF_2016/Simulations/CanESM4.1/model_output/pert_run/end/org_file/\
rsut_Amon_DevAM4-1_piControl-end_r1i1p1_600101-610012.nc", output=outfile_seasmean_SON)





################## CDO calculation ends here ####################################################################################################
#################################################################################################################################################



filename_SM = outfile_seasonmean 
filename_subdomain = outfile_subdomain
filename_fldmean = outfile_subdomain_fldmean
file_annual_mean_subdomain = outfile_subdomain1
file_60_90N = outfile_60_90N_fldmean

file_DJF = outfile_seasmean_DJF
file_MAM = outfile_seasmean_MAM
file_JJA = outfile_seasmean_JJA
file_SON = outfile_seasmean_SON



nc_file_SM = Dataset(filename_SM)
nc_file_subdomain = Dataset(filename_subdomain)
nc_file_fldmean = Dataset(filename_fldmean)
nc_file_annual_mean_subdomain = Dataset(file_annual_mean_subdomain)
nc_file_60_90N = Dataset(file_60_90N)

nc_file_DJF =  Dataset(file_DJF)
nc_file_MAM =  Dataset(file_MAM)
nc_file_JJA =  Dataset(file_JJA)
nc_file_SON =  Dataset(file_SON)


lats=nc_file_SM.variables['lat'][:]
lons=nc_file_SM.variables['lon'][:]
lon, lat=np.meshgrid([lons],[lats])

#######################################

############   TO PLOT A POINT AT THE CENTER OF GRID BOX ####################
############ USING MID-POINT OF THE GRID BOX ##############################

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

lon_mid = np.zeros(128)
lat_mid = np.zeros(16)


for i in my_range(1,15,1):
#    lat_mid[0] = 0
#    lat_mid[i] = (lats[i-1] + lats[i])/2    
#    print ((lats[i-1] + lats[i])/2.)
#    print (i-1)
    lat_mid[i] = ((lats[i-1] + lats[i])/2.)

lat_mid[[0]] = 46.04472653
    

for i in my_range(1,127,1):
    lon_mid[i] = ((lons[i-1] + lons[i])/2.)

lon_x, lat_y = np.meshgrid([lon_mid],[lat_mid])



############ ANNUAL MEAN OVER SUB-DOMAIN ###########
ANN_MEAN = nc_file_annual_mean_subdomain.variables[var_name][:]
ANN_MEAN1 = np.average(ANN_MEAN,axis=0)

############# SUBDOMAIN TIMESERIES ###########
subdomain_TS = nc_file_fldmean.variables[var_name][:]
TS_subdom = np.squeeze(subdomain_TS)
smoothed_TS = np.convolve(TS_subdom, np.ones(10)/10, mode='valid')
Global_avg = round(np.average(TS_subdom[-30:]),2)
print 'Last 30 Years global value is '
print var_name, 'PERT =', Global_avg

############  60N to 90N ###########################

sub_60_90N1 = nc_file_60_90N.variables[var_name][-30:]   #for last 30 years data
sub_60_90N = np.squeeze(sub_60_90N1)


############ SEASONAL MEAN OVER GLOBAL DOMAIN ###########
## Calculated to analyze statistical significance in the seasonal distribution #######

DJF_s = nc_file_DJF.variables[var_name][:]
MAM_s = nc_file_MAM.variables[var_name][:]
JJA_s = nc_file_JJA.variables[var_name][:]
SON_s = nc_file_SON.variables[var_name][:]


DJF = np.average(DJF_s,axis=0)
MAM = np.average(MAM_s,axis=0)
JJA = np.average(JJA_s,axis=0)
SON = np.average(SON_s,axis=0)



##########   PLOTTING BEGINS HERE ##########################################
###############################################################################
fig=plt.figure(2001,figsize=(13,8))
fig.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.1,hspace=0.1)

ax = fig.add_subplot(221)
x, y = m_golbe(lon, lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
cax = m_golbe.pcolor(x,y,DJF,shading='flat',cmap=my_cmap,vmin=10,vmax=300)
#cax = m_golbe.pcolor(x,y,T_DJF,shading='flat',cmap=my_cmap)
m_golbe.drawcoastlines(linewidth=0.5)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
#figure_title = '1 - Ratio GE 0.75  \
#                0 - Ratio LT 0.75'
figure_title = ' '

ax2  = plt.subplot(2,2,1)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)

plt.text(0.90, 1.15,'Outgoing SW flux (TOA)', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'a) DJF', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

############## ############ ##########
ax = fig.add_subplot(222)
my_cmap = matplotlib.cm.get_cmap('rainbow')
cax = m_golbe.pcolor(x,y,MAM,shading='flat',cmap=my_cmap,vmin=10,vmax=300)
#cax = m.pcolor(x,y,T_MAM,shading='flat',cmap=my_cmap)
m_golbe.drawcoastlines(linewidth=0.5)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
figure_title = ' '

ax2  = plt.subplot(2,2,2)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)

#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'b) MAM', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')



############## ############ ##########
ax = fig.add_subplot(223)
my_cmap = matplotlib.cm.get_cmap('rainbow')
cax = m_golbe.pcolor(x,y,JJA,shading='flat',cmap=my_cmap,vmin=10,vmax=300)
#cax = m.pcolor(x,y,T_MAM,shading='flat',cmap=my_cmap)
m_golbe.drawcoastlines(linewidth=0.5)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
figure_title = ' '

ax2  = plt.subplot(2,2,3)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)

#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'c) JJA', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')


############## ############ ##########
ax = fig.add_subplot(224)
my_cmap = matplotlib.cm.get_cmap('rainbow')
cax = m_golbe.pcolor(x,y,SON,shading='flat',cmap=my_cmap,vmin=10,vmax=300)
#cax = m.pcolor(x,y,T_MAM,shading='flat',cmap=my_cmap)
m_golbe.drawcoastlines(linewidth=0.5)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
figure_title = ' '

ax2  = plt.subplot(2,2,4)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)

#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'d) SON', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')


#plt.savefig('Figure_'+TS1+'_'+TS2+'_tas.png',dpi=150)
plt.savefig('Figure_801_rsut_g',dpi=150)

####################




fig=plt.figure(2002,figsize=(10,5))
plt.plot(YEAR, TS_subdom, linewidth=0.5, color='k',label='Outgoing SW flux (averaged over the sub-domain)')
plt.plot(YEAR_moving_avg, smoothed_TS, linewidth=2, color='darkgreen',label= '10-yrs moving average')




ax = plt.gca()
plt.ylabel('Outgoing SW flux at TOA [W/m$^{2}$]',fontsize=18, fontweight='normal')
plt.xlabel('Model year',fontsize=18, fontweight='normal')
plt.tick_params(which='major', width=2, length=7)
plt.tick_params(direction='out', top=False, right=False) # Turn ticks out
plt.legend(loc='upper left',frameon=False,fontsize='medium',ncol=1)

ax.set_xlim([6001,6100])
ax.set_xticks(np.linspace(6000,6100,6))
plt.tick_params(which='major', width=2, length=7)
plt.tick_params(direction='out', top=False, right=False) # Turn ticks out
plt.legend(loc='lower left',frameon=False,fontsize='medium',ncol=1)


plt.savefig('Figure_802_rsut_g.png',dpi=150)



#####################################################################################


fig=plt.figure(3003,figsize=(12,10))
fig.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.2,hspace=0.1)
ax = fig.add_subplot(221)
x, y = m_golbe(lon, lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
#cax = m_golbe.pcolor(x,y,mean_annual_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r,vmin=-1.,vmax=1.)
#cax = m_golbe.pcolor(x,y,mean_annual_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r)
cax = m_golbe.pcolor(x,y,ANN_MEAN1,shading='flat',cmap=my_cmap,vmin=10,vmax=300 )
#m_golbe.scatter(xx,yy,P_value_ANNUAL_sig_test_T,marker='D',color='k',linewidth=2.,vmin=1,vmax=1)
m_golbe.drawcoastlines(linewidth=1.0)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
#figure_title = '1 - Ratio GE 0.75  \
#                0 - Ratio LT 0.75'
figure_title = ' '

ax2  = plt.subplot(2,2,1)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)
#cbar.set_ticks([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
#cbar.ax.set_yticklabels([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])

#plt.text(0.90, 1.19,'Statistical Distribution', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.1,'a) Mean annual (Outgoing SW Radiation)', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

plt.savefig('Figure_803_rsut_g.png',dpi=150)
###########################################################################
