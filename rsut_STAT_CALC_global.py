# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:57:01 2017

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

 
import cdo 
from cdo import *
cdo = Cdo()
cdo.debug = True



import colorbar_cmap


import rsut_CTRL_global
import rsut_PERT_global


matplotlib.rcParams.update({'font.size': 15,'font.family':'Times New Roman'})


###############################################################################
ALPHA = 0.05


def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

###############################################################################

# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
#m_golbe = Basemap(projection='robin',lon_0=0,resolution='c')
## For some reason 'robin' map projection did not work very well, did not nicely 
## show meridional structure and also showed discontinuity from east of Greenwich 

#m_golbe = Basemap(projection='mill',lon_0=180)

m_golbe = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')

##############################################################################
lats= rsut_CTRL_global.lats
lons= rsut_CTRL_global.lons 
lon, lat=np.meshgrid([lons],[lats])


############ SEASONAL MEAN CONTROL ###########
DJF1 = rsut_CTRL_global.DJF 
MAM1 = rsut_CTRL_global.MAM
JJA1 = rsut_CTRL_global.JJA
SON1 = rsut_CTRL_global.SON


############ SEASONAL MEAN PERTRUBED ###########
DJF2 = rsut_PERT_global.DJF
MAM2 = rsut_PERT_global.MAM
JJA2 = rsut_PERT_global.JJA
SON2 = rsut_PERT_global.SON




DJF_diff = DJF2 - DJF1
MAM_diff = MAM2 - MAM1
JJA_diff = JJA2 - JJA1
SON_diff = SON2 - SON1

############ 60 to 90N MEAN  ######################################

sub_60_90N_CTRL = rsut_CTRL_global.sub_60_90N
sub_60_90N_PERT = rsut_PERT_global.sub_60_90N

###############################################################################
# This is needed to carry out significance test for seasonal mean data 
DJF_PERT = rsut_PERT_global.DJF_s
MAM_PERT = rsut_PERT_global.MAM_s
JJA_PERT = rsut_PERT_global.JJA_s
SON_PERT = rsut_PERT_global.SON_s

DJF_CTRL = rsut_CTRL_global.DJF_s
MAM_CTRL = rsut_CTRL_global.MAM_s
JJA_CTRL = rsut_CTRL_global.JJA_s
SON_CTRL = rsut_CTRL_global.SON_s

############## STATISTICAL SIGNIFICANCE TEST FOR SEASONAL DISTRIBUTION ##################################
############### DJF #########################################

##### F-TEST         ###############################

F_statics_DJF = np.zeros((64,128))
F_value_DJF = np.zeros((64,128))


for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        F_statics_DJF[i,j], F_value_DJF[i,j] = scipy.stats.levene(DJF_PERT[:,i,j], DJF_CTRL[:,i,j])
#        F_statics[i,j], F_value[i,j] = scipy.stats.bartlett(RF_ANNUAL_TOA_spdn[:,i,j], RF_ANNUAL_TOA_spdn_pertb[:,i,j])

################  T -test ##################
           
for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (F_value_DJF[i,j] <= ALPHA):
            T_value_DJF, P_value_DJF = scipy.stats.mstats.ttest_ind(DJF_PERT,DJF_CTRL, equal_var=True)
            
        else:
            T_value_DJF, P_value_DJF = scipy.stats.mstats.ttest_ind(DJF_PERT,DJF_CTRL, equal_var=False)  


# to plot 1 and 0 for statistical significance test   

P_DJF = np.zeros((64,128))

for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (abs(P_value_DJF[i,j]) <= ALPHA):
            P_DJF[i,j] = 1   #"There are significant difference between two simulations."
        else:
            P_DJF[i,j] = 0  # "There are no significant difference between two simulations."



################### MAM ######################################

F_statics_MAM = np.zeros((64,128))
F_value_MAM = np.zeros((64,128))


for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        F_statics_MAM[i,j], F_value_MAM[i,j] = scipy.stats.levene(MAM_PERT[:,i,j], MAM_CTRL[:,i,j])
#        F_statics[i,j], F_value[i,j] = scipy.stats.bartlett(RF_ANNUAL_TOA_spdn[:,i,j], RF_ANNUAL_TOA_spdn_pertb[:,i,j])

################  T -test ##################
           
for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (F_value_MAM[i,j] <= ALPHA):
            T_value_MAM, P_value_MAM = scipy.stats.mstats.ttest_ind(MAM_PERT,MAM_CTRL, equal_var=True)
            
        else:
            T_value_MAM, P_value_MAM = scipy.stats.mstats.ttest_ind(MAM_PERT,MAM_CTRL, equal_var=False)  


# to plot 1 and 0 for statistical significance test   

P_MAM = np.zeros((64,128))

for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (abs(P_value_MAM[i,j]) <= ALPHA):
            P_MAM[i,j] = 1   #"There are significant difference between two simulations."
        else:
            P_MAM[i,j] = 0  # "There are no significant difference between two simulations."


#################################  JJA ############################################

F_statics_JJA = np.zeros((64,128))
F_value_JJA = np.zeros((64,128))


for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        F_statics_JJA[i,j], F_value_JJA[i,j] = scipy.stats.levene(JJA_PERT[:,i,j], JJA_CTRL[:,i,j])
#        F_statics[i,j], F_value[i,j] = scipy.stats.bartlett(RF_ANNUAL_TOA_spdn[:,i,j], RF_ANNUAL_TOA_spdn_pertb[:,i,j])

################  T -test ##################
           
for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (F_value_JJA[i,j] <= ALPHA):
            T_value_JJA, P_value_JJA = scipy.stats.mstats.ttest_ind(JJA_PERT,JJA_CTRL, equal_var=True)
            
        else:
            T_value_JJA, P_value_JJA = scipy.stats.mstats.ttest_ind(JJA_PERT,JJA_CTRL, equal_var=False)  


# to plot 1 and 0 for statistical significance test   

P_JJA = np.zeros((64,128))

for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (abs(P_value_JJA[i,j]) <= ALPHA):
            P_JJA[i,j] = 1   #"There are significant difference between two simulations."
        else:
            P_JJA[i,j] = 0  # "There are no significant difference between two simulations."


######################## SON ###########################################


F_statics_SON = np.zeros((64,128))
F_value_SON = np.zeros((64,128))


for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        F_statics_SON[i,j], F_value_SON[i,j] = scipy.stats.levene(SON_PERT[:,i,j], SON_CTRL[:,i,j])
#        F_statics[i,j], F_value[i,j] = scipy.stats.bartlett(RF_ANNUAL_TOA_spdn[:,i,j], RF_ANNUAL_TOA_spdn_pertb[:,i,j])

################  T -test ##################
           
for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (F_value_SON[i,j] <= ALPHA):
            T_value_SON, P_value_SON = scipy.stats.mstats.ttest_ind(SON_PERT,SON_CTRL, equal_var=True)
            
        else:
            T_value_SON, P_value_SON = scipy.stats.mstats.ttest_ind(SON_PERT,SON_CTRL, equal_var=False)  


# to plot 1 and 0 for statistical significance test   

P_SON = np.zeros((64,128))

for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (abs(P_value_SON[i,j]) <= ALPHA):
            P_SON[i,j] = 1   #"There are significant difference between two simulations."
        else:
            P_SON[i,j] = 0  # "There are no significant difference between two simulations."


##### MEAN ANNUAL CONTROL #########
annual_mean_temp1 = rsut_CTRL_global.ANN_MEAN
mean_annual1 = np.average(annual_mean_temp1,axis=0)



##### MEAN ANNUAL PERTURBED #########
annual_mean_temp2 = rsut_PERT_global.ANN_MEAN
mean_annual2 = np.average(annual_mean_temp2,axis=0)

#####################################

mean_annual_diff = mean_annual2 - mean_annual1


##### F-TEST         ###############################

F_statics = np.zeros((64,128))
F_value = np.zeros((64,128))


for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        F_statics[i,j], F_value[i,j] = scipy.stats.levene(annual_mean_temp1[:,i,j], annual_mean_temp2[:,i,j])
#        F_statics[i,j], F_value[i,j] = scipy.stats.bartlett(RF_ANNUAL_TOA_spdn[:,i,j], RF_ANNUAL_TOA_spdn_pertb[:,i,j])

################  T -test ##################

mm = annual_mean_temp1
nn = annual_mean_temp2
           
for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (F_value[i,j] <= ALPHA):
            T_value, P_value = scipy.stats.mstats.ttest_ind(mm,nn, equal_var=True)
            
        else:
            T_value, P_value = scipy.stats.mstats.ttest_ind(mm,nn, equal_var=False)  


# to plot 1 and 0 for statistical significance test   

P_value_ANNUAL_sig_test_T = np.zeros((64,128))

for i in my_range(0,63,1):
    for j in my_range(0,127,1):
        if (abs(P_value[i,j]) <= ALPHA):
            P_value_ANNUAL_sig_test_T[i,j] = 1   #"There are significant difference between two simulations."
        else:
            P_value_ANNUAL_sig_test_T[i,j] = 0  # "There are no significant difference between two simulations."





##############################################################################
lon_x=rsut_PERT_global.lon_x
lat_y=rsut_PERT_global.lat_y
xx, yy = m_golbe(lon_x, lat_y)





############### PLOTTING BEGINS HERE ##########################################

###############################################################################

fig=plt.figure(4001,figsize=(13,9))
fig.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.1,hspace=0.1)
ax = fig.add_subplot(221)
x, y = m_golbe(lon, lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
#cax = m_golbe.pcolor(x,y,DJF_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r,vmin=-1.,vmax=1.)
#cax = m_golbe.pcolor(x,y,DJF_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r)
cax = m_golbe.pcolor(x,y,DJF_diff,shading='flat',cmap=my_cmap,vmin=-30.,vmax=30. )
m_golbe.scatter(x,y,P_DJF,marker='D',color='k',linewidth=0.25,vmin=1,vmax=1)
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
cbar.set_ticks([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
cbar.ax.set_yticklabels([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])

plt.text(0.50, 1.15,'RSUT Difference: CanESM4.1 [Perturb - Control]', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'a) DJF', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

############## ############ ##########
ax = fig.add_subplot(222)
my_cmap = matplotlib.cm.get_cmap('rainbow')
#cax = m_golbe.pcolor(x,y,MAM_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r,vmin=-1.,vmax=1.)
#cax = m_golbe.pcolor(x,y,MAM_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r)
cax = m_golbe.pcolor(x,y,MAM_diff,shading='flat',cmap=my_cmap,vmin=-30.,vmax=30. )
m_golbe.scatter(x,y,P_MAM,marker='D',color='k',linewidth=0.25,vmin=1,vmax=1)
m_golbe.drawcoastlines(linewidth=1.0)
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
cbar.set_ticks([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
cbar.ax.set_yticklabels([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])

#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'b) MAM', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')



############## ############ ##########
ax = fig.add_subplot(223)
my_cmap = matplotlib.cm.get_cmap('rainbow')
#cax = m_golbe.pcolor(x,y,JJA_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r,vmin=-1.,vmax=1.)
#cax = m_golbe.pcolor(x,y,JJA_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r)
cax = m_golbe.pcolor(x,y,JJA_diff,shading='flat',cmap=my_cmap,vmin=-30.,vmax=30. )
m_golbe.scatter(x,y,P_JJA,marker='D',color='k',linewidth=0.25,vmin=1,vmax=1)
m_golbe.drawcoastlines(linewidth=1.0)
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
cbar.set_ticks([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
cbar.ax.set_yticklabels([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])

#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'c) JJA', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')


############## ############ ##########
ax = fig.add_subplot(224)
my_cmap = matplotlib.cm.get_cmap('rainbow')
#cax = m_golbe.pcolor(x,y,SON_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r,vmin=-1.,vmax=1.)
#cax = m_golbe.pcolor(x,y,SON_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r)
cax = m_golbe.pcolor(x,y,SON_diff,shading='flat',cmap=my_cmap,vmin=-30.,vmax=30. )
m_golbe.scatter(x,y,P_SON,marker='D',color='k',linewidth=0.25,vmin=1,vmax=1)
m_golbe.drawcoastlines(linewidth=1.0)
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
cbar.set_ticks([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
cbar.ax.set_yticklabels([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])

#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'d) SON', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')


plt.savefig('Figure_101_rsut_g.png',dpi=150)

##############################################################################

fig=plt.figure(4002,figsize=(30,15))
fig.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.2,hspace=0.1)
ax = fig.add_subplot(221)
x, y = m_golbe(lon, lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
#cax = m_golbe.pcolor(x,y,mean_annual_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r,vmin=-1.,vmax=1.)
#cax = m_golbe.pcolor(x,y,mean_annual_diff,shading='flat',cmap=colorbar_cmap.my_cmap_bias_r)
cax = m_golbe.pcolor(x,y,mean_annual_diff,shading='flat',cmap=my_cmap,vmin=-30.,vmax=30. )
#m_golbe.scatter(xx,yy,P_value_ANNUAL_sig_test_T,marker='D',color='k',linewidth=2.,vmin=1,vmax=1)
m_golbe.scatter(x,y,P_value_ANNUAL_sig_test_T,marker='D',color='k',linewidth=0.5,vmin=1,vmax=1)
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
cbar.set_ticks([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
cbar.ax.set_yticklabels([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])

plt.text(0.90, 1.19,'Statistical Distribution', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.05,'a) Mean annual difference (Outgoing SW Radiation)', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

###########################################################################
ax = fig.add_subplot(223)
x, y = m_golbe(lon, lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
cax = m_golbe.pcolor(x,y,P_value,shading='flat',cmap=my_cmap,vmin=0.,vmax=1.)
#cax = m.pcolor(x,y,P_value,shading='flat',cmap=my_cmap_bias_r)
m_golbe.drawcoastlines(linewidth=1.0)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
#figure_title = '1 - Ratio GE 0.75  \
#                0 - Ratio LT 0.75'
figure_title = ' '

ax2  = plt.subplot(2,2,3)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)

#plt.text(0.90, 1.19,'Statistical Distribution', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
plt.text(0.0, 1.03,'b) P-value', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

############## ############ ##########
ax = fig.add_subplot(224)
my_cmap = matplotlib.cm.get_cmap('rainbow')
cax = m_golbe.pcolor(x,y,P_value_ANNUAL_sig_test_T,shading='flat',cmap=colorbar_cmap.my_cmap_stat,vmin=0,vmax=1)
#cax = m.pcolor(x,y,MAM_diff,shading='flat',cmap=my_cmap_bias_r)
m_golbe.drawcoastlines(linewidth=1.0)
#m_golbe.drawcountries(linewidth=1.)
#m_golbe.drawstates(linewidth=1., linestyle='solid', color='k')
m_golbe.drawparallels(np.arange(-90.,91.,30.))
m_golbe.drawmeridians(np.arange(0.,360.,60.))
plt.gcf().subplots_adjust(wspace=0.1,hspace=0.19)
figure_title = '1 - Siginificant \
                0 - Not Siginificant'

ax2  = plt.subplot(2,2,4)       # this code increase gap between title and figure 
plt.text(0.5, 1.08, figure_title,
         horizontalalignment='center',
         fontsize=15,
         transform = ax2.transAxes)

cbar = m_golbe.colorbar(cax, pad = 0.2, location='right',filled=True)
cbar.set_ticks([0, 1])
#plt.text(0.65, 1.3,'Surface temperature', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=25,color='k')
#plt.text(0.0, 1.05,'b) Significance at ${\\alpha}$ = ', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

plt.text(0.0, 1.03,'c) Significance at ${\\alpha}$ = ', horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')
plt.text(0.65, 1.03,ALPHA, horizontalalignment='left',verticalalignment='center', transform=ax2.transAxes,fontsize=20,color='k')

plt.savefig('Figure_102_rsut_g.png',dpi=150)



################################################

smoothed_TS1 = rsut_CTRL_global.smoothed_TS
TS_subdom1 = rsut_CTRL_global.TS_subdom
TS_mean1 = "%.2f" % round(TS_subdom1.mean(axis=0),2)
str_leg1 = ''.join(['10-yrs moving average (Control run, Mean =', TS_mean1, ')']) 

smoothed_TS2 = rsut_PERT_global.smoothed_TS
TS_subdom2 = rsut_PERT_global.TS_subdom
TS_mean2 = "%.2f" % round(TS_subdom2.mean(axis=0),2)
str_leg2 = ''.join(['10-yrs moving average (Canadian shipping x 10$^{6}$, Mean =', TS_mean2, ')']) 

YEAR = rsut_CTRL_global.YEAR
YEAR_moving_avg = rsut_CTRL_global.YEAR_moving_avg


fig=plt.figure(4003,figsize=(10,5))
fig.subplots_adjust(left=None,bottom=0.25,right=None,top=None,wspace=0.2,hspace=0.1)
plt.plot(YEAR_moving_avg, smoothed_TS1, linewidth=2, color='darkgreen',label= str_leg1)
plt.plot(YEAR_moving_avg, smoothed_TS2, linewidth=2, color='darkblue',label= str_leg2)
plt.plot(YEAR, TS_subdom1, linewidth=0.5, color='k',label='Downward flux (averaged over the sub-domain)')
plt.plot(YEAR, TS_subdom2, linewidth=0.5, color='k')

ax = plt.gca()
plt.ylabel('Outgoing SW Radiation at TOA [W/m$^{2}$]',fontsize=18, fontweight='normal')
plt.xlabel('Model year',fontsize=18, fontweight='normal')
plt.tick_params(which='major', width=2, length=7)
plt.tick_params(direction='out', top=False, right=False) # Turn ticks out
plt.legend(loc='upper left',frameon=False,fontsize='medium',ncol=1)

ax.set_xlim([6001,6100])
ax.set_xticks(np.linspace(6000,6100,6))
plt.tick_params(which='major', width=2, length=7)
plt.tick_params(direction='out', top=False, right=False) # Turn ticks out
plt.legend(loc='upper right',frameon=False,fontsize='medium',ncol=1)
#plt.gcf().subplots_adjust(bottom=0.25) 

plt.savefig('Figure_103_rsut_g.png',dpi=150)


############# 60 - 90N mean weighted ##########

sub_mean1 = "%.2f" % round(sub_60_90N_CTRL.mean(axis=0),2)
str_leg1 = ''.join(['Control run (Mean = ', sub_mean1, ')']) 


sub_mean2 = "%.2f" % round(sub_60_90N_PERT.mean(axis=0),2)
str_leg2 = ''.join(['Canadian shipping x 10$^{6}$ (Mean = ', sub_mean2, ')']) 

YEAR_sub = YEAR = np.linspace(6071,6100,30)


fig=plt.figure(4004,figsize=(10,5))
fig.subplots_adjust(left=None,bottom=0.25,right=None,top=None,wspace=0.2,hspace=0.1)
#plt.plot(YEAR_moving_avg, smoothed_TS1, linewidth=2, color='darkgreen',label= str_leg1)
#plt.plot(YEAR_moving_avg, smoothed_TS2, linewidth=2, color='darkblue',label= str_leg2)
plt.plot(YEAR_sub, sub_60_90N_CTRL, linewidth=2, color='darkgreen',label=str_leg1)
plt.plot(YEAR_sub, sub_60_90N_PERT, linewidth=2, color='darkblue',label=str_leg2)

ax = plt.gca()
plt.ylabel('Outgoing SW Radiation at TOA [W/m$^{2}$]',fontsize=18, fontweight='normal')
plt.xlabel('Model year',fontsize=18, fontweight='normal')
plt.tick_params(which='major', width=2, length=7)
plt.tick_params(direction='out', top=False, right=False) # Turn ticks out
plt.legend(loc='upper left',frameon=False,fontsize='medium',ncol=1)

ax.set_xlim([6071,6100])
ax.set_xticks(np.linspace(6070,6100,7))
plt.tick_params(which='major', width=2, length=7)
plt.tick_params(direction='out', top=False, right=False) # Turn ticks out
plt.legend(loc='center right',frameon=False,fontsize='medium',ncol=1)
#plt.gcf().subplots_adjust(bottom=0.25) 
plt.text(0.40, 1.05,'60$^{o}$N - 90$^{o}$N', horizontalalignment='left',verticalalignment='center', transform=ax.transAxes,fontsize=18,color='k')
plt.savefig('Figure_104_rsut_60-90N.png',dpi=150)







