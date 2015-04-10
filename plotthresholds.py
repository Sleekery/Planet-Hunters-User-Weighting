from astropy.io import ascii
from astropy.table import Column
import numpy as np
import matplotlib.pyplot as plt

# Reading in thresholds data
samplethresholdsfile='80weightedsamplethreshold.dat'
weighteddata=ascii.read(samplethresholdsfile,header_start=0,data_start=1,delimiter='&')

samplethresholdsfile='80nonweightedsamplethreshold.dat'
nonweighteddata=ascii.read(samplethresholdsfile,header_start=0,data_start=1,delimiter='&')

samplethresholdsfile='80weightedsamplethreshold0.2comp.dat'
weighteddata02=ascii.read(samplethresholdsfile,header_start=0,data_start=1,delimiter='&')

samplethresholdsfile='80weightedsamplethreshold0.5comp.dat'
weighteddata05=ascii.read(samplethresholdsfile,header_start=0,data_start=1,delimiter='&')

fig=plt.figure()
fig.subplots_adjust(right=0.875,top=0.995,left=0.145,bottom=0.095)

# Left axis is the number of synthetic transits properly detected
ax1=fig.add_subplot(111)
ax1.plot(nonweighteddata['thresholdup'],nonweighteddata['syntces'],color='k',ls='-',label='Nonweighted')
ax1.plot(weighteddata['thresholdup'],weighteddata['syntces'],color='r',ls='--',label='0.1 weighted')
ax1.plot(weighteddata02['thresholdup'],weighteddata02['syntces'],color='g',ls='-.',label='0.2 weighted')
ax1.plot(weighteddata05['thresholdup'],weighteddata05['syntces'],color='b',ls=':',label='0.5 weighted')

# Right axis is the fraction of false positives
ax2=ax1.twinx()
ax2.plot(1,1,color='white',label='Top = Synthetics')
ax2.plot(1,1,color='white',label='Bottom = FPs')

ax2.plot(nonweighteddata['thresholdup'],nonweighteddata['fppercent'],color='k',ls='-',label='Nonweighted')
ax2.plot(weighteddata['thresholdup'],weighteddata['fppercent'],color='r',ls='--',label='0.1 weighted')
ax2.plot(weighteddata02['thresholdup'],weighteddata02['fppercent'],color='g',ls='-.',label='0.2 weighted')
ax2.plot(weighteddata05['thresholdup'],weighteddata05['fppercent'],color='b',ls=':',label='0.5 weighted')

# Getting the labels of ax2 onto ax1
#ax1.plot(1,1,color='k',ls='-',label='Nonweighted FPs')
#ax1.plot(1,1,color='r',ls='--',label='0.1 weighted FPs')
#ax1.plot(1,1,color='g',ls='-.',label='0.2 weighted FPs')
#ax1.plot(1,1,color='b',ls=':',label='0.5 weighted FPs')

# Setting the axes up
ax1.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax1.set_xlabel('Upper threshold value (fraction of users marking an event)')
ax1.set_ylabel('Number of successfully detected synthetics')
ax2.set_ylabel('Fraction of TCEs that are false positives',rotation=270,labelpad=20)
ax2.grid(True)
ax2.legend(loc=1,fontsize=12,handlelength=3.6)
ax2.set_ylim(0.0,0.76)
plt.show()

weighteddata['floatsyntces']=0.0
nonweighteddata['floatsyntces']=0.0
for i in range(len(weighteddata)):
	weighteddata['floatsyntces'][i]=float(weighteddata['syntces'][i])
	nonweighteddata['floatsyntces'][i]=float(nonweighteddata['syntces'][i])

'''
# Comparing the weighted and nonweighted values directly
fig=plt.figure()
fig.subplots_adjust(right=0.875,top=0.995,left=0.145,bottom=0.095)
ax1=fig.add_subplot(111)
ax1.plot(weighteddata['thresholdup'],weighteddata['floatsyntces']/nonweighteddata['floatsyntces'],color='b',ls='-',label='Number of synthetics')
ax2=ax1.twinx()
ax2.plot(weighteddata['thresholdup'],weighteddata['fppercent']/nonweighteddata['fppercent'],color='r',ls='-',label='False positives')
ax1.plot(1,1,color='r',ls='-',label='False positives')
ax1.set_xlabel('Upper threshold value (fraction of users marking an event)')
ax1.set_ylabel('Weighted number of success synthetic recoveries over nonweighted')
ax2.set_ylabel('Weighted fraction of false positives over nonweighted',rotation=270,labelpad=20)
ax1.legend(loc='best')
ax2.grid(True)
plt.show()
'''











