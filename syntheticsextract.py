### Extract syntheticid information from each FITS file ###

from astropy.io import fits
from astropy.io import ascii
import glob as glob
import numpy as np

fwrite=open('Mdwarfstransittimes.dat','w')
#fwrite=open('tmp.dat','w')
fwrite.write('KIC'.rjust(11)+'synid'.rjust(7)+'plrad'.rjust(10)+'plper'.rjust(10)+'plphase'.rjust(10)+'synpixmin'.rjust(11)+'synpixmax'.rjust(11)+'synxmin'.rjust(21)+'synxmax'.rjust(21)+'\n')

# List of all synthetic FITS files
allfits=glob.glob('./Mdwarfsynthetics/synthetic_*.fits')

for i in range(len(allfits)):
#for i in [856]:
	ffits = fits.open(allfits[i])
	synid=int(allfits[i].split('_')[1].split('.')[0])
	synxmin,synxmax=[],[]
	# Read FITS header for synthetic information
	kic=ffits[0].header['KEPLERID']
	plrad=float(ffits[0].header['PLRAD'])
	plper=float(ffits[0].header['PLPER'])
	plphase=float(ffits[0].header['PLPHASE'])
	synpixmin=int(ffits[0].header['TRANELS'].split(':')[0])
	synpixmax=int(ffits[0].header['TRANELS'].split(':')[1])+1 # +1 fudge factor
	tmpxmin=ffits[1].data[synpixmin][0]
	tmpxmax=ffits[1].data[synpixmax][0]
	# Tests to see if synpixmin or synpixmax is a np.nan
	# If so, take the two data points n steps to the left and right, increasing n by 1 
	# each time, and check to see if they're both != np.nan
	# 
	if np.isnan(tmpxmin) == True:
		j=0
		xdown,xup=np.nan,np.nan
		while np.isnan(xdown) == True or np.isnan(xup) == True:
			j+=1
			xdown=ffits[1].data[synpixmin-j][0]
			xup=ffits[1].data[synpixmin+j][0]
		tmpxmin=np.median([xdown,xup])
	if np.isnan(tmpxmax) == True:
		j=0
		xdown,xup=np.nan,np.nan
		while np.isnan(xdown) == True or np.isnan(xup) == True:
			j+=1
			xdown=ffits[1].data[synpixmax-j][0]
			xup=ffits[1].data[synpixmax+j][0]
		tmpxmax=np.median([xdown,xup])
	# Counts how many nan's are in the transit
	nonancount=0
	for j in range(synpixmin,synpixmax+1):  
		if ffits[1].data[j][0]>0:
			nonancount+=1
	# Only counts transits with fewer than 50% nan's
	if float(nonancount)>0.5*float(len(range(synpixmin,synpixmax+1))):
		synxmin.append(tmpxmin)
		synxmax.append(tmpxmax)
	# Adds synthetic start and end times if there are multiple synthetics in the light curve after making sure they are not a majority nan's
	while tmpxmin+plper < ffits[1].data[-1][0]:
		tmpxmin=tmpxmin+plper
		if tmpxmax+plper < ffits[1].data[-1][0]:
			tmpxmax=tmpxmax+plper
		else:
			tmpxmax=ffits[1].data[-1][0]
		tmpfits=ffits[1].data[np.where(ffits[1].data['TIME']>tmpxmin)]
		tmpfits=tmpfits[np.where(tmpfits['TIME']<tmpxmax)]
		# Only counts transits with fewer than 50% nan's
		if float(len(tmpfits)) > 0.5*float(len(range(synpixmin,synpixmax+1))):
			synxmin.append(tmpxmin)
			synxmax.append(tmpxmax)
	strsynxmin=str(synxmin).replace(' ','').replace('[','').replace(']','')
	strsynxmax=str(synxmax).replace(' ','').replace('[','').replace(']','')
	fwrite.write(str(kic).rjust(11)+str(synid).rjust(7)+str(plrad).rjust(10)+str(plper).rjust(10)+str(plphase).rjust(10)+str(synpixmin).rjust(11)+str(synpixmax).rjust(11)+'   '+strsynxmin+'   '+strsynxmax+'\n')
	if i%20==0:
		print '{0:0.2f}'.format(float(i)/float(len(allfits))*100.0)+'% completed.'

fwrite.close()





