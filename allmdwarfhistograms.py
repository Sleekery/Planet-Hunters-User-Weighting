from astropy.io import ascii
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter
from matplotlib import rcParams

plot=False
writethefile=True

def topercent(y, position):
   # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)[:-2]
    # The percent symbol needs escaping in latex
    if rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'

# Obtaining list of M-dwarf stars in the S-giants, >=13 Q's of data
mdwarfkepid=[]
i=0
with open('phmdwarfs.csv') as infile:
	for line in infile:
		if i>0:
			mdwarfkepid.append(int(line.split(',')[2]))
		i+=1

# Obtaining list of M-dwarf stars in the S-giants, <13 Q's of data 
i=0
with open('phmdwarfsLT13q.csv') as infile:
	for line in infile:
		if i>0:
			mdwarfkepid.append(int(line.split(',')[2]))
		i+=1

tolerance=0.01
grid = gs.GridSpec(2,1,height_ratios=[2.5,1],hspace=0)
count=0
mdwarfkepid=mdwarfkepid[3642:]
# NEED TO DO kepid=2439250, i=3641
for i in range(len(mdwarfkepid)):
	# Have to take into account how there are multiple quarters for each Mdwarf
	# Need to add 0's to the filename because all filenames are 9 characters long
	if len(str(mdwarfkepid[i]))==6:
		filelist=glob.glob('./lightcurvedata/'+'000'+str(mdwarfkepid[i])+'/*_llc.fits')
	elif len(str(mdwarfkepid[i]))==7:
		filelist=glob.glob('./lightcurvedata/'+'00'+str(mdwarfkepid[i])+'/*_llc.fits')
	elif len(str(mdwarfkepid[i]))==8:
		filelist=glob.glob('./lightcurvedata/'+'0'+str(mdwarfkepid[i])+'/*_llc.fits')
	else:
		print 'Length of kepid != 6, 7, or 8'
		break 
	filelist.sort()
	# Creating a list of all times and fluxes listed in the FITS files
	time,flux,plottime=[],[],[]
	for j in range(len(filelist)):
		ffits=fits.open(filelist[j])
		time.extend(ffits[1].data['TIME'])
		flux.extend(ffits[1].data['PDCSAP_FLUX'])
		# Creating a list of all times with gaps filled in
		if j == 0:
			plottime.extend(ffits[1].data['TIME'])
		else:
			plottime.extend(np.arange(plottime[-1]+0.0204338629846,ffits[1].data['TIME'][0]-0.015,0.0204338629846))
			plottime.extend(ffits[1].data['TIME'])
	time=np.asarray(time)
	flux=np.asarray(flux)
	plottime=np.asarray(plottime)
	# Finding all classifications of the given mdwarfkepid
	userxmin,userxmax,allclassid,markedclassid=[],[],[],[]
	j=0
	with open('2015-02-03_planet_hunter_classifications.csv') as infile:
	    for line in infile:
	    	tmpline=line.replace('\",\"','&')[1:-1]
	    	linesplit=tmpline.split('&')
	    	#print linesplit[7], mdwarfkepid[i]
	    	if linesplit[7] == str(mdwarfkepid[i]):
	    		allclassid.append(linesplit[0])
	    	if linesplit[7] == str(mdwarfkepid[i]) and len(linesplit[14])>0:
	    		markedclassid.append(linesplit[0])
	    		userxmin.append(float(linesplit[14]))
	    		userxmax.append(float(linesplit[15][:-1]))
	    	if j > np.inf:
	    		break
	    	j+=1		
	# Finding all points marked by any user and avoiding where users marked the same data points twice
	timeuser,fluxuser=np.asarray([]),np.asarray([])
	setmarkedclassid=list(set(markedclassid))
	for j in range(len(setmarkedclassid)):
		tmptimeuser,tmpfluxuser=np.asarray([]),np.asarray([])
		for k in range(len(markedclassid)):
			if markedclassid[k]==setmarkedclassid[j]:
				tmptimeuser=np.append(tmptimeuser,time[np.where((time>=(userxmin[k]-tolerance)) & (time <= (userxmax[k]+tolerance)))])
				tmpfluxuser=np.append(tmpfluxuser,flux[np.where((time>=(userxmin[k]-tolerance)) & (time <= (userxmax[k]+tolerance)))])
		timeuser=np.append(timeuser,np.asarray(list(set(tmptimeuser))))
		fluxuser=np.append(fluxuser,np.asarray(list(set(tmpfluxuser))))
	# Filling in the np.nans in the time array for the histogram
	histtime=np.append(time,time[-1]+tolerance)
	for j in range(len(histtime)):
		if histtime[j] > 0:
			pass
		elif j > 1:
			histtime[j]=histtime[j-1]+(histtime[j-1]-histtime[j-2])
		else: 
			histtime[j]=histtime[j-1]+0.0204338629846
	hist=np.histogram(timeuser,histtime)
	times=list(hist[1][:-1])
	histcount=list(hist[0])
	if writethefile==True:
		filewrite=open('./allmdwarfhistograms/'+str(mdwarfkepid[i])+'.dat','w')
		filewrite.write(str(mdwarfkepid[i])+'&'+str(len(set(allclassid)))+'&'+str(times)[1:-1]+'&'+str(histcount)[1:-1]+' \n')
		filewrite.close()
	if plot==True:
		ms=5
		fig=plt.figure()
		fig.set_size_inches(6.75,5.0)
		fig.subplots_adjust(right=0.875,top=0.995,left=0.145,bottom=0.095)
		ax1=fig.add_subplot(grid[0])
		ax1.set_ylabel('Normalized flux')
		ax2=fig.add_subplot(grid[1])
		ax2.set_xlabel('BJD-2454833 (days)')
		ax2.set_ylabel('Count')
		#ax1.set_title('Synthetic ID = '+syntheticid)
		# Setting up the color scheme
		syncolor='r'
		usercolor='b'
		if usercolor=='r':
			cmap='Reds'
		elif usercolor=='b':
			cmap='Blues'
		# Setting up the grid for the heat map
		ax1.plot(time,flux,ls='none',marker='o',mec='none',color='k')
		plothisttime=np.append(plottime,plottime[-1]+tolerance)
		y0,y1=ax1.get_ylim()[0],ax1.get_ylim()[1]
		dy=ax1.get_ylim()[1]-ax1.get_ylim()[0]
		x0,x1=plothisttime[0],plothisttime[-2]
		dx=plothisttime[1]-plothisttime[0]
		y, x = np.mgrid[slice(y0,y1+dy,dy),slice(x0,x1+dx,dx)]
		x,y=x[0:2],y[0:2]
		plothist=np.histogram(timeuser,x[0])
		z=np.asarray([list(plothist[0])])
		print len(x[0]),len(y[0]),len(z[0])
		z_min, z_max = 0, z.max()
		ax1.pcolor(x, y, z, cmap=cmap,vmin=z_min,vmax=z_max)
		ax2.plot(plothist[1][:-1],plothist[0],color=usercolor)
		ax2.fill_between(plothist[1][:-1], 0, plothist[0],color=usercolor)
		ax2.yaxis.grid(True)
		# Setting up the limits and ticks
		xlimoffset=1.0
		ax1.set_xlim(np.nanmin(time)-xlimoffset,np.nanmax(time)+xlimoffset)
		ax2.set_xlim(np.nanmin(time)-xlimoffset,np.nanmax(time)+xlimoffset)
		ax1.set_ylim(y0,y1)
		ax1.axes.get_xaxis().set_ticks([])
		y1formatter = ScalarFormatter(useOffset=False)
		ax1.yaxis.set_major_formatter(y1formatter)
		ax1.yaxis.set_major_locator(MaxNLocator(prune=None))
		ticklabels=ax1.get_yticks().tolist()
		ax1.yaxis.set_major_locator(MaxNLocator(prune='both'))
		ax2.yaxis.set_major_locator(MaxNLocator(prune='both'))
		# Setting up percentage ylabel on twin axis		
		ax22 = ax2.twinx()
		ax22.set_ylim(0,ax2.get_ylim()[1]/float(len(set(allclassid))))
		ax22.set_ylabel('Percent of all\nclassifications',rotation=270,labelpad=30)
		y2formatter = FuncFormatter(topercent)
		ax22.yaxis.set_major_formatter(y2formatter)
		ax22.yaxis.set_major_locator(MaxNLocator(prune='lower'))
		ax11 = ax1.twinx()
		ax11.set_ylabel(str(mdwarfkepid[i]),rotation=0	,fontsize=9,labelpad=30)
		ax11.axes.get_yaxis().set_ticks([])
		plt.show()
	count+=1
	print '{0:0.2f}'.format(100.0*float(count)/float(len(mdwarfkepid)))+'% completed. Kepid = '+str(mdwarfkepid[i])+'. i = '+str(i)+'.'



















