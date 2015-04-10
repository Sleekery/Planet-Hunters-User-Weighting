from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter
from matplotlib import rcParams

count=0
iterange=range(1413,3793)
#iterange=[2258]
for j in iterange:
	syntheticid=str(j)
	#oneclassid='545fd8b4415ac156fe000908'  # Select one specific classid to look at
	oneclassid=''                          # Not that
	if len(oneclassid)>0:                       # Checks whether correct syntheticid was
		with open('matching.dat') as infile:    # chosen and if not, switches to the 
			for line in infile:                 # syntheticid that corresponds to the given 
				linesplit=line.split('&')       # classid
				if linesplit[29]==oneclassid:
					if syntheticid==linesplit[11].split('_')[0]:
						synidprint='syntheticid correctly corresponds to given classid'
					else:
						synidprint='syntheticid changed to correct value for given classid'
						syntheticid=linesplit[11].split('_')[0]
					break
	
	# Extracing the time and flux from the given syntheticid
	ffits=fits.open('./Mdwarfsynthetics/synthetic_'+syntheticid+'.fits')  
	time=ffits[1].data['TIME']
	flux=ffits[1].data['PDCSAP_FLUX']
	
	# Extracting important information from 'matching.'dat'
	i=0
	transitid,userxmin,userxmax,synxmin,synxmax,classid=[],[],[],[],[],[]
	synpixmin,synpixmax=[],[]
	with open('matching.dat') as infile:
		for line in infile:
			linesplit=line.split('&')
			if linesplit[11].split('_')[0] == syntheticid:
				transitid.append(linesplit[11])
				userxmin.append(float(linesplit[20]))
				userxmax.append(float(linesplit[21]))
				synxmin.append(float(linesplit[16]))
				synxmax.append(float(linesplit[17]))
				synpixmin.append(float(linesplit[14]))
				synpixmax.append(float(linesplit[15]))
				classid.append(linesplit[29])
				period='{0:0.2f}'.format(float(linesplit[6]))
				prad='{0:0.2f}'.format(float(linesplit[7]))
				srad='{0:0.2f}'.format(float(linesplit[8]))
				kepmag='{0:0.2f}'.format(float(linesplit[9]))
				kepid=linesplit[0]
			if i > np.inf:
				break
			i+=1
	
	
	# Setting up figure stuff
	ms=5
	grid = gs.GridSpec(2,1,height_ratios=[2.5,1],hspace=0)
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
	
	tolerance=0.01  # Necessary to get the end points of the marked transits
	setsynxmin=list(set(synxmin))  # Getting the unique set of synxmins
	setsynxmax=list(set(synxmax))  # Getting the unique set of synxmaxes
	setsynxmin.sort()
	setsynxmax.sort()
	timesyn,fluxsyn=np.asarray([]),np.asarray([])
	for i in range(len(setsynxmin)):  # Finding all points in all synthetic transists
		timesyn=np.append(timesyn,time[np.where((time>=(setsynxmin[i]-tolerance)) & (time <= (	setsynxmax[i]+tolerance)))])
		fluxsyn=np.append(fluxsyn,flux[np.where((time>=(setsynxmin[i]-tolerance)) & (time <= (	setsynxmax[i]+tolerance)))])
	
	# Plotting all data an then overplotting the synthetic poitns
	ax1.plot(time,flux,ls='none',marker='o',mec='none',color='k')
	ax1.plot(timesyn,fluxsyn,ls='none',marker='o',mec='none',color=syncolor)
	
	# If I want to look at one specific classid, this finds all points within the user markings
	if len(oneclassid)>0:  
		userxmin=np.asarray(userxmin)[np.where(np.asarray(classid)==oneclassid)]
		userxmin=userxmin[np.where(userxmin>0)]
		userxmax=np.asarray(userxmax)[np.where(np.asarray(classid)==oneclassid)]
		userxmax=userxmax[np.where(userxmax>0)]
		timeuser,fluxuser=np.asarray([]),np.asarray([])
		for i in range(len(userxmin)):
			timeuser=np.append(timeuser,time[np.where((time>=(userxmin[i]-tolerance)) & (time <=	 (userxmax[i]+tolerance)))])
			fluxuser=np.append(fluxuser,flux[np.where((time>=(userxmin[i]-tolerance)) & (time <=	 (userxmax[i]+tolerance)))])
		ax.plot(timeuser,fluxuser,ls='none',marker='o',mec='none',color=usercolor,markersize=ms)
		print synidprint
	
	# Finding all classifications of the given syntheticid
	userxmin,userxmax,allclassid,markedclassid=[],[],[],[]
	i=0
	with open('2015-02-03_planet_hunter_classifications.csv') as infile:
	    for line in infile:
	    	tmpline=line.replace('\",\"','&')[1:-1]
	    	linesplit=tmpline.split('&')
	    	if linesplit[10] == syntheticid:
	    		allclassid.append(linesplit[0])
	    	if linesplit[10] == syntheticid and len(linesplit[14])>0:
	 	    	markedclassid.append(linesplit[0])
	    		userxmin.append(float(linesplit[14]))
	    		userxmax.append(float(linesplit[15][:-1]))
	    	if i > np.inf:
	    		break
	    	i+=1
	# Finding all points marked by any user and avoiding where users marked the same data points twice
	timeuser,fluxuser=np.asarray([]),np.asarray([])
	setmarkedclassid=list(set(markedclassid))
	for i in range(len(setmarkedclassid)):
		tmptimeuser,tmpfluxuser=np.asarray([]),np.asarray([])
		for k in range(len(markedclassid)):
			if markedclassid[k]==setmarkedclassid[i]:
				tmptimeuser=np.append(tmptimeuser,time[np.where((time>=(userxmin[k]-tolerance)) & (time <= (userxmax[k]+tolerance)))])
				tmpfluxuser=np.append(tmpfluxuser,flux[np.where((time>=(userxmin[k]-tolerance)) & (time <= (userxmax[k]+tolerance)))])
		timeuser=np.append(timeuser,np.asarray(list(set(tmptimeuser))))
		fluxuser=np.append(fluxuser,np.asarray(list(set(tmpfluxuser))))
	
	#ax.plot(timeuser,fluxuser,ls='none',marker='o',mec='none',color='r',markersize=ms)
	# Filling in the np.nans in the time array for the histogram
	histtime=np.append(time,time[-1]+tolerance)
	for i in range(len(histtime)):
		if histtime[i] > 0:
			pass
		elif i > 1:
			histtime[i]=histtime[i-1]+(histtime[i-1]-histtime[i-2])
		else: 
			histtime[i]=histtime[i-1]+0.0204338629846
	
	hist=np.histogram(timeuser,histtime)
	ax2.plot(hist[1][:-1],hist[0],color=usercolor)
	ax2.fill_between(hist[1][:-1], 0, hist[0],color=usercolor)
	
	# Setting up the grid for the heat map
	y0,y1=ax1.get_ylim()[0],ax1.get_ylim()[1]
	dy=ax1.get_ylim()[1]-ax1.get_ylim()[0]
	x0,x1=histtime[0],histtime[-2]
	dx=histtime[1]-histtime[0]
	y, x = np.mgrid[slice(y0,y1+dy,dy),slice(x0,x1+dx,dx)]
	x,y=x[0:2],y[0:2]
	z=np.asarray([list(hist[0])])
	z_min, z_max = 0, z.max()
	ax1.pcolor(x, y, z, cmap=cmap,vmin=z_min,vmax=z_max)
	ax2.yaxis.grid(True)
	
	# Setting up the limits and ticks
	xlimoffset=0.2
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
	def topercent(y, position):
	    # Ignore the passed in position. This has the effect of scaling the default
	    # tick locations.
	    s = str(100 * y)[:-2]
	    # The percent symbol needs escaping in latex
	    if rcParams['text.usetex'] == True:
	        return s + r'$\%$'
	    else:
	        return s + '%'
	
	ax22 = ax2.twinx()
	ax22.set_ylim(0,ax2.get_ylim()[1]/float(len(set(allclassid))))
	ax22.set_ylabel('Percent of all\nclassifications',rotation=270,labelpad=30)
	y2formatter = FuncFormatter(topercent)
	ax22.yaxis.set_major_formatter(y2formatter)
	ax22.yaxis.set_major_locator(MaxNLocator(prune='lower'))

	ax11 = ax1.twinx()
	ax11.set_ylabel(kepid+'\nP='+period+'\nRp='+prad+'\nRs='+srad+'\nmag='+kepmag,rotation=0,fontsize=9,labelpad=30)
	ax11.axes.get_yaxis().set_ticks([])

	#plt.show()
	fig.savefig('./heatmap/'+syntheticid+'.eps')
	plt.close(fig)
	count+=1
	print '{0:0.2f}'.format(100.0*float(count)/float(len(iterange)))+'% completed. Syntheticid = '+syntheticid+'.'












