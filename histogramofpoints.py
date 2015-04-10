from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter
from matplotlib import rcParams

# Pick which user weighting file/scheme to use
userweightfile='80percentuserweighting0.5comp.dat'
userweights=ascii.read(userweightfile,data_start=0,delimiter='&',names=['username','upweight','downweight','combined','normupweight','normdownweight','normcombined','numclasses'])

# The random synthetics test set
count=0
testrandomsynthetics=True
weighting=True
if testrandomsynthetics==True:
	randomsynthetics=ascii.read('80randomsynthetics.dat',data_start=0,names=['syntheticid'])
	matching=ascii.read('matching.dat',header_start=0,data_start=1,delimiter='&')
	allsynthetics=list(set(matching['syntheticid']))
	# Create list of all non-tested synthetics
	nontestedtransits=[]
	for i in range(len(allsynthetics)):
		if allsynthetics[i] not in randomsynthetics['syntheticid']:
			nontestedtransits.append(int(allsynthetics[i].split('_')[0]))
	nontestedtransits=list(set(nontestedtransits))
	nontestedtransits.sort()
	if weighting==True:
		filewrite=open('80weightedhistogramofpoints0.5comp.dat','w')
	else:
		filewrite=open('80nonweightedhistogramofpoints.dat','w')
	iterange=nontestedtransits
else:
	if weighting==True:
		filewrite=open('weightedhistogramofpoints.dat','w')
	else:
		filewrite=open('nonweightedhistogramofpoints.dat','w')
	iterange=range(1413,3793)

for j in iterange:
	syntheticid=str(j) # The entire light curve
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
				transitid.append(linesplit[11]) # The individual transit
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
	tolerance=0.01  # Necessary to get the end points of the marked transits
	setsynxmin=list(set(synxmin))  # Getting the unique set of synxmins
	setsynxmax=list(set(synxmax))  # Getting the unique set of synxmaxes
	setsynxmin.sort()
	setsynxmax.sort()
	timesyn,fluxsyn=np.asarray([]),np.asarray([])
	for i in range(len(setsynxmin)):  # Finding all points in all synthetic transists
		timesyn=np.append(timesyn,time[np.where((time>=(setsynxmin[i]-tolerance)) & (time <= (setsynxmax[i]+tolerance)))])
		fluxsyn=np.append(fluxsyn,flux[np.where((time>=(setsynxmin[i]-tolerance)) & (time <= (setsynxmax[i]+tolerance)))])
	# Finding all classifications of the given syntheticid
	userxmin,userxmax,allclassid,markedclassid=[],[],[],[]
	allusernames,markedusernames=[],[]
	i=0
	with open('2015-02-03_planet_hunter_classifications.csv') as infile:
		for line in infile:
			tmpline=line.replace('\",\"','&')[1:-1]
			linesplit=tmpline.split('&')
			if linesplit[10] == syntheticid:
				allclassid.append(linesplit[0])
				allusernames.append(tmpline.split('&')[6])
			if linesplit[10] == syntheticid and len(linesplit[14])>0:
				markedusernames.append(tmpline.split('&')[6])
				markedclassid.append(linesplit[0])
				userxmin.append(float(linesplit[14]))
				userxmax.append(float(linesplit[15][:-1]))
			if i > np.inf:
				break
			i+=1
	# Removing the spaces at the end of some usernames that go missing when reading in files with astropy.ascii
	for i in range(len(allusernames)):
		if allusernames[i][-1]==' ':
			allusernames[i]=allusernames[i][:-1]
	for i in range(len(markedusernames)):
		if markedusernames[i][-1]==' ':
			markedusernames[i]=markedusernames[i][:-1]	
	# Finding all points marked by any user and avoiding where users marked the same data points twice
	weightedcount=np.zeros(len(time))
	setmarkedclassid=list(set(markedclassid))
	for i in range(len(setmarkedclassid)):
		timeuser=np.asarray([])
		username=allusernames[allclassid.index(setmarkedclassid[i])]
		userweight=userweights['normcombined'][np.where(userweights['username']==username)][0]
		for k in range(len(markedclassid)):
			if markedclassid[k]==setmarkedclassid[i]:
				timeuser=np.append(timeuser,time[np.where((time>=(userxmin[k]-tolerance)) & (time <= (userxmax[k]+tolerance)))])
		# Making sure each person can only increase the weightedcount for a specific time once, not multiple times in case they made multiple marks for the same time points
		timeuser=np.asarray(list(set(timeuser)))
		for i in range(len(timeuser)):
			if weighting==True:
				weightedcount[np.where(time==timeuser[i])]+=userweight
			else:
				weightedcount[np.where(time==timeuser[i])]+=1
	# Sum of weights of all unique people who saw this light curve
	totalweightedcount=0
	setallusernames=list(set(allusernames))
	for i in range(len(setallusernames)):
		userweight=userweights['normcombined'][np.where(userweights['username']==setallusernames[i])][0]
		if weighting==True:
			totalweightedcount+=userweight
		else:
			totalweightedcount+=1
	# Filling in the np.nans in the time array for the histogram
	for i in range(len(time)):
		if time[i] > 0:
			pass
		elif i > 1:
			time[i]=time[i-1]+(time[i-1]-time[i-2])
		else: 
			time[i]=time[i-1]+0.0204338629846
	time=list(time)
	weightedcount=list(weightedcount)
	filewrite.write(syntheticid+'&'+str(totalweightedcount)+'&'+str(time)[1:-1]+'&'+str(weightedcount)[1:-1]+' \n')
	count+=1
	print '{0:0.2f}'.format(100.0*float(count)/float(len(iterange)))+'% completed. Syntheticid = '+syntheticid+'.'

filewrite.close()










