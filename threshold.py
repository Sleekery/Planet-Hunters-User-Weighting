from astropy.io import ascii
from astropy.table import Column
import numpy as np
import matplotlib.pyplot as plt

# Determines whether I sample multiple thresholds and create 'samplethresholds.dat' or whether I do it for just threshold and create 'threshold.dat'
samplethresholds=True
filetype='80weighted'
inputhistogramfile=filetype+'histogramofpoints0.5comp.dat' # Input all the time
outputthresholdfile=filetype+'threshold.dat'  # Output when samplethresholds=False
sampleoutputthresholdfile=filetype+'samplethreshold0.5comp.dat'#Output when samplethresholds=True
setthresholdup=0.3
setthresholddown=setthresholdup-0.02

# Reading in the synthetics and making new columns
allsynthetics=ascii.read('allsynthetics.dat',header_start=0,data_start=1)
allsynthetics=allsynthetics[np.where(allsynthetics['synthetic_id']>1412)]
colsynxmid=Column(['-9999.99999999' for x in range(len(allsynthetics))],name='synxmid')
colsynduration=Column(['-0.99999999' for x in range(len(allsynthetics))],name='synduration')
colsynoverlap=Column(['-0.99999999' for x in range(len(allsynthetics))],name='synoverlap')
coluseroverlap=Column(['-0.99999999' for x in range(len(allsynthetics))],name='useroverlap')
coluserxmin=Column(['-9999.99999999' for x in range(len(allsynthetics))],name='userxmin')
coluserxmax=Column(['-9999.99999999' for x in range(len(allsynthetics))],name='userxmax')
coluserxmid=Column(['-0.99999999' for x in range(len(allsynthetics))],name='userxmid')
coluserduration=Column(['-0.99999999' for x in range(len(allsynthetics))],name='userduration')
colmiddiff=Column(['-0.99999999' for x in range(len(allsynthetics))],name='middiff')
colfps=Column(['-1' for x in range(len(allsynthetics))],name='fps')
colthresholddown=Column([setthresholddown for x in range(len(allsynthetics))],name='thresholddown')
colthresholdup=Column([setthresholdup for x in range(len(allsynthetics))],name='thresholdup')
allsynthetics.add_column(colsynxmid)
allsynthetics.add_column(colsynduration)
allsynthetics.add_column(colsynoverlap)
allsynthetics.add_column(coluseroverlap)
allsynthetics.add_column(coluserxmin)
allsynthetics.add_column(coluserxmax)
allsynthetics.add_column(coluserxmid)
allsynthetics.add_column(coluserduration)
allsynthetics.add_column(colmiddiff)
allsynthetics.add_column(colfps)
allsynthetics.add_column(colthresholddown)
allsynthetics.add_column(colthresholdup)

# Overlap function
def funcoverlap(synx, userx):
	'''Returns the total length of the overlapping areas'''
	return max(0, min(synx[1], userx[1]) - max(synx[0], userx[0]))

if samplethresholds==True:
	threshdiff=2
	thresholdup=np.asarray(range(threshdiff+1,101))/100.0
	thresholddown=np.asarray(range(1,101-threshdiff))/100.0
	thresholdup.tolist()
	thresholddown.tolist()
	filewrite=open(sampleoutputthresholdfile,'w')
	filewrite.write('thresholddown&thresholdup&totaltces&syntces&fps&fppercent \n')
else: 
	thresholdup=[setthresholdup]
	thresholddown=[setthresholddown]

for l in range(len(thresholdup)):
	count,i,tces,syntces=0,0,0,0
	with open(inputhistogramfile) as infile:  # Reading the histograms in
		for line in infile:
			linesplit=line.split('&')	
			syntheticid=int(linesplit[0])
			numclassid=float(linesplit[1])
			time=linesplit[2].split(',')
			histcount=linesplit[3].split(',')
			for j in range(len(time)):  # Converting time and histcount to lists of floats
				time[j]=float(time[j])
				histcount[j]=float(histcount[j])
			time,histcount=np.asarray(time),np.asarray(histcount)
			userxmin,userxmax,userxmid=[],[],[]
			for j in range(len(histcount)):
				if j==0:  # Special case for j=0 to avoid [j-1] being out of index
					if histcount[j]/numclassid>=thresholdup[l]:
						userxmin.append(time[j])
				# Determines if histcount crosses upper threshold going up and that we are not currently in a transit
				elif histcount[j]/numclassid>=thresholdup[l] and histcount[j-1]/numclassid<	thresholdup[l] and len(userxmin)==len(userxmax):
					userxmin.append(time[j])
				# Determines if we're crossing the upper threshold going down and are currently in a transit
				if j!=0 and histcount[j]/numclassid<thresholdup[l] and histcount[j-1]/numclassid >= thresholdup[l] and len(userxmin) == (len(userxmax)+1):
					mostrecentcrossing=time[j-1]
				# Determines if histcount crosses lower threshold going down and that we are currently in a transit
				if j!=0 and histcount[j]/numclassid<thresholddown[l] and histcount[j-1]/numclassid >= thresholddown[l] and len(userxmin) == (len(userxmax)+1):
					userxmax.append(mostrecentcrossing)
				# Closes a transit on the last time point if it's above the threshold and we're currently in a transit
				if j==(len(histcount)-1) and histcount[j]/numclassid>=thresholddown[l] and len(	userxmin) == (len(userxmax)+1):
					userxmax.append(time[j])
			for j in range(len(userxmin)):
				userxmid.append(np.median([userxmin[j],userxmax[j]]))
			tces+=len(userxmid)  # Keeping track of total number of TCEs
			# Finds the relevant synthetic transits
			synthetics=allsynthetics[np.where(allsynthetics['synthetic_id']==syntheticid)]
			# Initiating count for total number of detected synthetic transits and saves the index to put that number in each row
			tmpsyntces=0  
			tmpi1=i
			for j in range(len(synthetics)):
				# Finds various values for the synthetic transits
				synx=[synthetics['synxmin'][j],synthetics['synxmax'][j]]
				synduration=synx[1]-synx[0]
				synxmid=np.median(synx)
				allsynthetics['synduration'][i]=synduration
				allsynthetics['synxmid'][i]=synxmid
				synoverlap,useroverlap,userduration=[],[],[]
				for k in range(len(userxmid)):
					# Cycles through each TCE and determines duration and overlap
					userx=[userxmin[k],userxmax[k]]
					userduration.append(userx[1]-userx[0])
					if userduration[-1] < 0.0001: # In case of one point TCEs
						userduration[-1] = 0.01
					overlap=funcoverlap(synx,userx)
					synoverlap.append(overlap/synduration)
					useroverlap.append(overlap/userduration[k])
				if len(synoverlap)>0 and max(synoverlap)>0.5:
					# Only counting overlaps of 50% or more
					index=synoverlap.index(max(synoverlap))
					# Extracting index of overlapping TCE
					allsynthetics['synoverlap'][i]=synoverlap[index]
					allsynthetics['useroverlap'][i]=useroverlap[index]
					allsynthetics['userxmin'][i]=userxmin[index]
					allsynthetics['userxmax'][i]=userxmax[index]
					allsynthetics['userxmid'][i]=userxmid[index]
					allsynthetics['userduration'][i]=userduration[index]
					allsynthetics['middiff'][i]=np.abs(synxmid-userxmid[index])
					syntces+=1
					tmpsyntces+=1
				else: # If overlap < 0.5, fill relevant values with np.nan
					allsynthetics['synoverlap'][i]=allsynthetics['useroverlap'][i]=allsynthetics['userxmin'][i]=allsynthetics['userxmax'][i]=allsynthetics['userxmid'][i]=allsynthetics['userduration'][i]=allsynthetics['middiff'][i]=np.nan
				i+=1
			tmpi2=i
			for j in range(tmpi1,tmpi2): 
				# Determining the number of FPs
				allsynthetics['fps'][j]=max([len(userxmid)-tmpsyntces,0])
			if count > np.inf:
				break
			count+=1
			if samplethresholds != True:
				print '{0:0.2f}'.format(100.0*float(count)/float(2380))+'% completed. Syntheticid = '+str(syntheticid)+'.'
	fps=tces-syntces
	if tces>0:
		fpspercent=float(fps)/float(tces)
	elif tces==0 and fps==0:
		fpspercent=np.nan
	else:
		fpspercent=np.inf
	if samplethresholds==True:
		filewrite.write(str(thresholddown[l])+'&'+str(thresholdup[l])+'&'+str(tces)+'&'+str(syntces)+'&'+str(fps)+'&'+str(fpspercent)+' \n')
		print '{0:0.2f}'.format(100.0*float(l+1)/float(len(thresholdup)))+'% completed.'
	else:
		print 'Total TCEs = '+str(tces)+'.'
		print 'Total synthetics classified as TCEs = '+str(syntces)+'.'
		print 'Total potential false positives = '+str(fps)+' ('+'{0:0.2f}'.format(fpspercent*100.0)+'%).'
		ascii.write(allsynthetics,outputthresholdfile)











