# Goes through each transitid and finds all classids of that syntheticid
# For each unique classid
	# If no user markings overlap transitid, outputs np.nan (failure to mark transitid) 
	# If user marking(s) does overlap transitid, outputs closest user markings to transitid as measured by midpoint difference (successful marking of transitid)

import numpy as np
from astropy.io import ascii
from astropy.table import Column

### classifications
#"classification_id","created_at","data_location","quarter","start_time","subject_id","user_name","kepler_id","known_transits","synthetic","synthetic_id","kepler_2","xMinRelative","xMaxRelative","xMinGlobal","xMaxGlobal"

def funcoverlap(synx, userx):
	'''Returns the total length of the overlapping areas'''
	return max(0, min(synx[1], userx[1]) - max(synx[0], userx[0]))

# Reading in the synthetic files
synthetics=ascii.read('allsynthetics.dat',header_start=0,data_start=1)
synthetics=synthetics[np.where(synthetics['synthetic_id']>1412)]
filewrite=open('matching.dat','w')
filewrite.write('kepid&fits&i&j&k&l&period&prad&srad&kepmag&activity&syntheticid&transitid&plphase&synpixmin&synpixmax&synxmin&synxmax&synmidpoint&synduration&userxmin&userxmax&usermidpoint&userduration&synoverlap&useroverlap&midpointdiff&username&quarter&classid&createdat&datalocation&starttime&subjectid&knowntransits&syntheticbool&kepler2&xminrelative&xmaxrelative\n')

lline,usernames,kic,datasynthetictimes,syntheticid,xminglobal,xmaxglobal,quarters,classid=[],[],[],[],[],[],[],[],[]
createdat,datalocation,starttime,subjectid,knowntransits,syntheticbool,kepler2,xminrelative,xmaxrelative=[],[],[],[],[],[],[],[],[]
i=0
with open('2015-02-03_planet_hunter_classifications.csv') as infile:
    for line in infile:
    	tmpline=line.replace('\",\"','&')[1:-1]
    	if len(tmpline.split('&')[10]) and i>0 and int(tmpline.split('&')[10])>13:
    		lline.append(tmpline.split('&'))
	    	usernames.append(tmpline.split('&')[6])
	    	classid.append(tmpline.split('&')[0])
	    	quarters.append(tmpline.split('&')[3])
    		kic.append(tmpline.split('&')[7])
    		datasynthetictimes.append(tmpline.split('&')[8])
	    	syntheticid.append(int(tmpline.split('&')[10]))
    		xminglobal.append(tmpline.split('&')[14])
    		xmaxglobal.append(tmpline.split('&')[15][:-1])
    		createdat.append(tmpline.split('&')[1])
    		datalocation.append(tmpline.split('&')[2])
    		starttime.append(tmpline.split('&')[4])
    		subjectid.append(tmpline.split('&')[5])
    		knowntransits.append(tmpline.split('&')[8])
    		syntheticbool.append(tmpline.split('&')[9])
    		kepler2.append(tmpline.split('&')[11])
    		xminrelative.append(tmpline.split('&')[12])
    		xmaxrelative.append(tmpline.split('&')[13])    		
    	if i > np.inf:
    		break
    	i+=1

print 'Completed classification extraction.'

# Setting up a list of data arrays containing all of the information from the classification data file
# More needed data:
usernames,classid,quarters,kic,datasynthetictimes,syntheticid,xminglobal,xmaxglobal=np.asarray(usernames),np.asarray(classid),np.asarray(quarters),np.asarray(kic),np.asarray(datasynthetictimes),np.asarray(syntheticid),np.asarray(xminglobal),np.asarray(xmaxglobal)
# Less needed data
createdat,datalocation,starttime,subjectid,knowntransits,syntheticbool,kepler2,xminrelative,xmaxrelative=np.asarray(createdat),np.asarray(datalocation),np.asarray(starttime),np.asarray(subjectid),np.asarray(knowntransits),np.asarray(syntheticbool),np.asarray(kepler2),np.asarray(xminrelative),np.asarray(xmaxrelative)
# All data
data=[usernames,classid,quarters,kic,datasynthetictimes,syntheticid,xminglobal,xmaxglobal,createdat,datalocation,starttime,subjectid,knowntransits,syntheticbool,kepler2,xminrelative,xmaxrelative]

print 'Completed setting up list of data arrays.'

############################### data index to variable ##################################
#             0              1          2             3                   4            5
#     usernames        classid   quarters          kics  datasynthetictimes  syntheticid    
#             6              7          8             9                  10           11
#    xminglobal     xmaxglobal  createdat  datalocation           starttime    subjectid
#            12             13         14            15                  16
# knowntransits  syntheticbool    kepler2  xminrelative        xmaxrelative

for i in range(len(synthetics)):
#for i in [0]:
	# Extracting the information about the synthetic
	synxmin,synxmax=synthetics['synxmin'][i],synthetics['synxmax'][i]
	synduration=synxmax-synxmin
	synmidpoint=np.average([synxmin,synxmax]) 
	onedata=data[:]  # Copying to an independent list
	# Keeping only the classifications of the current synthetic_id
	# 6 is the last index to prevent it from screwing up the np.where
	for j in [0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,5]:  
		onedata[j]=onedata[j][np.where(onedata[5]==synthetics['synthetic_id'][i])]
	setclassid=list(set(onedata[1]))  # Creating a unique list of classid's
	for j in range(len(setclassid)):  # Iterating through the unique list of classid's
	#for j in [4]:
		classdata=onedata[:]  # Copying to an independent list
		# Keeping only the classifications of the current classid
		# 2 is the last index to prevent it from screwing up the np.where
		for k in [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1]: 
			classdata[k]=classdata[k][np.where(classdata[1]==setclassid[j])]
		numclassifications=len(classdata[1]) # Number of user markings in current classid
		# Extracting all unchanging information
		username=classdata[0][0]
		quarter=' '+classdata[2][0]  # Space added in front so Excel doesn't read in a date
		classid=classdata[1][0]
		createdat=classdata[8][0]
		datalocation=classdata[9][0]
		starttime=classdata[10][0]
		subjectid=classdata[11][0]
		knowntransits=classdata[12][0]
		syntheticbool=classdata[13][0]
		kepler2=classdata[14][0]
		# Setting up overlap and midpoint lists
		alloverlap,alluseroverlap,allmidpointdiff,alluserxmin,alluserxmax,alluserduration,allusermidpoint,alluserxminrel,alluserxmaxrel=[],[],[],[],[],[],[],[],[]
		# Iterating over number of user markkings made on current classid
		for k in range(numclassifications):  
			if len(classdata[6][k])>0:  # True if user made a marking
				# Computing user markings, overlap, and midpointdiff
				tmpuserxmin,tmpuserxmax=float(classdata[6][k]),float(classdata[7][k])
				tmpuserduration=tmpuserxmax-tmpuserxmin
				alluserxmin.append(tmpuserxmin)
				alluserxmax.append(tmpuserxmax)
				tmpuserxminrel=float(classdata[15][k])
				tmpuserxmaxrel=float(classdata[16][k])
				alluserxminrel.append(tmpuserxminrel)
				alluserxmaxrel.append(tmpuserxmaxrel)
				alluserduration.append(tmpuserxmax-tmpuserxmin)
				alloverlap.append(funcoverlap([tmpuserxmin,tmpuserxmax],[synxmin,synxmax])/synduration)  
				alluseroverlap.append(funcoverlap([tmpuserxmin,tmpuserxmax],[synxmin,synxmax])/tmpuserduration)
				tmpusermidpoint=np.average([tmpuserxmin,tmpuserxmax])
				allusermidpoint.append(tmpusermidpoint)
				allmidpointdiff.append(np.abs(synmidpoint-tmpusermidpoint))
			else:  # In case user made no marking
				pass
		# Setting up arrays to use for np.where
		alloverlaparray,allmidpointdiffarray=np.asarray(alloverlap),np.asarray(allmidpointdiff)
		allmidpointwithoverlap=allmidpointdiffarray[np.where(alloverlaparray>0)]
		if len(allmidpointwithoverlap)>0:  # True if there is a marking with overlap
			# Finding index of the closest user marked midpoint with overlap
			usermarkindex=allmidpointdiff.index(np.min(allmidpointwithoverlap))
			synoverlap=alloverlap[usermarkindex]
			useroverlap=alluseroverlap[usermarkindex]
			midpointdiff=allmidpointdiff[usermarkindex]
			userxmin=alluserxmin[usermarkindex]
			userxmax=alluserxmax[usermarkindex]
			userduration=alluserduration[usermarkindex]
			usermidpoint=allusermidpoint[usermarkindex]
			xminrelative=alluserxminrel[usermarkindex]
			xmaxrelative=alluserxmaxrel[usermarkindex]
		else:  # In case of no marking or no overlap
			synoverlap,useroverlap,midpointdiff,userxmin,userxmax,userduration,usermidpoint,usermarkindex,xminrelative,xmaxrelative=np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
		line=[]
		for value in synthetics[i]:
			line.append(value)
		for value in [synmidpoint,synduration,userxmin,userxmax,usermidpoint,userduration,synoverlap,useroverlap,midpointdiff,username,quarter,classid,createdat,datalocation,starttime,subjectid,knowntransits,syntheticbool,kepler2,xminrelative,xmaxrelative]:
			line.append(value)
		for k in range(len(line)):
			if k < len(line)-1:
				filewrite.write(str(line[k])+'&')
			else:
				filewrite.write(str(line[k])+'\n')
	if i%10==0:
		print '{0:0.2f}'.format(100*float(i)/float(len(synthetics)))+'% completed.'

print 'Completed matching synthetics to user markings.'
filewrite.close()



### Problems
# What's the false positive rate?
	# How often are they marking non-transits?
# How many of them are accidentally marking synthetics?
	# What if they mark everything in the light curve as a transit?
# If most points are marked, then it has to rise above the background level for that light curve to qualify as a detection

# *** Dealt with in the alternative formulation using thresholds ('threshold.py')
















