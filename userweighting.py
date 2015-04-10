from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter
from matplotlib import rcParams
from astropy.table import Table

# Setting up the maximum amount a single correct classification can increase your score
# Max increase in score = 1.0/completenesscutoff
completenesscutoff=0.1

# Reading in some percentage of syntheticids so (e.g., 80%) to test how the user weighting does on the rest of the syntheticids (e.g, 20%)
testrandomsynthetics=True
if testrandomsynthetics==True:
	filewritename='80percentuserweighting.dat'
	randomsynthetics=ascii.read('80randomsynthetics.dat',data_start=0,names=['syntheticid'])
else:
	filewritename='userweighting.dat'

# Read in the data that I need to check whether users correctly marked synthetics
matching=ascii.read('matching.dat',header_start=0,data_start=1,delimiter='&')

# Setting up the user list and initializing the user weights
allusers=list(set(matching['username']))
userweights=Table({'username':allusers,'upweight':np.ones(len(allusers)),'downweight':np.ones(len(allusers)),'combined':np.ones(len(allusers)),'normupweight':np.zeros(len(allusers)),'normdownweight':np.zeros(len(allusers)),'normcombined':np.zeros(len(allusers)),'numclasses':np.zeros(len(allusers))},names=['username','upweight','downweight','combined','normupweight','normdownweight','normcombined','numclasses'])


# Starting the upweighting portion
for i in range(len(allusers)):
	# Find all classifications of specific user and which syntheticids they classified
	userclasses=matching[np.where(matching['username']==allusers[i])]
	# If using training data, remove non-training syntheticids 
	if testrandomsynthetics == True:
		tmpsyntheticids=list(set(userclasses['syntheticid']))
		syntheticids=[]
		for j in range(len(tmpsyntheticids)):
			if tmpsyntheticids[j] in randomsynthetics['syntheticid']:
				syntheticids.append(tmpsyntheticids[j])
	else:
		syntheticids=list(set(userclasses['syntheticid']))
	# Go through all classifications of each syntheticid and count the total number detected
	for j in range(len(syntheticids)):
		# Next for loop just in case some user saw the same synthetic light curve
		for k in range(len(userclasses[np.where(userclasses['syntheticid']==syntheticids[j])]['synoverlap'])):
			# If user correctly classified a synthetic, then we increase his weight
			if userclasses[np.where(userclasses['syntheticid']==syntheticids[j])]['synoverlap'][k]>0.5: 
				syntheticidclasses=matching[np.where(matching['syntheticid']==syntheticids[j])]
				numidentified=syntheticidclasses[np.where(syntheticidclasses['synoverlap']>0.5)]
				# Find increase in user weight, maximum increase = (1.0/completenesscutoff)-1
				syntheticidcompleteness=float(len(numidentified))/float(len(syntheticidclasses))
				syntheticidcompleteness=max([syntheticidcompleteness,completenesscutoff])
				userweights[i]['upweight']+=1.0/syntheticidcompleteness-1.0
		userweights[i]['numclasses']+=1
	print '{0:0.2f}'.format(100.0*float(i+1)/float(len(allusers)))+'% completed.  '+allusers[i]+' weight = '+str(userweights[i]['upweight'])+'. Numclasses = '+str(int(userweights[i]['numclasses']))+'. i = '+str(i)+'.'


lline,usernames,kic,datasynthetictimes,transitid,xminglobal,xmaxglobal,quarters,classid=[],[],[],[],[],[],[],[],[]
createdat,datalocation,starttime,subjectid,knowntransits,syntheticbool,kepler2,xminrelative,xmaxrelative=[],[],[],[],[],[],[],[],[]
i=0
print 'Loading Planet Hunters database'
# Extracting all classifications of synthetics
with open('2015-02-03_planet_hunter_classifications.csv') as infile:
    for line in infile:
    	tmpline=line.replace('\",\"','&')[1:-1]
    	if len(tmpline.split('&')[10]) and i>0:
    		lline.append(tmpline.split('&'))
	    	usernames.append(tmpline.split('&')[6])
	    	classid.append(tmpline.split('&')[0])
	    	quarters.append(tmpline.split('&')[3])
    		kic.append(tmpline.split('&')[7])
    		datasynthetictimes.append(tmpline.split('&')[8])
	    	transitid.append(int(tmpline.split('&')[10]))
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

print 'Loading of Planet Hunters data completed.'


# Starting the downweighting portion
# Getting a list of unique classid's
setclassid=list(set(classid))
classid,xminglobal,xmaxglobal,usernames=np.asarray(classid),np.asarray(xminglobal),np.asarray(xmaxglobal),np.asarray(usernames)

for i in range(len(usernames)):
	if usernames[i][-1]==' ':
		usernames[i]=usernames[i][:-1]

randomtransitids=[]
for i in range(len(randomsynthetics)):
	randomtransitids.append(randomsynthetics[i][0].split('_')[0])

for i in range(len(setclassid)):
	# Always is true is we're not testing with randoms
	# If we are testing with randoms, only true if the transitid in the list of randoms
	if testrandomsynthetics==False or str(transitid[i]) in randomtransitids:
		# Getting list of markings for each unique classid
		classxmin=xminglobal[np.where(classid==setclassid[i])]
		username=usernames[np.where(classid==setclassid[i])][0]
		countwrong=0
		# If user didn't marking anything, don't do anything
		if len(classxmin)==1 and len(classxmin[0])==0:
			pass
		# If userdid mark something:
		else:
			classxmax=xmaxglobal[np.where(classid==setclassid[i])]
			syntheticmarkings=matching[np.where(matching['classid']==setclassid[i])]
			# For each userxmin, check for a corresponding correctly marked synthetic
			# If it exists, break out of for loop
			# If it doesn't exist, increase countwrong and user's downweight
			for j in range(len(classxmin)):
				isitcorrect=False
				for k in range(len(syntheticmarkings)):
					if np.abs(syntheticmarkings[k]['userxmin']-float(classxmin[j]))<0.000001 and np.abs(syntheticmarkings[k]['userxmax']-float(classxmax[j]))<0.000001 and syntheticmarkings[k]['synoverlap']>0.5:
						isitcorrect=True
						break
				if isitcorrect == False and username in userweights['username']:
					countwrong+=1
					userweights[np.where(userweights['username']==username)[0][0]]['downweight']+=1
		print '{0:0.2f}'.format(100.0*float(i+1)/float(len(setclassid)))+'% completed. '+username+' countwrong = '+str(countwrong)+'. i = '+str(i)+'.'

# Finding the highest weight possible if users marked every single transit and no false positives
syntheticids=list(set(matching['syntheticid']))
highestweight=1
for i in range(len(syntheticids)):
	synthetic=matching[np.where(matching['syntheticid']==syntheticids[i])]
	numidentified=synthetic[np.where(synthetic['synoverlap']>0.5)]
	syntheticidcompleteness=float(len(numidentified))/float(len(synthetic))
	syntheticidcompleteness=max([syntheticidcompleteness,completenesscutoff])
	highestweight+=1.0/syntheticidcompleteness-1.0

print 'Highest weight possible = '+str(highestweight)+'.'

# Normalizing upweights and downweights to 1
userweights['normupweight']=userweights['upweight']/np.average(userweights['upweight'])
userweights['normdownweight']=userweights['downweight']/np.average(userweights['downweight'])
# Combining the upweights and downweights
userweights['combined']=userweights['normupweight']/userweights['normdownweight']
userweights['normcombined']=userweights['combined']/np.average(userweights['combined'])

# Writing output
filewrite=open(filewritename,'w')
for i in range(len(userweights)):
	filewrite.write(userweights[i][0]+' & '+str(userweights[i][1])+' & '+str(userweights[i][2])+' & '+str(userweights[i][3])+' & '+str(userweights[i][4])+' & '+str(userweights[i][5])+' & '+str(userweights[i][6])+' & '+str(userweights[i][7])+' \n')

filewrite.close()





# Counting the number of registered users and not-logged-in users who have classified synthetics and then counting the total number of registered and not-logged-in users
usernames=[]
i=0
print 'Loading Planet Hunters database'
# Extracting all classifications of synthetics
with open('2015-02-03_planet_hunter_classifications.csv') as infile:
    for line in infile:
    	tmpline=line.replace('\",\"','&')[1:-1]
    	if len(tmpline.split('&')[10]) and i>0:
	    	usernames.append(tmpline.split('&')[6])    		
    	i+=1

synsetusernames=list(set(usernames))
usernames=[]
i=0
with open('2015-02-03_planet_hunter_classifications.csv') as infile:
    for line in infile:
    	tmpline=line.replace('\",\"','&')[1:-1]
    	usernames.append(tmpline.split('&')[6])   		
    	i+=1

allsetusernames=list(set(usernames))
countnotloggedinsyns,countregisteredsyns=0,0
for i in range(len(synsetusernames)):
	if 'not-logged-in-' in synsetusernames[i]:
		countnotloggedinsyns+=1
	else:
		countregisteredsyns+=1

countnotloggedinall,countregisteredall=0,0
for i in range(len(allsetusernames)):
	if 'not-logged-in-' in allsetusernames[i]:
		countnotloggedinall+=1
	else:
		countregisteredall+=1

print 'Number of registered users who have classified synthetics = '+str(countregisteredsyns)+'.' # 5525
print 'Number of registered users total = '+str(countregisteredall)+'.' # 7152
print 'Number of registered users who haven\'t classified synthetics = '+str(countregisteredall-countregisteredsyns)+'.' # 1627
print 'Number of not-logged-in users who have classified synthetics = '+str(countnotloggedinsyns)+'.' # 4849
print 'Number of not-logged-in users total = '+str(countnotloggedinall)+'.' # 9437
print 'Number of not-logged-in users who haven\'t classified synthetics = '+str(countnotloggedinall-countnotloggedinsyns)+'.' # 4588


# Counting the number of non-synthetic classifications done by people with >=x number of synthetic classificiations
usernames=[]
i=0
with open('2015-02-03_planet_hunter_classifications.csv') as infile:
    for line in infile:
    	tmpline=line.replace('\",\"','&')[1:-1]
    	if len(tmpline.split('&')[10])==0:
	    	usernames.append(tmpline.split('&')[6])    		
    	i+=1

userweights=ascii.read('userweighting.dat',data_start=0,delimiter='&',names=['username','upweight','downweight','combined','normupweight','normdownweight','normcombined','numclasses'])
num5=np.sum(userweights['numclasses'][np.where(userweights['numclasses']>5)])/np.sum(userweights['numclasses'])
num10=np.sum(userweights['numclasses'][np.where(userweights['numclasses']>10)])/np.sum(userweights['numclasses'])
num5users=userweights['username'][np.where(userweights['numclasses']>=5)]
num10users=userweights['username'][np.where(userweights['numclasses']>=10)]

numclasses5users,numclasses10users=0,0
for i in range(len(usernames)):
	if usernames[i] in num5users:
		numclasses5users+=1
	if usernames[i] in num10users:
		numclasses10users+=1

percent5users=float(numclasses5users)/float(len(usernames))
percent10users=float(numclasses10users)/float(len(usernames))

print 'Percent of classifications of non-synthetics by users with 5+ synthetic classifications = '+str(percent5users) # 88.6%
print 'Percent of classifications by users with 10+ synthetic classifications = '+str(percent10users) # 82.0%





