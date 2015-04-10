from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter
from matplotlib import rcParams
from astropy.table import Table
import random

# Loading file of all synthetics
matching=ascii.read('matching.dat',header_start=0,data_start=1,delimiter='&')

# Shuffling the list of all transitids (light curves) and keeping only some of them
alltransitids=list(set(matching['transitid']))
random.shuffle(alltransitids)
percent=0.8
transitidspercent=alltransitids[:int(round(len(alltransitids)*0.8))+1]
for i in range(len(transitidspercent)):
	transitidspercent[i]=str(transitidspercent[i])

# Determining which syntheticids match with the transits we're keeping
allsyntheticids=list(set(matching['syntheticid']))
syntheticidspercent=[]
for i in range(len(allsyntheticids)):
	if allsyntheticids[i].split('_')[0] in transitidspercent:
		syntheticidspercent.append(allsyntheticids[i])

print 'Fraction of transitids used for user weighting = '+str(float(len(transitidspercent))/float(len(alltransitids)))+'.'
print 'Fraction of syntheticids used for user weighting = '+str(float(len(syntheticidspercent))/float(len(allsyntheticids)))+'.'

filewrite=open('80randomsynthetics.dat','w')
for i in range(len(syntheticidspercent)):
	filewrite.write(syntheticidspercent[i]+'\n')

filewrite.close()











