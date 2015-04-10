### Extracts synthetic information for each individual transitid ###

import numpy as np
from astropy.io import ascii
from astropy.table import Column

# Load the syntheticids and add transitid column
synthetics=ascii.read('Mdwarfsynthetics.csv',header_start=0,data_start=1)
transitcolumn=Column(['0000_00' for x in range(len(synthetics))],name='transitid')
synthetics.add_column(transitcolumn,index=11)
data=ascii.read('Mdwarfsynthetics.csv',header_start=0,data_start=1)
data.add_column(transitcolumn,index=11)

# Creating empty ascii file
for i in range(len(data))[::-1]:
	data.remove_row(i)

# Iterating through each transitid of each syntheticid
for i in range(len(synthetics)):
#for i in range(1):
	allsynxmin=synthetics['synxmin'][i].split(',')
	allsynxmax=synthetics['synxmax'][i].split(',')
	for j in range(len(allsynxmin)):
		synxmin=float(allsynxmin[j])
		synxmax=float(allsynxmax[j])
		line=synthetics[i]
		line['synxmin']=synxmin
		line['synxmax']=synxmax
		transitid=str(synthetics['synthetic_id'][i])+'_'+str(j)
		line['transitid']=transitid
		data.add_row(line)

ascii.write(data,'allsynthetics.dat')

