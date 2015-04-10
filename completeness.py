from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plotboxy=True
plotsmooth=False
# Read in the data from the matching algorithm
#allclass=ascii.read('matching.dat',header_start=0,data_start=1,delimiter='&')
allclass=ascii.read('threshold.dat',header_start=0,data_start=1,delimiter=' ')
allclass=allclass[np.where(allclass['kepmag']<15.3765)]
detected=allclass[np.where(allclass['synoverlap']>0.5)]


# Defining bin edges for radius and period
rbins=np.asarray([1.0,2.0,3.0,6.0,12.0])
pbins=np.asarray([0.5,50.0,100.0,200.0,300.0,500.0,700.0,1000.0])

rtext=[1.2,2.2,4.2,8.8]
ptext=[5,55,110,210,360,560,820]
# Finding the 2D histograms for all synthetics and then again for the detected synthetics
histall=np.histogram2d(allclass['period'],allclass['prad'],bins=[pbins,rbins])
histdetected=np.histogram2d(detected['period'],detected['prad'],bins=[pbins,rbins])
completeness=np.transpose(histdetected[0]/histall[0]*100)

#histall=np.transpose(histall)
# Setting up the grid to plot
X, Y = np.meshgrid(pbins, rbins)
histallcount=np.transpose(histall[0])
histdetectedcount=np.transpose(histdetected[0])
if plotboxy==True:
	xlabel='Period (days)'
	ylabel=r'Radius ($R_{\oplus}$)'
	cbarlabel='Percent of synthetics detected above threshold = '+'{0:0.0f}'.format(allclass['thresholdup'][0]*100)+'%'
	# Plotting in boxy form
	fig=plt.figure()
	fig.set_size_inches(6.75,5.0)
	fig.subplots_adjust(right=1.0,top=0.98,left=0.08,bottom=0.095)
	ax=fig.add_subplot(111)
	ax.pcolormesh(X, Y, completeness,cmap='Greens')
	cbar=plt.colorbar(ax.pcolormesh(X, Y, completeness,cmap='Greens'),ticks=range(0,110,10))
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	cbar.set_label(cbarlabel, rotation=270,labelpad=20)
	cbar.set_ticks=[0.5]
	for i in range(len(rtext)):
		if i<2:
			textcolor='k'
		else:
			textcolor='white'
		for j in range(len(ptext)):
			if j<2:
				fontsize=5
			else: 
				fontsize=9
			if j==0:
				y=rtext[i]+0.4
			elif j==1:
				y=rtext[i]-0.15
			else:
				y=rtext[i]
			ax.text(ptext[j],y,' '+'{0:0.0f}'.format(completeness[i][j])+'% \n'+'{0:0.0f}'.format(histdetectedcount[i][j])+'/'+'{0:0.0f}'.format(histallcount[i][j]),fontsize=fontsize,color=textcolor)
	#plt.show()
	ax.text(380,0.4,'Bright sample')
	fig.savefig('completenessbright.eps')

# Plotting in smowoth, interpolated form
if plotsmooth==True:
	fig=plt.figure()
	fig.set_size_inches(6.75,5.0)
	ax=fig.add_subplot(111)
	im = mpl.image.NonUniformImage(ax, interpolation='bilinear',cmap='Greens')
	pcenters = pbins[:-1] + 0.5 * (pbins[1:] - pbins[:-1])
	rcenters = rbins[:-1] + 0.5 * (rbins[1:] - rbins[:-1])
	im.set_data(pcenters, rcenters, completeness)
	ax.images.append(im)
	cbar=plt.colorbar(im)
	ax.set_xlim(pbins[0], pbins[-1])
	ax.set_ylim(rbins[0], rbins[-1])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	cbar.set_label(cbarlabel, rotation=270,labelpad=20)
	plt.show()



plotboxy=True
plotsmooth=False
# Read in the data from the matching algorithm
#allclass=ascii.read('matching.dat',header_start=0,data_start=1,delimiter='&')
allclass=ascii.read('threshold.dat',header_start=0,data_start=1,delimiter=' ')
allclass=allclass[np.where(allclass['kepmag']>15.3765)]
detected=allclass[np.where(allclass['synoverlap']>0.5)]


# Defining bin edges for radius and period
rbins=np.asarray([2.0,3.0,6.0,12.0])
pbins=np.asarray([0.5,50.0,100.0,200.0,300.0,500.0,700.0,1000.0])

rtext=[2.2,4.2,8.8]
ptext=[5,55,110,210,360,560,820]
# Finding the 2D histograms for all synthetics and then again for the detected synthetics
histall=np.histogram2d(allclass['period'],allclass['prad'],bins=[pbins,rbins])
histdetected=np.histogram2d(detected['period'],detected['prad'],bins=[pbins,rbins])
completeness=np.transpose(histdetected[0]/histall[0]*100)

#histall=np.transpose(histall)
# Setting up the grid to plot
X, Y = np.meshgrid(pbins, rbins)
histallcount=np.transpose(histall[0])
histdetectedcount=np.transpose(histdetected[0])
if plotboxy==True:
	xlabel='Period (days)'
	ylabel=r'Radius ($R_{\oplus}$)'
	cbarlabel='Percent of synthetics detected above threshold = '+'{0:0.0f}'.format(allclass['thresholdup'][0]*100)+'%'
	# Plotting in boxy form
	fig=plt.figure()
	fig.set_size_inches(6.75,5.0)
	fig.subplots_adjust(right=1.0,top=0.98,left=0.08,bottom=0.095)
	ax=fig.add_subplot(111)
	ax.pcolormesh(X, Y, completeness,cmap='Greens')
	cbar=plt.colorbar(ax.pcolormesh(X, Y, completeness,cmap='Greens'),ticks=range(0,110,10))
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	cbar.set_label(cbarlabel, rotation=270,labelpad=20)
	cbar.set_ticks=[0.5]
	for i in range(len(rtext)):
		if i<1:
			textcolor='k'
		else:
			textcolor='white'
		for j in range(len(ptext)):
			if j<2:
				fontsize=5
			else: 
				fontsize=9
			if j==0:
				y=rtext[i]+0.4
			elif j==1:
				y=rtext[i]-0.15
			else:
				y=rtext[i]
			ax.text(ptext[j],y,' '+'{0:0.0f}'.format(completeness[i][j])+'% \n'+'{0:0.0f}'.format(histdetectedcount[i][j])+'/'+'{0:0.0f}'.format(histallcount[i][j]),fontsize=fontsize,color=textcolor)
	#plt.show()
	ax.set_ylim(0,12)
	ax.text(392,0.8,'Faint sample')
	fig.savefig('completenessfaint.eps')

# Plotting in smowoth, interpolated form
if plotsmooth==True:
	fig=plt.figure()
	fig.set_size_inches(6.75,5.0)
	ax=fig.add_subplot(111)
	im = mpl.image.NonUniformImage(ax, interpolation='bilinear',cmap='Greens')
	pcenters = pbins[:-1] + 0.5 * (pbins[1:] - pbins[:-1])
	rcenters = rbins[:-1] + 0.5 * (rbins[1:] - rbins[:-1])
	im.set_data(pcenters, rcenters, completeness)
	ax.images.append(im)
	cbar=plt.colorbar(im)
	ax.set_xlim(pbins[0], pbins[-1])
	ax.set_ylim(rbins[0], rbins[-1])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	cbar.set_label(cbarlabel, rotation=270,labelpad=20)
	plt.show()

