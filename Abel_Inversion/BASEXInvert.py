# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 19:06:42 2017

@author: Hillslab
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path
import numpy as np
import matplotlib.pyplot as plt
import abel


# User-defined parameters:
numberOfImages = 125

#Home folder. Should have raw_image and inverted_image subfolders
folder = './'

# Specify the center in y,x (vert,horiz) format
center = (193,204)

#NOTE: To run this you need a folder called "bases" in the working directory

from scipy.misc import toimage

# Specify the path to the file

for i in range(1,numberOfImages+1):
    if i<10:
        filename = os.path.join(folder + 'raw_image/raw_image_00' +  str(i) + '.tif')
    elif i<100:
        filename = os.path.join(folder + 'raw_image/raw_image_0' +  str(i) + '.tif')
    else: 
        filename = os.path.join(folder + 'raw_image/raw_image_' +  str(i) + '.tif')

# Name the output files
#output_image = filename[:-4] + '_Abel_transform.png'
#output_text  = filename[:-4] + '_speeds.txt'
#output_plot  = filename[:-4] + '_comparison.pdf'

# Step 1: Load an image file as a numpy array
    print('Loading ' + filename)
    raw_data = plt.imread(filename).astype('float64')

# Step 2: perform the BASEX transform!
    print('Performing the inverse Abel transform:')

    recon = abel.Transform(raw_data, direction='inverse', method='basex',
                       center=center, transform_options=dict(basis_dir='bases'),
                       verbose=True).transform


    if i<10:
        outname = folder + 'inverted_image/inverted_image_00' +  str(i) + '.tif'
    elif i<100:
        outname = folder + 'inverted_image/inverted_image_0' +  str(i) + '.tif'
    else: 
        outname = folder + 'inverted_image/inverted_image_' +  str(i) + '.tif'
        
    #imsave(outname, recon)
    toimage(recon, cmin=0, cmax=1600).save(outname)

    #np.savetxt(outname,recon,fmt='%.4f',delimiter='\t', newline='\n')
#file = open('180_1_transform.txt','w') 
#file.write(str(recon))
#file.close()

# speeds = abel.tools.vmi.angular_integration(recon)

# Set up some axes
    fig = plt.figure(figsize=(15,4))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)

# Plot the raw data
    im1 = ax1.imshow(raw_data,origin='lower',aspect='auto')
    fig.colorbar(im1,ax=ax1,fraction=.1,shrink=0.9,pad=0.03)
    ax1.set_xlabel('x (pixels)')
    ax1.set_ylabel('y (pixels)')

# Plot the 2D transform
#im2 = ax2.imshow(recon,origin='lower',aspect='auto',clim=(0,2000))
    im2 = ax2.imshow(recon,origin='lower',aspect='auto')
    fig.colorbar(im2,ax=ax2,fraction=.1,shrink=0.9,pad=0.03)
    ax2.set_xlabel('x (pixels)')
    ax2.set_ylabel('y (pixels)')

# Plot the 1D speed distribution

#ax3.plot(*speeds)
#ax3.set_xlabel('Speed (pixel)')
#ax3.set_ylabel('Yield (log)')
#ax3.set_yscale('log')
#ax3.set_ylim(1e2,1e5)

# Prettify the plot a little bit:
    plt.subplots_adjust(left=0.06,bottom=0.17,right=0.95,top=0.89,wspace=0.35,hspace=0.37)

# Show the plots
    plt.show()
