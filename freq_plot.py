#!/usr/bin/env python
'''
Ks type frequency histogram

prototype
'''

import sys
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path

# variables
bin_num=30
image_name = "hist.png"


fig = plt.figure(1, (8,8), dpi=300)
ax = axes([0.1,0.1,.8,.8])

fin=file(sys.argv[1]) # data in two-column tab delimited file
data=[[row.strip('\n').split('\t')[0] for row in fin]]
fin.seek(0)
data.append([row.strip('\n').split('\t')[1] for row in fin])

x=[int(i) for i in data[1]]

# calc frequencies
n, bins = np.histogram(x, bin_num)
nsum = float(sum(n))
n=map(lambda x: x/nsum, n)

# get the corners of the rectangles for the histogram
left = np.array(bins[:-1])
right = np.array(bins[1:])
bottom = np.zeros(len(left))
top = bottom + n

# we need a (numrects x numsides x 2) numpy array for the path helper
# function to build a compound path
XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

# get the Path object
barpath = path.Path.make_compound_path_from_polys(XY)

# make a patch out of it
patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
ax.add_patch(patch)

# update the view
ax.set_xlim(left[0], right[-1])
ax.set_ylim(bottom.min(), top.max()+.05)
ax.set_xticks(np.arange(left[0], right[-1],int((-left[0]+right[-1])/bin_num*3)))

ax.set_xlabel('size')
ax.set_ylabel('frequency')
ax.set_title("Size distribution")

ax.grid(True)

##plt.show()
plt.savefig(image_name,dpi=300)
print >>sys.stderr, "print image to %s" % image_name

