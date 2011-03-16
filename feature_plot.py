#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog <NCBI.tbl> <fasta_file> <output_prefix>

extract features in sliding window and plot.

currently implemented is a method to extract 4 features (genes, pseudogenes, transposable elements, gypsy elements) from NCBI tbl files.
"""

import sys
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import ticker

window = 100000 # adjust here

def count_features(fin):
    counts=[[],[],[],[]]
    genes=[]
    pseudogenes=[]
    TEs=[]
    gypsy=[]
    end=0
    wn=1
    wflag=1
    old_locus=""
    cur_locus=""
    locus_flag=1
    for row in file(fin):
        row=row.strip('\n').split("\t")
        if len(row)==3 and row[2]=="gene": #start of new locus
            old_locus=cur_locus
            locus_flag=0
        if len(row)>3 and locus_flag==0 and "locus_tag"==row[3]:
            cur_locus=row[4].strip('\n')
            locus_flag=1
        try: float(row[1])
        except: pass
        else: end=int(row[1])
        if end > 0+window*wn:
            counts[0].append(len(genes))
            counts[1].append(len(pseudogenes))
            counts[2].append(len(TEs))
            counts[3].append(len(gypsy))
            wflag=0
            genes=[]
            pseudogenes=[]
            TEs=[]
            gypsy=[]
            wn+=1
        else:
            if len(row)==3 and row[2]=="CDS" and (cur_locus not in genes):
                genes.append(cur_locus)
                wflag=1
            if len(row)==5 and row[3]=="pseudo" and (cur_locus not in pseudogenes):
                pseudogenes.append(cur_locus)
                wflag=1
            if "Transposable element gene" in ' '.join(row) and (cur_locus not in TEs):
                TEs.append(cur_locus)
                wflag=1
            if "gypsy" in ' '.join(row) and (cur_locus not in gypsy):
                gypsy.append(cur_locus)
                wflag=1

    if wflag:
        counts[0].append(len(genes))
        counts[1].append(len(pseudogenes))
        counts[2].append(len(TEs))
        counts[3].append(len(gypsy))
    print len(genes), len(pseudogenes), len(TEs), len(gypsy)
    print len(counts)
    length=end
    return length, counts


def count_GC(fin):
    counts=[]
    fin=open(fin,'r')
    fin.readline()
    line = fin.readline()
    abin=""
    wflag=1
    length=0
    print line
    while line:
        length+=len(line)-1
        if len(abin)+len(line)-1 <= window:
            abin+=line.strip("\n")
        else:
            left=(window-len(abin))
            abin+=line.strip("\n")[0:(window-len(abin))]
            gc=abin.count("G")+abin.count("g")+abin.count("C")+abin.count("c")
            counts.append(float(gc)/len(abin))
            wflag=0
            abin = line.strip("\n")[left:]
            if len(abin)>0: wflag=1
        line = fin.readline()
    if wflag:
        gc=abin.count("G")+abin.count("g")+abin.count("C")+abin.count("c")
        counts.append(float(gc)/len(abin)) # normalize the last one
    print length
    return length, counts
    

def draw_feature(ax, feature, length, data, fout, offset):
# suited to plot up to 5 features
    col="rgbymg"
    ax1=fig.add_axes([0.1,offset,.8,.10])
    x = np.arange((0+window/2.0), length+window/2.0-1, window)
    print len(x),len(data), len(x)==len(data)
    ax1.fill_between(x, data, y2=0, color=col[int(offset/.16)-1])
    if feature=="GC":
        ylim(.3,.45) # hard coded
    xlim(0,length)
    xlabel(feature)
    ylabel("counts")
    

if __name__ == '__main__':
        
    usage = "usage: %prog NCBI.tbl fasta output_prefix"
    parser = OptionParser(usage)
    args=parser.parse_args()[1]
    try:
        fin = args[0]
        ffa = args[1]
        fout = args[2]
    except:
        sys.exit(parser.print_help())

    fig = plt.figure(1, (8,8), dpi=300)
    ax = axes([0.1,0.1,.8,.8])
    features=["gene","pseudogene","transposable element","gypsy"] # order is hard coded now, need to change
    length, counts = count_features(fin)
    offset=0
    for feature, data in zip(features, counts):
        offset+=.16
        draw_feature(ax, feature, length, data, fout, offset)

    length, data = count_GC(ffa)
    offset+=.16
    draw_feature(ax, "GC", length, data, fout, offset)
    ax2=fig.add_axes([.1,.95,.8,.05])
    ax2.text(.5,0,"Features plot",size="large", ha="center", zorder=10)
    ax2.set_axis_off()

    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_axis_off()

    #plt.show()
    image_name = fout+".png"
    plt.savefig(image_name,dpi=300)
    print >>sys.stderr, "print image to %s" % image_name
