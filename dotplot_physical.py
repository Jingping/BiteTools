#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""

%prog --qbed <qbed> -s <blastn/lastz_m8_file> -topn <top_n_hits_to_plot>
This script draws dotplot according to physical positions (such as in blastn or lastz results) of subject.
Query does need to be cds seqs with a bed file available.
The blastn/lastz run is a cds against whole genome sequence.

"""

## dotplot_physical.py
## list sort should work for up to a dozen of Mb items in most cases

import sys,os, itertools
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from bed_utils import Bed, BlastLine
from My_utils import multikeySort


def merge_blast(blast_file):
    '''using the query as reference, merge blast hits that overlap with each other
        or are less than a given distance (in nt) apart.'''
    Dm = 100 # merging distance (in nt)
    blasts=[]
    with open(blast_file,'r') as f:
        for line in f:
            line=BlastLine(line)
            blasts.append([line.query, line.subject, line.qstart, line.qstop, line.sstart, line.sstop, line.evalue])
    blasts = multikeySort(blasts, ['0', '1', '4', '5'])
    merged=[]
    cur_hit = []
    for i, b in enumerate(blasts):
        if i==0 or b[0]!= cur_hit[0] or b[1]!= cur_hit[1] :
            cur_hit = b
        elif b[0] == cur_hit[0] and b[1]== cur_hit[1] and b[4]-cur_hit[5] > Dm:
            merged.append(cur_hit)
            cur_hit = b
        else:
            cur_hit = [b[0], b[1], min(cur_hit[2], b[2]), max(cur_hit[3], b[3]), cur_hit[4], b[5], min(cur_hit[6], b[6])]
    return multikeySort(merged, ['0', '6'])


def get_breaks(bed):
    '''get chromosome break positions'''
    simple_bed = [(b.seqid, i) for (i, b) in enumerate(bed.beds)]
    for seqid, ranks in itertools.groupby(simple_bed, key=lambda x:x[0]):
        ranks = list(ranks)
        yield seqid, ranks[0][1], ranks[-1][1]


def dotplot(ax, anchors, qbed, sbed, topn, axes):
    # modified from Haibao Tang's original function
    '''function to plot merged topn anchors.'''
    _ = lambda x: r"$\rm{%s}$" % x.replace(" ", r"\ ")
    get_order = lambda bed: dict((f['accn'], (i, f)) for (i, f) in enumerate(bed))
    get_len = lambda bed: sum([f.end for (i, f) in enumerate(bed)])
    qbed = Bed(qbed); sbed = Bed(sbed)
    xmax, ymax = len(qbed), get_len(sbed)
    print xmax, ymax

    qorder=get_order(qbed)
    chr_len = [f.end for (i, f) in enumerate(sbed)]
    
	# get topn hits
    data = []
    cur_q=""
    for anchor in anchors:
        if anchor[0] != cur_q: n=1
        if n>topn: continue
        try: qgene=qorder[anchor[0].split('.')[0]]
        except: continue    
        if 'r' in anchor[0]: continue
        anchor[4]+=sum(chr_len[0:(int(anchor[1].lstrip('chr0'))-1)])
        if qgene[0] < xmax and anchor[4] < ymax:
            data.append((qgene[0] , anchor[4]))
        cur_q=anchor[0]
        n+=1

    print 'data length: ',len(data)
    
    x, y = zip(*data)
    x, y = np.array(x, 'f')/xmax, np.array(y, 'f')/ymax
    ax.scatter(x, y, c='b', s=.5, lw=0, alpha=.8)
    
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    xchr_labels, ychr_labels = [], []
    cbreaks = {}
    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed):
        if "random" in seqid: continue
        cbreaks[("query", seqid)] = (beg, end)
        xchr_labels.append((seqid, (beg + end)/2))
        x, y= np.array([beg, beg], 'f')/xmax, [0,1]
        ax.plot(x, y, "-", color='g', alpha=.8, zorder=10)
    ax.add_patch(Rectangle((.998,0), .002,1, lw=.2, color='g', fc='g', fill=True, alpha=.8,zorder=10))

    get_breaks_subject = lambda bed: [[f.accn, f.start, f.end] for (i, f) in enumerate(bed)]
    chr_cum=0
    for items in get_breaks_subject(sbed):
        seqid, beg, end = items
        beg+=chr_cum
        end+=chr_cum
        if "random" in seqid: continue
        cbreaks[("subject", seqid)] = (beg, end)
        ychr_labels.append((seqid, (beg + end)/2))
        x, y= [0, 1], np.array([beg, beg], 'f')/ymax
        ax.plot(x, y, "-", color='g', alpha=.8,zorder=10)
        chr_cum = end
    ax.add_patch(Rectangle((0,1), 1, .002, lw=.2, color='g', fc='g', fill=True, alpha=.8,zorder=10))

    
    # plot the chromosome labels
    for label, pos in xchr_labels:
        x, y= pos*1./xmax-.015, 1.02
        if label[0]=="0": label = label.replace("0", "")
        ax.text(x, y, _("%s" % label))
    for label, pos in ychr_labels:
        x, y= -.065, pos*1./ymax-.0065
        if label[0]=="0": label = label.replace("0", "")
        ax.text(x, y, _("%s" % label))
    # plot axis labels
    ax.text(.5, 1.06, _("%s" % axes[0]))
    ax.text(-.1, .5, _("%s" % axes[1]), rotation="vertical")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("--qbed", dest="qbed", help="path to qbed")
    p.add_argument("--sbed", dest="sbed", help="path to sbed")
    p.add_argument("-s", dest="blast_file", help="path to blast file")
    p.add_argument("--topn", dest="topn", type=int, help="top n hits to plot (default=3)", default=3)
    args = p.parse_args()
    qbed, sbed, blast_file, topn = args.qbed, args.sbed, args.blast_file, args.topn
    if len(args._get_kwargs())<3:
        print __doc__
        p.print_help()

    anchors=merge_blast(blast_file)
    print "{0} anchors".format(len(anchors))

    fig = plt.figure(1, (8,8), dpi=300)
    ax = axes([0.1,0.1,.8,.8])

    dotplot(ax, anchors, qbed, sbed, topn, axes=blast_file.split("_")[0:2])

    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_axis_off()

    #plt.show()
    image_name = "".join(blast_file.split(".")[:-1])+"_Top"+str(topn)+".png"
    plt.savefig(image_name,dpi=300)
    print >>sys.stderr, "print image to %s" % image_name


