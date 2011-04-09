#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# block_retriever.py

"""
%prog <.align_file> <output_file> 

retrieve blocks from MCscan output .aligns file

output format:
block_ID qstart_gene qend_gene sstart_gene send_gene block_size score

"""


import os,sys
from optparse import OptionParser
from itertools import groupby


def retrieve_blocks(fi):
    '''function to parse MCSCAN .align file for geenrating block summary.'''
    with open(fi,'r') as f:
        for block_id, block in groupby(f, lambda x: x.split()[2].rstrip(':') if "Alignment" in x else x.split('-')[0].strip()): # hard coded
            block=list(block)
            try: size = block[0].split("N=")[1].split()[0]
            except IndexError: continue
            score = block[0].split("score=")[1].split()[0]
            flag=0 # initialize
            for align in block[1:]:
                align = align.strip('\n').split('\t')
                end1, end2 = align[1:3]
                if (not flag):
                    flag=1
                    start1, start2 = align[1:3]
            if start1 > end1:
                start1,end1=end1,start1
            if start2 > end2:
                start2,end2=end2,start2
            yield (block_id, start1, end1, start2, end2, size, score)
                
                 
if __name__ == "__main__":
    
    usage = "usage: %prog <.align_file> <output_file>"
    parser = OptionParser(usage)
    args=parser.parse_args()[1]
    try:
        fi = args[0]
        fo = args[1]
    except:
        sys.exit(parser.print_help())

    fo = open(fo,'w')
    for block in retrieve_blocks(fi):
        fo.write('\t'.join(block)+'\n')
    fo.close()
    
    
