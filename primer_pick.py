#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog <fasta_file>

Design primers for given fasta nucleotide sequences using Primer3 (http://primer3.sourceforge.net/releases.php)
"""

from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Application import generic_run
from Bio.Emboss.Primer3 import *
import os,sys

if len(sys.argv) != 1:
    print >>sys.stderr, "Usage: %prog template.fasta"
    sys.exit()
    
open_fasta = open(sys.argv[1], "r")
fw = open(sys.argv[1]+".pr3", "w")
fw.write("seq_id\tovergo_start\tforward_seq\tforward_start\tforward_length\tforward_tm\tforward_gc\treverse_seq\treverse_start\treverse_length\treverse_tm\treverse_gc\tinput_seq_length\tPCR_product_length\n")

for rec in SeqIO.parse(open_fasta,"fasta"):
    print rec.id
    fw2 = open('in.fas','w')
    SeqIO.write(rec,'in.fas','fasta')
    fw2.close()
    primer_cl = Primer3Commandline(sequence="in.fas",auto=True)
    primer_cl.outfile = "out.pr3"
    primer_cl.numreturn = 3
#    primer_cl.target = str(overgo1s)+","+str(overgo1e/3+overgo2s*2/3) # can specify here the region that requires inclusion in the product
    primer_cl.osize = 20
    primer_cl.maxsize = 26
    primer_cl.otm = 58
    primer_cl.mintm = 52
    primer_cl.mingc = 35
    primer_cl.maxgc = 75
    primer_cl.psizeopt = 200
    primer_cl.prange = "100-400"

    result, messages, errors = generic_run(primer_cl)
    print result

    try: 
        open_outfile = file("out.pr3", "r")
    except: pass
    else:    
        primer_record = read(open_outfile)
        for primer in primer_record.primers:
            product_len = -primer.forward_start+primer.reverse_start+primer.reverse_length
            fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                     (rec.id,overgos[rec.id][0][0],primer.forward_seq,primer.forward_start,primer.forward_length,primer.forward_tm,primer.forward_gc,primer.reverse_seq,primer.reverse_start,primer.reverse_length,primer.reverse_tm,primer.reverse_gc,len(rec.seq),product_len))
        open_outfile.close()
    try: os.system("rm out.pr3")
    except: pass

fw.close()
open_fasta.close()
