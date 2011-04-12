#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# codes from Haibao Tang with some modifications.

"""
%prog gff_file [options] > <output>

convert the gff format into the bed format 
"""

import sys 

try:
    from BCBio.GFF import GFFParser
except:
    print >>sys.stderr, "BCBio.GFF module not found, try download and install from`" \
            "<http://github.com/chapmanb/bcbb/tree/master/gff/BCBio/>"
    sys.exit(1)


def gff_to_bed(gff_file, bed_fh=sys.stdout, cds=True, species=None, rename=False):

    parser = GFFParser()
    seqids = parser.parse(gff_file, None)

    cur_chr = None
    cur_gene_order = 0
    for seqid in seqids:
        for feat in seqid.features:
            subf = feat.sub_features
            if feat.type in ("chromosome", "protein"): continue
            is_cds = any(f.type=="mRNA" or f.type=="CDS" for f in subf) and\
                    feat.type=="gene"
            if cds == is_cds:
                cur_gene_order +=1
                if species != None:
                    seqid_final = species+seqid.id[-2:] # hard coded
                else:
                    seqid_final = seqid.id
                if rename:
                    if seqid.id != cur_chr:
                        cur_gene_order = 1
                        cur_chr = seqid.id
                    gene_name = seqid_final+'g'+'0'*(5-len(str(cur_gene_order)))+str(cur_gene_order)
                else:
                    gene_name = feat.id
                    
                print >>bed_fh, "\t".join(str(x) for x in (seqid_final, int(str(feat.location.start))+1, \
                        feat.location.end, gene_name))  # +1 is hard coded to current BCBio.GFF


if __name__ == "__main__":

    import optparse,sys

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--noncoding", dest="cds", action="store_false", 
            default=True, help="extract coding features?")
    parser.add_option("--species", dest="species", default=None, help="2 letters prefix for species code")
    parser.add_option("--rename", dest="rename", default=False, help="rename genes according to their rank order")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        sys.exit(parser.print_help())

    gff_file = args[0]

    gff_to_bed(gff_file, cds=options.cds, species=options.species, rename=options.rename)

