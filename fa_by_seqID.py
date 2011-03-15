#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog <ID_file> <fasta_file> <biopython>(Y/N)

extract fasta sequences by providing an ID file with each ID as one line.
"""

## fasta_by_seqName.py
## extract fasta according to sequence identifiers stored in a separate file
## one seq a line

def main():
    import sys,os

    if len(sys.argv)==1:
        print "Usage: fasta_by_seqName.py <ID_file> <fasta_file> <biopython>(Y/N)"
        sys.exit()

    file_id = sys.argv[1]
    file_fa = sys.argv[2]

    try:
        fid = open(file_id, 'r')
        ffa = open(file_fa, 'r')
    except IOError:
        print 'cannot open files'
    else:
        print '(else)..file can be opened.'
        ffa.close()

    if sys.argv[3] not in 'YN':
        print 'Please use "Y" or "N" to indicate whether you have biopython package installed.'

    file_out='extracted.fa'
    while os.path.isfile(file_out):
        file_out2=raw_input('file already existed. enter another name for outputfile or enter ooo to overwrite: ')
        if file_out2!='ooo':
            file_out=file_out2
        else:
            break
    outputfile = open(file_out,'w')

    if sys.argv[3]=="N":
        name=fid.readline()
        while name != '':
            name=name.strip('\n')
            ffa = open(file_fa, 'r')
            line=ffa.readline()
            flag=0
            while line != '':
                if '>' in line:
                    if flag==1:
                        outputfile.write('\n')
                        break
                    if name in line:
                        outputfile.write(line)
                        flag=1
                if flag==1 and ('>' not in line):
                    outputfile.write(line.strip('\n'))
                line=ffa.readline()
            else:
                print >>sys.stderr, '%s is not found in the fasta file' % name
            ffa.close()
            name=fid.readline()
            
    else:
        from Bio import SeqIO
        ## load fa
        handle=open(file_fa,'rU')
        fa={}
        for record in SeqIO.parse(handle, "fasta") :
            fa[str(record.id).strip()]=str(record.seq)
        handle.close()    
        print len(fa)
        name=fid.readline().strip('\n')
        while name!= '':
            try:
                outputfile.write(">%s\n%s\n" % (name,fa[name]))
            except KeyError:
                print >>sys.stderr, '%s is not found in the fasta files' % (name)
            name=fid.readline().strip('\n')


    fid.close()
    outputfile.close()

if __name__ == "__main__":
    main()

